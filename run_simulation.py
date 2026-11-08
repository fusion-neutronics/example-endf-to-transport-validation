"""Step 2: simulate an OKTAVIAN sphere and compare with the measurement.

    python run_simulation.py                       # the Al sphere, TENDL-2025
    python run_simulation.py --case mn             # the Mn sphere instead
    python run_simulation.py --case mn --neutron-file my-Mn55.endf --label "new Mn55"
    python run_simulation.py --compare endf-b8.1 jeff-4.0

The spheres (--case), each a single-isotope sample with measured neutron and
photon leakage spectra:

    al  Al27  40 cm Type-II vessel, SUS-304
    mn  Mn55  61 cm Type-I vessel, SS-41 soft steel

The Co and Nb spheres are also single-isotope, but are left out until two
problems are understood (issue #1): the Co measurement reports 1.31 neutrons
leaking per source neutron from a 0.5 mean free path shell, more than its
physics allows, and for Nb every library misses by the same amount, which
points at the Type-IV vessel model rather than the data.

Without --neutron-file the evaluation under test is TENDL-2025, converted
from ENDF by convert_to_arrow.py. --neutron-file swaps in any other ENDF-6
evaluation of the sample nuclide, named in the plot by --label (required with
it). It is processed with NJOY here, and the conversion is cached by the
file's sha256, so a changed file is converted again. --compare overlays
reference libraries, fetched by yamc by keyword (tendl-2025, tendl-2017,
endf-b8.1, jeff-4.0, jendl-5.0, fendl-3.2d).

Only the sample's neutron data is under test. The vessel's neutron data
(TENDL-2025) and every element's photon interaction data (ENDF/B-VIII.1;
TENDL publishes no photoatomic data) are loaded by yamc keyword and are the
same for every case, so the differences between the curves come from the
sample's neutron evaluation alone.

Writes results/oktavian_<case>.png and results/oktavian_<case>.json. The JSON
holds, for every library, where its data came from (file path and sha256, or
keyword), the calculated spectra with their Monte Carlo errors, C/E per bin
and the goodness-of-fit numbers, plus the measurement itself.

Everything comes from JAERI-M 94-014 (Maekawa et al., 1994), downloaded at run
time from the FENDL-2 benchmark sublibrary and checked against sha256 pins:

* Geometry: the vessel descriptions and Table 4.1 of n-leak/report.txt. For
  the Al sphere they match the authors' own MCNP deck (Fig. 4.21) surface for
  surface. The Type-I vessel has a 5.1 cm re-entrant hole to the centre whose
  tube thickness the report does not give, so the hole is modelled as empty.
  The sample masses in Table 4.1 agree with these geometries to within 3%.
* Source neutrons: for Al the SI1/SP1 histogram of the MCNP deck, which
  resolves the D-T peak more finely; for Mn its Table 4.2.
* Target gammas: Table 3.1 of g-leak/report.txt, the same for every sphere.
  The tritium target emits 0.0862 gammas per source neutron itself; they are
  transported through the sphere in a separate photon-only run and added.
* Measured leakage: n-leak/<sphere>.dat (per unit lethargy per source
  neutron) and g-leak/<sphere>.dat (per MeV per source neutron).

The files are scanned report pages. The parsers undo the OCR substitutions
they can undo safely (O for 0, l for 1, a stray ~ or comma) and repair lost
exponent signs only where the bin edges must join up. A value that cannot be
read safely is left out of the comparison and reported, never guessed.

One known difference from the experiment: yamc runs to steady state, while
the gamma measurement gated out photons arriving later than about 65 ns, so
slow capture gammas are included here.

Leakage is tallied as track-length flux in a void shell from the vessel out to
R_BIG with a vacuum boundary there. Leaking particles cross it once and never
return, each with path length R_BIG - R to within R/R_BIG, so dividing by
R_BIG gives the leakage current per source neutron.
"""

import argparse
import hashlib
import json
import math
import pathlib
import re
import shutil
import urllib.request

import matplotlib
import yamc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).parent

# ---------------------------------------------------------------------------
# The spheres. To add another, take its density from Table 4.1, its vessel
# type and source table from n-leak/report.txt, and pin its two .dat files.
# ---------------------------------------------------------------------------

SPHERES = {
    "al": {"nuclide": "Al27", "element": "Al", "density": 1.223, "vessel": "Type-II",
           "steel": "SUS-304", "source": "MCNP deck",
           "neutrons": ("n-leak/al40.dat",
                        "3fdefe48a9a45157b7b1109fd387e7b1841460cadc77ca4f13022c56007e4aaf"),
           "photons": ("g-leak/al40.dat",
                       "bfba2b0b0e164f04328165dbcbaa856d23095230cea850ca881541eb3288437f")},
    "mn": {"nuclide": "Mn55", "element": "Mn", "density": 4.37, "vessel": "Type-I",
           "steel": "SS-41",
           "source": "Table 4.2 Neutron source spectrum for LiF, Mn, Cu, Mo and W.",
           "neutrons": ("n-leak/mn60.dat",
                        "51aaf22e1a735248e5e2156b4df025ffd39498fb5a0cc55d83995b80104f5464"),
           "photons": ("g-leak/mn60.dat",
                       "09d2dd28dfd0936a08dce6f15e4ec526222f118d88069be87999a11e0fd293fb")},
}

# Vessel geometry (cm). Type-II is the MCNP deck of Fig. 4.21: a 20 cm central
# void, 0.2 cm walls, and an 11 cm duct (a tube of radius 5.55 to 5.75) along
# +x from x = 8.32. Type-I (60 cm inside, 0.5 cm walls) has no central void,
# only a 5.1 cm hole along +x to the centre, modelled as empty since the
# report gives no tube thickness.
VESSELS = {
    "Type-II": {"void": 10.0, "inner": 10.2, "outer": 19.75, "vessel": 19.95,
                "bore": 5.55, "tube": 5.75, "duct_start": 8.32},
    "Type-I": {"void": None, "inner": None, "outer": 30.0, "vessel": 30.5,
               "bore": 2.55, "tube": 2.55, "duct_start": 0.0},
}

# Vessel steels: mass fractions and density. SUS-304 is card M2 of the MCNP
# deck. SS-41 (JIS general structural steel) is taken as iron.
STEELS = {
    "SUS-304": ({"Cr": 0.185, "Fe": 0.704, "Ni": 0.111}, 7.824),
    "SS-41": ({"Fe": 1.0}, 7.85),
}

STEEL_LIBRARY = "tendl-2025"   # the vessel's neutron data, the same for every case
PHOTON_LIBRARY = "endf-b8.1"   # photon interaction data, the same for every case

# ---------------------------------------------------------------------------

R_BIG = 1.0e5
PARTICLES = 10_000_000

# Display names for the library keywords yamc can fetch.
LIBRARIES = {
    "tendl-2025": "TENDL-2025",
    "tendl-2017": "TENDL-2017",
    "endf-b8.1": "ENDF/B-VIII.1",
    "jeff-4.0": "JEFF-4.0",
    "jendl-5.0": "JENDL-5.0",
    "fendl-3.2d": "FENDL-3.2d",
}

# The measurements (FENDL-2 benchmark sublibrary, NEA mirror), pinned by sha256.
OKTAVIAN = ("https://www.oecd-nea.org/dbforms/data/eva/evatapes/fendl_20/"
            "validation/benchmarks/jaerim94014/oktavian/")
REPORTS = {
    "n-leak/report.txt": "39c07f467e6208532cc96f151160be9df94b6e6ea3fb3c7f68725e867e963b91",
    "g-leak/report.txt": "30969b2ddcf928e7b37ef48487545ff1f27063c0edc512f311b909f863bd24c0",
}
TARGET_PHOTONS_PER_NEUTRON = 0.0862  # Table 3.1, the same target for every sphere


# ---------------------------------------------------------------------------
# Reading the scanned report
# ---------------------------------------------------------------------------

def download(name, sha256):
    path = HERE / "downloads" / "oktavian" / name
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(OKTAVIAN + name, headers={"User-Agent": "Mozilla/5.0"})
        path.write_bytes(urllib.request.urlopen(request).read())
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha256:
        raise SystemExit(f"{path} does not match its sha256; delete it and rerun")
    return path.read_text().replace("\r", "")


def number(token):
    """A number from the scan, with the substitutions that are safe to undo."""
    return float(token.replace("O", "0").replace("l", "1").replace("~", "")
                 .replace(",", ".").replace("Ee", "E+"))


def flips(token):
    """The value as read, and with its exponent sign flipped either way."""
    value = number(token)
    match = re.search(r"E[-+]?(\d+)$", token)
    k = int(match.group(1)) if match else 0
    return [value] if k == 0 else [value, value * 10 ** (2 * k), value / 10 ** (2 * k)]


def read_bins(text, name, need_error=True):
    """The (upper, lower, value, error) rows of a scanned spectrum, descending.

    Edges are repaired only where they must join up: each bin's upper edge
    is the previous bin's lower edge, so a lost exponent sign is put back and
    a missing row is an error, not a silently widened bin. A value (or, with
    need_error, an error) that cannot be read is None, and that bin is left
    out of the comparison. The files mark an empty bin with 1E-38 or 0.
    """
    rows = []
    for line in text.splitlines():
        tokens = line.split()
        if len(tokens) < 3 or not re.match(r"[~\d.]", tokens[0]):
            continue
        try:
            first, second = number(tokens[0]), number(tokens[1])
        except ValueError:
            continue
        # Some files list the lower edge first and some the upper.
        up_tok, lo_tok = (tokens[0], tokens[1]) if first > second else (tokens[1], tokens[0])
        rows.append((up_tok, lo_tok, tokens[2:]))
    if len(rows) > 1 and number(rows[0][0]) < number(rows[1][0]):
        rows.reverse()  # ascending in the file
    out, unreadable, empty = [], [], []
    for up_tok, lo_tok, rest in rows:
        if out:
            upper = out[-1][1]
            printed = min(flips(up_tok), key=lambda c: abs(math.log(c / upper)))
            if not math.isclose(printed, upper, rel_tol=2e-3):
                raise SystemExit(f"{name}: the bin below {upper:.5g} MeV starts at {up_tok}: "
                                 "a row is missing or misread")
        else:
            upper = number(up_tok)
        lower = min(flips(lo_tok), key=lambda c: abs(math.log(c * 1.08 / upper)))
        values = []
        for token in (rest + [None, None])[:2]:
            try:
                values.append(None if token is None else number(token))
            except ValueError:
                values.append(None)
        value, error = values
        where = f"{lower:.4g}-{upper:.4g} MeV"
        if value is None or (need_error and error is None):
            unreadable.append(f"{where} ({' '.join(rest)})")
            value = None
        elif value < 1e-30:
            empty.append(where)
            value = None if need_error else 0.0
        out.append((upper, lower, value, error))
    if unreadable:
        print(f"  {name}: {len(unreadable)} bins unreadable in the scan, left out: "
              + "; ".join(unreadable))
    if empty and need_error:
        print(f"  {name}: {len(empty)} bins have no measurement (1E-38 or 0), left out: "
              + ", ".join(empty))
    return out


def source_spectrum(sphere):
    """Ascending source bin edges (eV) and the fraction of neutrons in each."""
    text = download("n-leak/report.txt", REPORTS["n-leak/report.txt"])
    if sphere["source"] == "MCNP deck":
        # The SI1/SP1 histogram of Fig. 4.21. One SP1 entry, 8.581E-03 at
        # 0.50-0.57 MeV, is ten times its neighbours and makes the
        # probabilities sum to 1.007; read as 8.581E-04 they sum to 0.9993, so
        # that is taken as a digit lost in the scan.
        values = lambda s: [number(t) for t in s.split()]
        edges = values(text.split("SI1")[1].split("C    *** SOURCE")[0])
        probs = values(text.split("SP1")[1].split("C     *****")[0])[1:]
        probs = [8.581e-4 if p == 8.581e-3 else p for p in probs]
        return [e * 1e6 for e in edges], [p / sum(probs) for p in probs]
    # A table of lethargy flux: each bin's share is the flux times ln(up/low).
    body = re.split(r"\nTable 4\.\d|\nFig", text.split(sphere["source"])[-1])[0]
    bins = read_bins(body, sphere["source"], need_error=False)[::-1]
    if any(v is None for _, _, v, _ in bins):
        raise SystemExit(f"{sphere['source']}: a source bin is unreadable in the scan")
    probs = [v * math.log(up / lo) for up, lo, v, _ in bins]
    edges = [bins[0][1]] + [up for up, _, _, _ in bins]
    return [e * 1e6 for e in edges], [p / sum(probs) for p in probs]


def target_photons():
    """Table 3.1: the target's own gammas, 0.1 MeV groups given by upper edge."""
    text = download("g-leak/report.txt", REPORTS["g-leak/report.txt"])
    body = text.split("Table 3.1")[-1].split("Total number")[0]
    table = []
    for line in body.splitlines():
        tokens = line.split()
        if len(tokens) == 2:
            try:
                table.append((number(tokens[0]), number(tokens[1])))
            except ValueError:
                pass
    edges = [(table[0][0] - 0.1) * 1e6] + [upper * 1e6 for upper, _ in table]
    return edges, [p / sum(p for _, p in table) for _, p in table]


def measured(sphere, particle):
    """Ascending edges (eV), values and errors; None marks an unusable bin."""
    name, sha256 = sphere[f"{particle}s"]
    bins = read_bins(download(name, sha256), name)[::-1]
    edges = [bins[0][1] * 1e6] + [up * 1e6 for up, _, _, _ in bins]
    return edges, [b[2] for b in bins], [b[3] for b in bins]


# ---------------------------------------------------------------------------
# Where each case's neutron data comes from
# ---------------------------------------------------------------------------

def sha256_of(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def convert_user_file(endf_path, nuclide):
    """Convert an ENDF-6 neutron evaluation to Arrow, cached by its sha256."""
    endf_path = pathlib.Path(endf_path).expanduser().resolve()
    if not endf_path.is_file():
        raise SystemExit(f"--neutron-file {endf_path} does not exist")
    out_dir = HERE / "data" / "user" / sha256_of(endf_path)[:16]
    existing = sorted(out_dir.glob("*.arrow"))
    if not existing:
        print(f"converting {endf_path} (NJOY, 294 K)")
        written = yamc.convert_neutron_transport(
            str(endf_path), str(out_dir), njoy_exec=shutil.which("njoy") or "njoy",
            temperatures=[294.0], library=f"user:{endf_path.name}")
        existing = [pathlib.Path(written)]
    found = existing[0].name.removesuffix(".arrow")
    if found != nuclide:
        raise SystemExit(f"{endf_path.name} is an evaluation of {found}, but this sphere "
                         f"is {nuclide}. Pass an evaluation of {nuclide}, or choose the "
                         "sphere with --case.")
    return existing[0]


def cases_from_args(args, sphere):
    """The list of (label, source description, sample neutron data)."""
    nuclide = sphere["nuclide"]
    cases = []
    if args.neutron_file:
        arrow = convert_user_file(args.neutron_file, nuclide)
        path = pathlib.Path(args.neutron_file).expanduser().resolve()
        cases.append((args.label, {"file": str(path), "sha256": sha256_of(path)},
                      {nuclide: str(arrow)}))
    else:
        arrow = HERE / "data" / "tendl-2025" / "neutron" / f"{nuclide}.arrow"
        if not arrow.exists():
            raise SystemExit(f"no converted data for {nuclide}: run convert_to_arrow.py "
                             "first, or pass --neutron-file")
        cases.append(("TENDL-2025", {"library": "tendl-2025", "arrow": str(arrow),
                                     "converted_by": "convert_to_arrow.py"},
                      {nuclide: str(arrow)}))
    for keyword in args.compare:
        cases.append((LIBRARIES.get(keyword, keyword), {"library": keyword},
                      {nuclide: keyword}))
    return cases


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

def build_model(sphere, neutron_data):
    """The sample in its vessel. Returns the geometry and the leakage cell."""
    v = VESSELS[sphere["vessel"]]
    sample = yamc.Material({sphere["nuclide"]: 1.0}, sphere["density"], temperature=294,
                           name=sphere["nuclide"])
    sample.read_nuclear_data(neutron_data, photon_data={sphere["element"]: PHOTON_LIBRARY})
    fractions, density = STEELS[sphere["steel"]]
    steel = yamc.Material(fractions, density, fraction_type="mass", temperature=294,
                          name=sphere["steel"])
    natural = yamc.data.element_nuclides()
    steel.read_nuclear_data({n: STEEL_LIBRARY for e in fractions for n in natural[e]},
                            photon_data={e: PHOTON_LIBRARY for e in fractions})

    outer = yamc.Sphere(radius=v["outer"])
    vessel = yamc.Sphere(radius=v["vessel"])
    world = yamc.Sphere(radius=R_BIG, boundary="vacuum")
    bore = yamc.Cylinder("x", y0=0.0, z0=0.0, radius=v["bore"])
    tube = yamc.Cylinder("x", y0=0.0, z0=0.0, radius=v["tube"])
    plane = yamc.Plane("x", offset=v["duct_start"])
    behind, beside = plane.below, plane.above
    duct = beside & tube.below  # the duct and its tube, out through the vessel

    cavity = beside & bore.below & vessel.below
    steel_parts = vessel.below & outer.above & ~duct
    if v["tube"] > v["bore"]:
        steel_parts = steel_parts | (beside & tube.below & bore.above & vessel.below)
    if v["void"]:
        void, inner = yamc.Sphere(radius=v["void"]), yamc.Sphere(radius=v["inner"])
        cavity = cavity | (void.below & behind)
        steel_parts = steel_parts | (inner.below & void.above & behind)
        sample_region = (outer.below & inner.above & behind) | (outer.below & tube.above & beside)
    else:
        sample_region = outer.below & ~duct
    shell = yamc.Cell(vessel.above & world.below, name="leakage shell")
    cells = [
        yamc.Cell(cavity, name="cavity and duct"),
        yamc.Cell(steel_parts, material=steel, name="vessel"),
        yamc.Cell(sample_region, material=sample, name="sample"),
        shell,
    ]
    return yamc.Geometry(cells), shell


def per_bin(result, widths, radius, scale=1.0):
    """Tally mean and standard deviation as leakage per bin width."""
    path = R_BIG - radius / 2
    mean = [scale * m / path / w for m, w in zip(result.mean, widths)]
    std = [scale * d / path / w for d, w in zip(result.standard_deviation, widths)]
    return mean, std


def goodness_of_fit(calc, calc_err, meas, meas_err, widths):
    """Measures of how well a calculation matches the measurement.

    * chi2/N: the mean over bins of ((C - E) / sigma)^2, with sigma the
      experimental and Monte Carlo uncertainties added in quadrature. 1 means
      agreement within the quoted errors. The neutron errors are counting
      statistics only, so chi2/N is large even where C/E looks close to 1; it
      is most useful for ranking libraries on the same model.
    * Integral C/E: total calculated over total measured leakage, over the
      bins that have a measurement.
    * Bin deviation: the root mean square of ln(C/E) over the bins, as a
      percentage, so a bin 20% high and one 20% low count the same.
    * MC error: median and largest per-bin Monte Carlo relative error.
    Bins without a usable measurement, or with no calculated leakage, are skipped.
    """
    pairs = [(c, ce, e, ee) for c, ce, e, ee in zip(calc, calc_err, meas, meas_err)
             if e is not None and ee is not None and c > 0 and e > 0]
    rel = sorted(ce / c for c, ce, _, _ in pairs)
    measured_bins = [(c, e, w) for c, e, w in zip(calc, meas, widths) if e is not None]
    return {
        "chi2_per_bin": sum((c - e) ** 2 / (ce ** 2 + ee ** 2)
                            for c, ce, e, ee in pairs) / len(pairs),
        "integral_c_over_e": (sum(c * w for c, _, w in measured_bins)
                              / sum(e * w for _, e, w in measured_bins)),
        "bin_deviation_percent": 100 * math.sqrt(
            sum(math.log(c / e) ** 2 for c, _, e, _ in pairs) / len(pairs)),
        "mc_error_median_percent": 100 * rel[len(rel) // 2],
        "mc_error_max_percent": 100 * rel[-1],
        "bins": len(pairs),
    }


def plot(case, experiment, results, out):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1]})
    for col, particle in enumerate(("neutron", "photon")):
        exp = experiment[particle]
        mev = [e / 1e6 for e in exp["edges_ev"]]
        mid = [(a + b) / 2 for a, b in zip(mev, mev[1:])]
        top, bottom = axes[0][col], axes[1][col]
        nan = lambda xs: [math.nan if x is None else x for x in xs]
        top.errorbar(mid, nan(exp["value"]), yerr=nan(exp["error"]), fmt="o", ms=3,
                     color="black", label="OKTAVIAN measurement", zorder=5)
        lines = []
        for i, result in enumerate(results):
            calc, err = result[particle]["value"], result[particle]["error"]
            # The first case is the evaluation under test: drawn thick, red and
            # on top. The reference libraries are thinner and fainter behind it.
            under_test = i == 0
            colour = "C3" if under_test else ["C0", "C2", "C4", "C1", "C5", "C6"][(i - 1) % 6]
            style = ({"lw": 2.0, "alpha": 1.0, "zorder": 4} if under_test
                     else {"lw": 1.0, "alpha": 0.7, "zorder": 3})
            top.stairs(calc, mev, color=colour, label=result["label"], **style)
            # The Monte Carlo statistical error, +/- 1 sigma. If it is visibly
            # wide anywhere, run more histories (--particles).
            top.stairs([c + e for c, e in zip(calc, err)], mev,
                       baseline=[max(c - e, 0) for c, e in zip(calc, err)],
                       fill=True, color=colour, alpha=0.25, lw=0, zorder=style["zorder"])
            bottom.stairs(nan(result[particle]["c_over_e"]), mev, color=colour, **style)
            fit = result[particle]["fit"]
            lines.append(f"{result['label'][:28]:28s} {fit['chi2_per_bin']:9.3g} "
                         f"{fit['integral_c_over_e']:6.3f} {fit['bin_deviation_percent']:5.0f}%")
        units = "per unit lethargy" if particle == "neutron" else "per MeV"
        top.set(yscale="log", title=f"{particle.capitalize()} leakage, OKTAVIAN {case} sphere",
                ylabel=f"leakage {units}\nper source neutron")
        # Neutrons: top left, clear of the 14 MeV peak. Photons: top right.
        top.legend(fontsize=8, loc="upper left" if particle == "neutron" else "upper right")
        header = f"{'':28s} {'chi2/N':>9s} {'C/E':>6s} {'dev':>6s}"
        top.text(0.02, 0.03, "\n".join([header] + lines), transform=top.transAxes,
                 fontsize=7.5, family="monospace", va="bottom",
                 bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "grey"})
        bottom.axhline(1, color="grey", lw=0.8)
        bottom.set(ylim=(0, 2), ylabel="C/E", xlabel=f"{particle} energy (MeV)")
        if particle == "neutron":
            bottom.set_xscale("log")
    fig.text(0.01, 0.005, "Shaded bands: Monte Carlo 1 sigma.  C/E: integral of the "
             "leakage.  dev: rms of ln(C/E) over the bins.", fontsize=8, color="grey")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out, dpi=150)


def main():
    parser = argparse.ArgumentParser(
        description="Simulate an OKTAVIAN sphere and compare with the measurement.")
    parser.add_argument("--case", choices=sorted(SPHERES), default="al",
                        help="the sphere to simulate (default: al)")
    parser.add_argument("--neutron-file", type=pathlib.Path,
                        help="an ENDF-6 neutron evaluation of the sample nuclide to test; "
                             "default: TENDL-2025, converted by convert_to_arrow.py")
    parser.add_argument("--label", help="name of the --neutron-file evaluation in the "
                                        "legend, table and JSON (required with it)")
    parser.add_argument("--compare", nargs="*", default=[], metavar="LIBRARY",
                        help="reference libraries to overlay, by yamc keyword: "
                             + ", ".join(LIBRARIES))
    parser.add_argument("--particles", type=int, default=PARTICLES,
                        help=f"source neutrons per library (default: {PARTICLES:.0e})")
    parser.add_argument("--output-dir", type=pathlib.Path, default=HERE / "results")
    args = parser.parse_args()
    if args.neutron_file and not args.label:
        parser.error("--neutron-file needs --label, the name to show for it")

    sphere = SPHERES[args.case]
    radius = VESSELS[sphere["vessel"]]["vessel"]
    cases = cases_from_args(args, sphere)
    source_edges, source_probs = source_spectrum(sphere)
    n_edges, n_meas, n_err = measured(sphere, "neutron")
    g_edges, g_meas, g_err = measured(sphere, "photon")
    n_width = [math.log(hi / lo) for lo, hi in zip(n_edges, n_edges[1:])]  # lethargy
    g_width = [(hi - lo) / 1e6 for lo, hi in zip(g_edges, g_edges[1:])]  # MeV

    # The target's own gammas never touch the neutron data, so one run serves
    # every library.
    geometry, shell = build_model(sphere, cases[0][2])
    target_edges, target_probs = target_photons()
    target_tally = yamc.Tally(scores=["flux"], cells=shell, particle="photon",
                              energy_bins=g_edges)
    target_run = yamc.Model(
        geometry=geometry, tallies=[target_tally],
        source=yamc.PhotonSource(position=(0, 0, 0),
                                 energy=yamc.sources.Histogram(target_edges, target_probs)),
    ).simulate_transport(total_particles=max(args.particles // 5, 1000))
    g_target, g_target_err = per_bin(target_run[target_tally], g_width, radius,
                                     TARGET_PHOTONS_PER_NEUTRON)

    results = []
    for label, source, neutron_data in cases:
        print(f"\n=== {label}")
        geometry, shell = build_model(sphere, neutron_data)
        neutron_tally = yamc.Tally(scores=["flux"], cells=shell, particle="neutron",
                                   energy_bins=n_edges)
        photon_tally = yamc.Tally(scores=["flux"], cells=shell, particle="photon",
                                  energy_bins=g_edges)
        run = yamc.Model(
            geometry=geometry, tallies=[neutron_tally, photon_tally],
            source=yamc.NeutronSource(position=(0, 0, 0),
                                      energy=yamc.sources.Histogram(source_edges, source_probs)),
            transport_secondary_photons=True,
        ).simulate_transport(total_particles=args.particles)
        n_calc, n_calc_err = per_bin(run[neutron_tally], n_width, radius)
        g_neut, g_neut_err = per_bin(run[photon_tally], g_width, radius)
        g_calc = [a + b for a, b in zip(g_neut, g_target)]
        g_calc_err = [math.hypot(a, b) for a, b in zip(g_neut_err, g_target_err)]
        result = {"label": label, "source": source}
        for particle, calc, err, meas, meas_err, widths in (
                ("neutron", n_calc, n_calc_err, n_meas, n_err, n_width),
                ("photon", g_calc, g_calc_err, g_meas, g_err, g_width)):
            result[particle] = {
                "value": calc, "error": err,
                "c_over_e": [c / m if m else None for c, m in zip(calc, meas)],
                "fit": goodness_of_fit(calc, err, meas, meas_err, widths),
            }
        results.append(result)

    experiment = {
        "neutron": {"edges_ev": n_edges, "value": n_meas, "error": n_err,
                    "units": "leakage per unit lethargy per source neutron",
                    "file": sphere["neutrons"][0]},
        "photon": {"edges_ev": g_edges, "value": g_meas, "error": g_err,
                   "units": "leakage per MeV per source neutron",
                   "file": sphere["photons"][0]},
    }
    element = sphere["element"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    png = args.output_dir / f"oktavian_{args.case}.png"
    plot(element, experiment, results, png)
    fractions, density = STEELS[sphere["steel"]]
    report = {
        "benchmark": f"OKTAVIAN {element} sphere (JAERI-M 94-014)",
        "case": args.case,
        "sphere": sphere | {"geometry_cm": VESSELS[sphere["vessel"]]},
        "yamc_version": yamc.__version__,
        "source_neutrons_per_library": args.particles,
        "photon_data": f"{PHOTON_LIBRARY} by yamc keyword, the same for every library",
        "vessel_steel": {"name": sphere["steel"], "composition_by_mass": fractions,
                         "density": density, "neutron_data": STEEL_LIBRARY,
                         "note": "the same for every library"},
        "target_photons": {"value": g_target, "error": g_target_err,
                           "note": "included in every library's photon spectrum"},
        "experiment": experiment,
        "cases": results,
    }
    out_json = args.output_dir / f"oktavian_{args.case}.json"
    out_json.write_text(json.dumps(report, indent=1))

    print(f"\n{'':40s} {'chi2/N':>9s} {'C/E':>7s} {'dev':>6s} {'max MC err':>10s}")
    for result in results:
        for particle in ("neutron", "photon"):
            fit = result[particle]["fit"]
            print(f"{result['label'][:30]:30s} {particle:9s} {fit['chi2_per_bin']:9.3g} "
                  f"{fit['integral_c_over_e']:7.3f} {fit['bin_deviation_percent']:5.0f}% "
                  f"{fit['mc_error_max_percent']:9.1f}%")
            if fit["mc_error_max_percent"] > 5:
                print(f"  some {particle} bins have a Monte Carlo error above 5%: "
                      "raise --particles for a sharper comparison")
    print(f"\nwrote {png}\nwrote {out_json}")


if __name__ == "__main__":
    main()
