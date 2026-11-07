"""Step 2: simulate the OKTAVIAN Al sphere and compare with the measurement.

    python run_simulation.py
    python run_simulation.py --neutron-file my-evaluations/n-Al027.endf --label "TALYS test 3"
    python run_simulation.py --neutron-file my.endf --compare tendl-2025 endf-b8.1 jeff-4.0

Without options the evaluation under test is TENDL-2025 Al27, converted from
ENDF by convert_to_arrow.py. --neutron-file swaps in any other ENDF-6
evaluation of Al27 (a freshly produced one, for example): it is processed
with NJOY here, and the conversion is cached by the file's sha256, so editing
the file and rerunning converts it again. --compare overlays reference
libraries, fetched by yamc by keyword (tendl-2025, tendl-2017, endf-b8.1,
jeff-4.0, jendl-5.0, fendl-3.2d). Every case uses the same model, source and
ENDF/B-VIII.1 photon interaction data, so the differences between the curves
come from the neutron evaluation alone.

Writes results/oktavian_al.png and results/oktavian_al.json. The JSON holds,
for every case, where its data came from (file path and sha256, or library
keyword), the calculated spectra with their Monte Carlo errors, C/E per bin
and the goodness-of-fit numbers, plus the measurement itself, so another
program can read the comparison without parsing the console.

The model: a 14 MeV D-T source (the measured spectrum, from the authors' MCNP
deck) sits at the centre of a hollow aluminium sphere (inner radius 10.2 cm,
outer 19.75 cm, apparent density 1.223 g/cm3, with the beam duct along +x),
exactly as in the authors' MCNP deck. Neutrons and the photons they produce
are transported, and the neutron and photon spectra leaking out of the sphere
are compared with the measured ones. The gamma rays the tritium target emits
itself (0.0862 per source neutron, Table 3.1) also reach the detector, so they
are transported through the sphere in a separate photon-only run and added to
every case's photon spectrum.

Simplifications, to keep this to one isotope (Al27):
* The thin SUS-304 vessel walls (0.2 cm) are left out (treated as void).
  They add Cr, Fe and Ni and about 0.05 photons per source neutron to the
  leakage (roughly 13% of the measured photon leakage), so the calculated
  photon spectrum here is expected to sit low. This affects every library
  equally, so comparisons between libraries are unaffected.
* yamc runs to steady state, while the gamma measurement gated out photons
  arriving later than 69 ns, so slow capture gammas are included here.

Where each spectrum comes from (all JAERI-M 94-014, downloaded at run time):
* Source neutrons: the SI1/SP1 histogram of the authors' MCNP deck for this
  sphere, printed as Fig. 4.21 in n-leak/report.txt.
* Target gammas: Table 3.1 of g-leak/report.txt (the same target for every
  sphere).
* Measured neutron leakage: n-leak/al40.dat, per unit lethargy per source
  neutron.
* Measured photon leakage: g-leak/al40.dat, per MeV per source neutron.

To change the element or material, edit only the block marked "WHAT TO CHANGE
FOR ANOTHER SPHERE" below, and convert the new evaluations first by editing
NEUTRON_FILES and PHOTON_FILES in convert_to_arrow.py. Issue #1 lists every
OKTAVIAN sphere with its vessel, density and source table.

Leakage is tallied as track-length flux in a void shell from the sphere out to
R_BIG with a vacuum boundary there. Leaking particles cross it once and never
return, each with path length R_BIG - R to within R/R_BIG, so dividing by
R_BIG gives the leakage current per source neutron.
"""

import argparse
import hashlib
import json
import math
import pathlib
import shutil
import urllib.request

import matplotlib
import yamc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).parent

# ---------------------------------------------------------------------------
# WHAT TO CHANGE FOR ANOTHER SPHERE
# Everything specific to the aluminium sphere is here. For another OKTAVIAN
# sphere, take the numbers from Table 4.1 and the vessel descriptions in
# n-leak/report.txt (summarised in issue #1).
# ---------------------------------------------------------------------------

# The sample: nuclide name -> atom fraction, and the apparent (powder) density
# from Table 4.1. A natural element needs every isotope, for example
# {"Fe54": 0.05845, "Fe56": 0.91754, ...}, and each one must have been
# converted by convert_to_arrow.py first.
COMPOSITION = {"Al27": 1.0}
DENSITY = 1.223  # g/cm3
ELEMENTS = ["Al"]  # one photon data file per element
LABEL = "Al sphere"
CASE = "al"  # used in the output file names

# The sphere (cm). Al used the 40 cm Type-II vessel: a 20 cm central void,
# the sample from 10.2 to 19.75 cm, and a re-entrant beam duct of radius 5.75
# cm along +x starting at x = 8.32. Other vessel types have other radii.
INNER_RADIUS = 10.2
OUTER_RADIUS = 19.75
DUCT_RADIUS = 5.75
DUCT_START = 8.32

# The measured spectra for this sphere, and their sha256 pins. Find the file
# names in the FENDL-2 directory listing; compute a new pin with sha256sum
# after the first download and check the file by eye.
MEASURED_NEUTRONS = "n-leak/al40.dat"
MEASURED_PHOTONS = "g-leak/al40.dat"
MEASURED_PINS = {
    MEASURED_NEUTRONS: "3fdefe48a9a45157b7b1109fd387e7b1841460cadc77ca4f13022c56007e4aaf",
    MEASURED_PHOTONS: "bfba2b0b0e164f04328165dbcbaa856d23095230cea850ca881541eb3288437f",
}

# The source neutron spectrum. Al is the only sphere with an MCNP deck in the
# report, so the others must use their Table 4.2 to 4.6 entry instead (see
# source_spectrum below): Table 4.3 for Co, 4.2 for Mn, 4.6 for Nb, and so on.

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
FILES = {
    "n-leak/report.txt": "39c07f467e6208532cc96f151160be9df94b6e6ea3fb3c7f68725e867e963b91",
    "g-leak/report.txt": "30969b2ddcf928e7b37ef48487545ff1f27063c0edc512f311b909f863bd24c0",
    **MEASURED_PINS,
}
TARGET_PHOTONS_PER_NEUTRON = 0.0862  # Table 3.1, the same target for every sphere


def download(name):
    path = HERE / "downloads" / "oktavian" / name
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(OKTAVIAN + name, headers={"User-Agent": "Mozilla/5.0"})
        path.write_bytes(urllib.request.urlopen(request).read())
    if hashlib.sha256(path.read_bytes()).hexdigest() != FILES[name]:
        raise SystemExit(f"{path} does not match its sha256; delete it and rerun")
    return path.read_text()


def rows(text, columns=4):
    """Numeric rows of a scanned table, undoing the OCR's O-for-0 and l-for-1."""
    out = []
    for line in text.splitlines():
        tokens = [t.replace("O", "0").replace("l", "1").replace("Ee", "E+")
                  for t in line.split()]
        try:
            out.append([float(t) for t in tokens])
        except ValueError:
            continue
    return [r for r in out if len(r) == columns]


def source_spectrum():
    """The source neutron spectrum from the authors' MCNP deck (Fig. 4.21).

    Returns ascending bin edges (eV) and the fraction of source neutrons in
    each bin. This deck is for the Al sphere only. For another sphere, parse
    its Table 4.2 to 4.6 from n-leak/report.txt instead: columns lower edge,
    upper edge (MeV), lethargy flux and error, so each bin's probability is
    the lethargy flux times ln(upper / lower). Watch for exponent signs lost
    in the scan (1.3574E-01 printed for 1.3574E+01 in Table 4.4).

    The deck's SI1/SP1 histogram resolves the D-T peak in 0.25 MeV bins,
    where Table 4.4 lumps it into 13.6-15.0 and 15.0-16.6 MeV and so smears
    the peak. One SP1 entry, 8.581E-03 at 0.50-0.57 MeV, is ten times its
    neighbours and makes the probabilities sum to 1.007; read as 8.581E-04
    they sum to 0.9993, so that is taken as a digit lost in the scan.
    """
    text = download("n-leak/report.txt")
    numbers = lambda s: [float(t.replace("O", "0").replace("l", "1")) for t in s.split()]
    edges = numbers(text.split("SI1")[1].split("C    *** SOURCE")[0])  # MeV
    probs = numbers(text.split("SP1")[1].split("C     *****")[0])[1:]  # per bin
    probs = [8.581e-4 if p == 8.581e-3 else p for p in probs]
    return [e * 1e6 for e in edges], [p / sum(probs) for p in probs]


def target_photons():
    """Table 3.1: gamma spectrum of the tritium target, 0.1 MeV groups by upper edge.

    The target is the same for every sphere, so this needs no change.
    """
    text = download("g-leak/report.txt").split("Table 3.1")[-1].split("Total number")[0]
    table = rows(text, columns=2)  # upper energy (MeV), photons per source neutron
    edges = [(table[0][0] - 0.1) * 1e6] + [upper * 1e6 for upper, _ in table]
    return edges, [p / sum(p for _, p in table) for _, p in table]


def measured_neutrons():
    """Measured neutron leakage per unit lethargy per source neutron.

    From n-leak/<sphere>.dat: columns upper edge, lower edge (MeV), value and
    the counting error, in descending energy. Returned with ascending edges.
    """
    table = rows(download(MEASURED_NEUTRONS))[::-1]  # upper, lower, value, error
    return [table[0][1] * 1e6] + [r[0] * 1e6 for r in table], \
        [r[2] for r in table], [r[3] for r in table]


def measured_photons():
    """Measured photon leakage per MeV per source neutron.

    From g-leak/<sphere>.dat: columns lower edge, upper edge (MeV), value and
    error, in ascending energy.
    """
    table = rows(download(MEASURED_PHOTONS))  # lower, upper, value, error
    return [table[0][0] * 1e6] + [r[1] * 1e6 for r in table], \
        [r[2] for r in table], [r[3] for r in table]


# ---------------------------------------------------------------------------
# Where each case's neutron data comes from
# ---------------------------------------------------------------------------

def sha256_of(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def convert_user_file(endf_path):
    """Convert an ENDF-6 neutron evaluation to Arrow, cached by its sha256.

    Returns ``(nuclide, arrow_path)``. The cache key is the file's content,
    not its name, so a re-evaluated file with the same name is converted
    again rather than silently reusing the old conversion.
    """
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
    arrow = existing[0]
    nuclide = arrow.name.removesuffix(".arrow")
    if nuclide not in COMPOSITION:
        raise SystemExit(
            f"{endf_path.name} is an evaluation of {nuclide}, but the {LABEL} needs "
            f"{', '.join(COMPOSITION)}. Pass an evaluation of a nuclide in the sample.")
    return nuclide, arrow


def cases_from_args(args):
    """The list of (label, source description, nuclide -> data) to simulate."""
    cases = []
    if args.neutron_file:
        nuclide, arrow = convert_user_file(args.neutron_file)
        path = pathlib.Path(args.neutron_file).expanduser().resolve()
        cases.append((args.label, {"file": str(path), "sha256": sha256_of(path)},
                      {nuclide: str(arrow)}))
    else:
        default = HERE / "data" / "tendl-2025" / "neutron"
        missing = [n for n in COMPOSITION if not (default / f"{n}.arrow").exists()]
        if missing:
            raise SystemExit(f"no converted data for {', '.join(missing)} in {default}: "
                             "run convert_to_arrow.py first, or pass --neutron-file")
        cases.append(("TENDL-2025",
                      {"library": "tendl-2025", "converted_by": "convert_to_arrow.py",
                       "arrow": str(default)},
                      {n: str(default / f"{n}.arrow") for n in COMPOSITION}))
    for keyword in args.compare:
        cases.append((LIBRARIES.get(keyword, keyword), {"library": keyword},
                      {n: keyword for n in COMPOSITION}))
    return cases


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

def build_model(neutron_data):
    """The sphere filled with the sample, loaded from ``neutron_data``.

    ``neutron_data`` maps each nuclide to an Arrow path or a library keyword.
    The photon interaction data is always the ENDF/B-VIII.1 conversion from
    convert_to_arrow.py, so only the neutron data differs between cases.
    """
    sample = yamc.Material(COMPOSITION, DENSITY, temperature=294, name=LABEL)
    sample.read_nuclear_data(
        neutron_data,
        photon_data={e: str(HERE / f"data/endf-b8.1/photon/{e}.arrow") for e in ELEMENTS},
    )
    # The surfaces of the authors' MCNP deck (Fig. 4.21), vessel walls left out.
    inner = yamc.Sphere(radius=INNER_RADIUS)
    outer = yamc.Sphere(radius=OUTER_RADIUS)
    world = yamc.Sphere(radius=R_BIG, boundary="vacuum")
    duct = yamc.Cylinder("x", y0=0.0, z0=0.0, radius=DUCT_RADIUS)
    duct_start = yamc.Plane("x", offset=DUCT_START)
    in_duct = duct.below & duct_start.above
    sphere = yamc.Cell(inner.above & outer.below & ~in_duct, material=sample, name=LABEL)
    cavity = yamc.Cell(inner.below | (in_duct & outer.below), name="cavity and duct")
    shell = yamc.Cell(outer.above & world.below, name="leakage shell")
    return yamc.Geometry([sphere, cavity, shell]), shell


def per_bin(result, widths, scale=1.0):
    """Tally mean and standard deviation as leakage per bin width."""
    path = R_BIG - OUTER_RADIUS / 2
    mean = [scale * m / path / w for m, w in zip(result.mean, widths)]
    std = [scale * d / path / w for d, w in zip(result.standard_deviation, widths)]
    return mean, std


def goodness_of_fit(calc, calc_err, meas, meas_err, widths):
    """Measures of how well a calculation matches the measurement.

    * chi2/N: the mean over bins of ((C - E) / sigma)^2, with sigma the
      experimental and Monte Carlo uncertainties added in quadrature. 1 means
      agreement within the quoted errors. The neutron errors are counting
      statistics only, a fraction of a percent near the peak, so the neutron
      chi2/N is large even where C/E looks close to 1; it is most useful for
      ranking libraries on the same model.
    * Integral C/E: total calculated over total measured leakage in the range.
    * Bin deviation: the root mean square of ln(C/E) over the bins, as a
      percentage, so a bin 20% high and one 20% low count the same.
    * MC error: median and largest per-bin Monte Carlo relative error.
    Bins where either value is zero or negative are skipped.
    """
    pairs = [(c, ce, e, ee) for c, ce, e, ee in zip(calc, calc_err, meas, meas_err)
             if c > 0 and e > 0]
    rel = sorted(ce / c for c, ce, _, _ in pairs)
    return {
        "chi2_per_bin": sum((c - e) ** 2 / (ce ** 2 + ee ** 2)
                            for c, ce, e, ee in pairs) / len(pairs),
        "integral_c_over_e": (sum(c * w for c, w in zip(calc, widths))
                              / sum(e * w for e, w in zip(meas, widths))),
        "bin_deviation_percent": 100 * math.sqrt(
            sum(math.log(c / e) ** 2 for c, _, e, _ in pairs) / len(pairs)),
        "mc_error_median_percent": 100 * rel[len(rel) // 2],
        "mc_error_max_percent": 100 * rel[-1],
        "bins": len(pairs),
    }


def plot(experiment, results, out):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1]})
    for col, particle in enumerate(("neutron", "photon")):
        exp = experiment[particle]
        mev = [e / 1e6 for e in exp["edges_ev"]]
        mid = [(a + b) / 2 for a, b in zip(mev, mev[1:])]
        top, bottom = axes[0][col], axes[1][col]
        top.errorbar(mid, exp["value"], yerr=exp["error"], fmt="o", ms=3, color="black",
                     label="OKTAVIAN measurement", zorder=5)
        lines = []
        for i, case in enumerate(results):
            calc, err = case[particle]["value"], case[particle]["error"]
            # The first case is the evaluation under test: drawn thick, red and
            # on top. The reference libraries are thinner and fainter behind it.
            under_test = i == 0
            colour = "C3" if under_test else ["C0", "C2", "C4", "C1", "C5", "C6"][(i - 1) % 6]
            style = ({"lw": 2.0, "alpha": 1.0, "zorder": 4} if under_test
                     else {"lw": 1.0, "alpha": 0.7, "zorder": 3})
            top.stairs(calc, mev, color=colour, label=case["label"], **style)
            # The Monte Carlo statistical error, +/- 1 sigma. If it is visibly
            # wide anywhere, run more histories (--particles).
            top.stairs([c + e for c, e in zip(calc, err)], mev,
                       baseline=[max(c - e, 0) for c, e in zip(calc, err)],
                       fill=True, color=colour, alpha=0.25, lw=0, zorder=style["zorder"])
            ce = case[particle]["c_over_e"]
            bottom.stairs([math.nan if v is None else v for v in ce], mev, color=colour,
                          **style)
            fit = case[particle]["fit"]
            lines.append(f"{case['label'][:28]:28s} {fit['chi2_per_bin']:9.3g} "
                         f"{fit['integral_c_over_e']:6.3f} {fit['bin_deviation_percent']:5.0f}%")
        units = "per unit lethargy" if particle == "neutron" else "per MeV"
        top.set(yscale="log", title=f"{particle.capitalize()} leakage, OKTAVIAN {LABEL}",
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
        description="Simulate the OKTAVIAN Al sphere and compare with the measurement.")
    parser.add_argument("--neutron-file", type=pathlib.Path,
                        help="an ENDF-6 neutron evaluation to test, of a nuclide in the "
                             f"sample ({', '.join(COMPOSITION)}); default: the TENDL-2025 "
                             "file converted by convert_to_arrow.py")
    parser.add_argument("--label", help="name of the --neutron-file evaluation in the "
                                        "legend, table and JSON (required with it)")
    parser.add_argument("--compare", nargs="*", default=[], metavar="LIBRARY",
                        help="reference libraries to overlay, by yamc keyword: "
                             + ", ".join(LIBRARIES))
    parser.add_argument("--particles", type=int, default=PARTICLES,
                        help=f"source neutrons per case (default: {PARTICLES:.0e})")
    parser.add_argument("--output-dir", type=pathlib.Path, default=HERE / "results")
    args = parser.parse_args()
    if args.neutron_file and not args.label:
        parser.error("--neutron-file needs --label, the name to show for it")

    cases = cases_from_args(args)
    source_edges, source_probs = source_spectrum()
    n_edges, n_meas, n_err = measured_neutrons()
    g_edges, g_meas, g_err = measured_photons()
    n_width = [math.log(hi / lo) for lo, hi in zip(n_edges, n_edges[1:])]  # lethargy
    g_width = [(hi - lo) / 1e6 for lo, hi in zip(g_edges, g_edges[1:])]  # MeV

    # The target's own gammas never touch the neutron data, so one run serves
    # every case. It uses the first case's model; only its photon data matters.
    geometry, shell = build_model(cases[0][2])
    target_edges, target_probs = target_photons()
    target_tally = yamc.Tally(scores=["flux"], cells=shell, particle="photon",
                              energy_bins=g_edges)
    target_run = yamc.Model(
        geometry=geometry, tallies=[target_tally],
        source=yamc.PhotonSource(position=(0, 0, 0),
                                 energy=yamc.sources.Histogram(target_edges, target_probs)),
    ).simulate_transport(total_particles=max(args.particles // 5, 1000))
    g_target, g_target_err = per_bin(target_run[target_tally], g_width,
                                     TARGET_PHOTONS_PER_NEUTRON)

    results = []
    for label, source, neutron_data in cases:
        print(f"\n=== {label}")
        geometry, shell = build_model(neutron_data)
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
        n_calc, n_calc_err = per_bin(run[neutron_tally], n_width)
        g_neut, g_neut_err = per_bin(run[photon_tally], g_width)
        g_calc = [a + b for a, b in zip(g_neut, g_target)]
        g_calc_err = [math.hypot(a, b) for a, b in zip(g_neut_err, g_target_err)]
        case = {"label": label, "source": source}
        for particle, calc, err, meas, meas_err, widths in (
                ("neutron", n_calc, n_calc_err, n_meas, n_err, n_width),
                ("photon", g_calc, g_calc_err, g_meas, g_err, g_width)):
            case[particle] = {
                "value": calc, "error": err,
                "c_over_e": [c / m if m > 0 else None for c, m in zip(calc, meas)],
                "fit": goodness_of_fit(calc, err, meas, meas_err, widths),
            }
        results.append(case)

    experiment = {
        "neutron": {"edges_ev": n_edges, "value": n_meas, "error": n_err,
                    "units": "leakage per unit lethargy per source neutron",
                    "file": MEASURED_NEUTRONS},
        "photon": {"edges_ev": g_edges, "value": g_meas, "error": g_err,
                   "units": "leakage per MeV per source neutron",
                   "file": MEASURED_PHOTONS},
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    png = args.output_dir / f"oktavian_{CASE}.png"
    plot(experiment, results, png)
    report = {
        "benchmark": f"OKTAVIAN {LABEL} (JAERI-M 94-014)",
        "yamc_version": yamc.__version__,
        "source_neutrons_per_case": args.particles,
        "photon_data": "endf-b8.1 (converted by convert_to_arrow.py), the same for every case",
        "target_photons": {"value": g_target, "error": g_target_err,
                           "note": "included in every case's photon spectrum"},
        "experiment": experiment,
        "cases": results,
    }
    out_json = args.output_dir / f"oktavian_{CASE}.json"
    out_json.write_text(json.dumps(report, indent=1))

    print(f"\n{'':40s} {'chi2/N':>9s} {'C/E':>7s} {'dev':>6s} {'max MC err':>10s}")
    for case in results:
        for particle in ("neutron", "photon"):
            fit = case[particle]["fit"]
            print(f"{case['label'][:30]:30s} {particle:9s} {fit['chi2_per_bin']:9.3g} "
                  f"{fit['integral_c_over_e']:7.3f} {fit['bin_deviation_percent']:5.0f}% "
                  f"{fit['mc_error_max_percent']:9.1f}%")
            if fit["mc_error_max_percent"] > 5:
                print(f"  some {particle} bins have a Monte Carlo error above 5%: "
                      "raise --particles for a sharper comparison")
    print(f"\nwrote {png}\nwrote {out_json}")


if __name__ == "__main__":
    main()
