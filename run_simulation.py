"""Step 2: simulate the OKTAVIAN Al sphere and compare with the measurement.

A 14 MeV D-T source (the measured spectrum, from the authors' MCNP deck) sits
at the centre of a hollow aluminium sphere (inner radius 10.2 cm, outer 19.75
cm, apparent density 1.223 g/cm3, with the beam duct along +x), exactly as in
the authors' MCNP deck. Neutrons and the photons they produce are transported,
and the neutron and photon spectra leaking out of the sphere are compared with
the measured ones. The gamma rays the tritium target emits itself (0.0862 per
source neutron, Table 3.1) also reach the detector, so they are transported
through the sphere in a second, photon-only run and added.

Simplifications, to keep this to one isotope (Al27):
* The thin SUS-304 vessel walls (0.2 cm) are left out (treated as void).
  They add Cr, Fe and Ni and about 0.05 photons per source neutron to the
  leakage (roughly 13% of the measured photon leakage), so the calculated
  photon spectrum here is expected to sit low.
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
MEMBERS in convert_to_arrow.py. Issue #1 lists every OKTAVIAN sphere with its
vessel, density and source table.

Leakage is tallied as track-length flux in a void shell from the sphere out to
R_BIG with a vacuum boundary there. Leaking particles cross it once and never
return, each with path length R_BIG - R to within R/R_BIG, so dividing by
R_BIG gives the leakage current per source neutron.
"""

import hashlib
import math
import pathlib
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


# Material: neutron data from TENDL-2025 per nuclide, photon data from
# ENDF/B-VIII.1 per element, both as converted by convert_to_arrow.py.
sample = yamc.Material(COMPOSITION, DENSITY, temperature=294, name=LABEL)
sample.read_nuclear_data(
    {n: str(HERE / f"data/tendl-2025/neutron/{n}.arrow") for n in COMPOSITION},
    photon_data={e: str(HERE / f"data/endf-b8.1/photon/{e}.arrow") for e in ELEMENTS},
)

# Geometry: the surfaces of the authors' MCNP deck (Fig. 4.21), vessel walls
# left out.
inner = yamc.Sphere(radius=INNER_RADIUS)
outer = yamc.Sphere(radius=OUTER_RADIUS)
world = yamc.Sphere(radius=R_BIG, boundary="vacuum")
duct = yamc.Cylinder("x", y0=0.0, z0=0.0, radius=DUCT_RADIUS)
duct_start = yamc.Plane("x", offset=DUCT_START)
in_duct = duct.below & duct_start.above

sphere = yamc.Cell(inner.above & outer.below & ~in_duct, material=sample, name=LABEL)
cavity = yamc.Cell(inner.below | (in_duct & outer.below), name="cavity and duct")
shell = yamc.Cell(outer.above & world.below, name="leakage shell")
geometry = yamc.Geometry([sphere, cavity, shell])

# The source spectrum, and the measured spectra whose bins the tallies reuse,
# so calculation and measurement are compared bin for bin.
source_edges, source_probs = source_spectrum()
n_edges, n_meas, n_err = measured_neutrons()
g_edges, g_meas, g_err = measured_photons()

neutron_tally = yamc.Tally(scores=["flux"], cells=shell, particle="neutron",
                           energy_bins=n_edges)
photon_tally = yamc.Tally(scores=["flux"], cells=shell, particle="photon",
                          energy_bins=g_edges)
model = yamc.Model(
    geometry=geometry,
    tallies=[neutron_tally, photon_tally],
    source=yamc.NeutronSource(position=(0, 0, 0),
                              energy=yamc.sources.Histogram(source_edges, source_probs)),
    transport_secondary_photons=True,
)
results = model.simulate_transport(total_particles=PARTICLES)

target_edges, target_probs = target_photons()
target_tally = yamc.Tally(scores=["flux"], cells=shell, particle="photon",
                          energy_bins=g_edges)
target_model = yamc.Model(
    geometry=geometry,
    tallies=[target_tally],
    source=yamc.PhotonSource(position=(0, 0, 0),
                             energy=yamc.sources.Histogram(target_edges, target_probs)),
)
target_results = target_model.simulate_transport(total_particles=PARTICLES // 5)

path = R_BIG - OUTER_RADIUS / 2
n_width = [math.log(hi / lo) for lo, hi in zip(n_edges, n_edges[1:])]  # lethargy
g_width = [(hi - lo) / 1e6 for lo, hi in zip(g_edges, g_edges[1:])]  # MeV


def per_bin(result, widths, scale=1.0):
    """Tally mean and standard deviation as leakage per bin width."""
    mean = [scale * m / path / w for m, w in zip(result.mean, widths)]
    std = [scale * d / path / w for d, w in zip(result.standard_deviation, widths)]
    return mean, std


n_calc, n_calc_err = per_bin(results[neutron_tally], n_width)
g_target, g_target_err = per_bin(target_results[target_tally], g_width,
                                 TARGET_PHOTONS_PER_NEUTRON)
g_neut, g_neut_err = per_bin(results[photon_tally], g_width)
g_calc = [a + b for a, b in zip(g_neut, g_target)]
g_calc_err = [math.hypot(a, b) for a, b in zip(g_neut_err, g_target_err)]


def goodness_of_fit(calc, calc_err, meas, meas_err, widths):
    """Three measures of how well the calculation matches the measurement.

    * chi2/N: the mean over bins of ((C - E) / sigma)^2, with sigma the
      experimental and Monte Carlo uncertainties added in quadrature. 1 means
      agreement within the quoted errors. The neutron errors are counting
      statistics only, a fraction of a percent near the peak, so the neutron
      chi2/N is large even where C/E looks close to 1; read it alongside the
      other two.
    * Integral C/E: total calculated over total measured leakage in the range.
    * Bin deviation: the root mean square of ln(C/E) over the bins, as a
      percentage, so a bin 20% high and one 20% low count the same.
    Bins where either value is zero or negative are skipped.
    """
    pairs = [(c, ce, e, ee) for c, ce, e, ee in zip(calc, calc_err, meas, meas_err)
             if c > 0 and e > 0]
    chi2 = sum((c - e) ** 2 / (ce ** 2 + ee ** 2) for c, ce, e, ee in pairs) / len(pairs)
    integral = (sum(c * w for c, w in zip(calc, widths))
                / sum(e * w for e, w in zip(meas, widths)))
    rms = math.sqrt(sum(math.log(c / e) ** 2 for c, _, e, _ in pairs) / len(pairs))
    return {"chi2/N": chi2, "integral C/E": integral, "bin deviation": 100 * rms,
            "bins": len(pairs)}


n_fit = goodness_of_fit(n_calc, n_calc_err, n_meas, n_err, n_width)
g_fit = goodness_of_fit(g_calc, g_calc_err, g_meas, g_err, g_width)

# Plot: spectra on top, calculated / experiment below.
fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex="col",
                         gridspec_kw={"height_ratios": [3, 1]})
def mc_error(calc, calc_err):
    """Median and largest Monte Carlo relative error (1 sigma) over scored bins."""
    rel = sorted(e / c for c, e in zip(calc, calc_err) if c > 0)
    return rel[len(rel) // 2], rel[-1]


for col, (name, edges, calc, calc_err, meas, err, units, log_x, fit) in enumerate([
    ("Neutron", n_edges, n_calc, n_calc_err, n_meas, n_err, "per unit lethargy", True, n_fit),
    ("Photon", g_edges, g_calc, g_calc_err, g_meas, g_err, "per MeV", False, g_fit),
]):
    mev = [e / 1e6 for e in edges]
    mid = [(a + b) / 2 for a, b in zip(mev, mev[1:])]
    top, bottom = axes[0][col], axes[1][col]
    top.errorbar(mid, meas, yerr=err, fmt="o", ms=3, color="black", label="OKTAVIAN measurement")
    top.stairs(calc, mev, color=f"C{col}", lw=1.5, label=f"yamc + TENDL-2025 {', '.join(COMPOSITION)}")
    # The Monte Carlo statistical error, +/- 1 sigma, as a shaded band. If it is
    # visibly wide anywhere, run more histories (PARTICLES at the top).
    top.stairs([c + e for c, e in zip(calc, calc_err)], mev,
               baseline=[max(c - e, 0) for c, e in zip(calc, calc_err)],
               fill=True, color=f"C{col}", alpha=0.3, lw=0, label="Monte Carlo 1 sigma")
    top.set(yscale="log", title=f"{name} leakage, OKTAVIAN {LABEL}",
            ylabel=f"leakage {units}\nper source neutron")
    top.legend()
    median_err, max_err = mc_error(calc, calc_err)
    top.text(0.02, 0.03,
             f"chi2/N = {fit['chi2/N']:.3g} ({fit['bins']} bins)\n"
             f"integral C/E = {fit['integral C/E']:.3f}\n"
             f"bin deviation (rms ln C/E) = {fit['bin deviation']:.0f}%\n"
             f"MC error per bin: median {100 * median_err:.1f}%, max {100 * max_err:.1f}%",
             transform=top.transAxes, fontsize=9, va="bottom",
             bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "grey"})
    bottom.stairs([c / m if m > 0 else math.nan for c, m in zip(calc, meas)], mev,
                  color=f"C{col}")
    bottom.stairs([(c + e) / m if m > 0 else math.nan for c, e, m in zip(calc, calc_err, meas)],
                  mev, baseline=[(c - e) / m if m > 0 else math.nan
                                 for c, e, m in zip(calc, calc_err, meas)],
                  fill=True, color=f"C{col}", alpha=0.3, lw=0)
    bottom.axhline(1, color="grey", lw=0.8)
    bottom.set(ylim=(0, 2), ylabel="C/E", xlabel=f"{name.lower()} energy (MeV)")
    if log_x:
        bottom.set_xscale("log")
fig.tight_layout()
out = HERE / "results" / f"oktavian_{'_'.join(ELEMENTS).lower()}.png"
out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=150)

for name, fit, energies, calc, calc_err in (
        ("neutrons", n_fit, "0.1-20 MeV", n_calc, n_calc_err),
        ("photons", g_fit, "0.5-11 MeV", g_calc, g_calc_err)):
    median_err, max_err = mc_error(calc, calc_err)
    print(f"{name} {energies}: integral C/E {fit['integral C/E']:.3f}, "
          f"chi2/N {fit['chi2/N']:.3g} over {fit['bins']} bins, "
          f"bin deviation {fit['bin deviation']:.0f}%, "
          f"MC error per bin median {100 * median_err:.1f}% max {100 * max_err:.1f}%")
    if max_err > 0.05:
        print(f"  some {name} bins have a Monte Carlo error above 5%: "
              "increase PARTICLES for a sharper comparison")
print(f"wrote {out}")
