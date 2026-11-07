# ENDF to transport to validation

A minimal end-to-end example: one ENDF evaluation goes in, and a comparison
with a measured benchmark experiment comes out.

1. **ENDF in.** The TENDL-2025 evaluation of Al27 (aluminium has a single
   stable isotope). NJOY reconstructs the resonances and Doppler broadens to
   294 K, and `yamc.convert_neutron_transport` writes the Arrow data yamc
   transports with, including TENDL's photon production. The photons are
   transported with the ENDF/B-VIII.1 Al photoatomic and atomic relaxation
   evaluations (`yamc.convert_photon`), because TENDL carries no photoatomic
   data.
2. **Transport.** yamc models the OKTAVIAN aluminium sphere (Osaka
   University): a D-T neutron source at the centre of a 40 cm hollow sphere of
   aluminium powder, built from the experimenters' own MCNP deck. Neutrons and
   the photons they produce are transported, and the spectra leaking out of the
   sphere are tallied.
3. **Validation.** Both leakage spectra are plotted against the measurement,
   with calculated over experimental (C/E) ratios.

![OKTAVIAN Al comparison](results/oktavian_al.png)

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Everything is a prebuilt wheel: `yamc` from PyPI, and `njoy2016` from the
wheel index named in `requirements.txt`.

## Run

```bash
python convert_to_arrow.py   # ENDF -> Arrow, once
python run_simulation.py     # transport, comparison, results/oktavian_al.png
```

The transport takes about a second for a million source neutrons.

## Where the data comes from

Nothing third party is stored in this repository. Each input is downloaded on
your machine into `downloads/` and checked against a pinned sha256.

* **TENDL-2025** neutron sublibrary, `tendl.imperial.ac.uk`. TENDL publishes it
  as one archive of about 3.5 GB, so the first run downloads all of it and
  extracts only `n-Al027.tendl`.
* **ENDF/B-VIII.1** photoatomic and atomic relaxation sublibraries, NNDC.
* **OKTAVIAN measurements**, JAERI-M 94-014 (Maekawa et al., 1994, chapters 3
  and 4), from the FENDL-2 benchmark sublibrary as mirrored by the OECD/NEA.
  The files are scanned report pages. The script undoes the OCR substitutions
  (`O` for `0`, `l` for `1`), and one source spectrum entry that is
  unambiguously a dropped digit is corrected, with the reasoning given in
  `run_simulation.py`.

## Goodness of fit

Each panel of the plot, and the script's output, gives three numbers:

* **chi2/N**: the mean over bins of ((C - E) / sigma)^2, with the experimental
  and Monte Carlo uncertainties added in quadrature. 1 means agreement within
  the quoted errors. The measured neutron errors are counting statistics only
  (a fraction of a percent near the peak), and neither spectrum is folded with
  the detector's energy resolution, so expect values well above 1 even for a
  good evaluation. It is most useful for comparing two libraries on the same
  model.
* **integral C/E**: total calculated over total measured leakage.
* **bin deviation**: the root mean square of ln(C/E) over the bins, as a
  percentage.

## What the comparison includes, and what it leaves out

To keep this to one isotope, the model leaves out the 0.2 cm SUS-304 vessel
walls. Those walls would bring in Cr, Fe and Ni.

| | neutrons 0.1-20 MeV | neutrons 1-10 MeV | photons 0.5-11 MeV |
|---|---|---|---|
| this example (Al27 only) | 1.02 | 0.76 | 0.74 |
| with the steel walls added | 1.02 | 0.87 | 0.85 |

(C/E of the integrated leakage.) The walls matter, mostly through extra
inelastic scattering and gamma production, so the low photon and low neutron
continuum here are largely the simplified model, not TENDL. Running the same
model with ENDF/B-VIII.1 Al27 instead gives 1.02, 0.78 and 0.77, so the two
libraries are close for aluminium.

Also worth knowing:

* The tritium target's own gamma rays (0.0862 per source neutron, Table 3.1 of
  the report) reach the detector too. They are transported through the sphere
  in a second run and included (the dashed line).
* The gamma measurement counted photons arriving within 69 ns of the source
  pulse. yamc runs to steady state, so slow capture gammas are included.
* yamc has no surface current tally, so leakage is tallied as track-length
  flux in a void shell around the sphere, which is exact to 1e-4 here (see
  the docstring of `run_simulation.py`).

## Extending it

See issue #1 for other OKTAVIAN spheres and other benchmarks.
