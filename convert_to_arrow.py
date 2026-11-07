"""Step 1: convert the Al27 ENDF evaluations to Arrow.

The three evaluations are committed in nuclear_data/, so nothing is downloaded:

* nuclear_data/n-Al027.tendl: TENDL-2025 neutron evaluation of Al27
  (CC-BY-4.0, tendl.imperial.ac.uk). NJOY processes it to 294 K and
  yamc.convert_neutron_transport writes the cross sections, secondary
  distributions and photon production.
* nuclear_data/photoat-013_Al_000.endf and atom-013_Al_000.endf: the
  ENDF/B-VIII.1 Al photoatomic and atomic relaxation evaluations (NNDC). TENDL
  has no photoatomic data, so the gammas are transported with these.

WHAT TO CHANGE FOR ANOTHER MATERIAL: put its evaluations in nuclear_data/ and
list them below.
* Neutrons: one TENDL file per nuclide, named n-<Symbol><AAA>.tendl
  (n-Co059.tendl, n-Fe056.tendl). A natural element needs every isotope. Take
  them from the TENDL-2025 neutron archive,
  https://tendl.imperial.ac.uk/tendl_2025/tar_files/TENDL-n.tgz
* Photons: one photoatomic and one atomic relaxation file per element, named
  with the atomic number, photoat-<ZZZ>_<Symbol>_000.endf and
  atom-<ZZZ>_<Symbol>_000.endf, from
  https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/
The nuclides must match COMPOSITION in run_simulation.py and the elements must
match ELEMENTS there.
"""

import pathlib
import shutil

import yamc

HERE = pathlib.Path(__file__).parent
ENDF = HERE / "nuclear_data"

NEUTRON_FILES = {"Al27": "n-Al027.tendl"}  # nuclide: TENDL-2025 file
PHOTON_FILES = {"Al": ("photoat-013_Al_000.endf", "atom-013_Al_000.endf")}  # element: files

neutron_dir = HERE / "data" / "tendl-2025" / "neutron"
photon_dir = HERE / "data" / "endf-b8.1" / "photon"

for nuclide, filename in NEUTRON_FILES.items():
    if not (neutron_dir / f"{nuclide}.arrow").exists():
        print(f"converting TENDL-2025 {nuclide} (NJOY)")
        yamc.convert_neutron_transport(str(ENDF / filename), str(neutron_dir),
                                       njoy_exec=shutil.which("njoy") or "njoy",
                                       temperatures=[294.0], library="tendl-2025")

for element, (photoatomic, relaxation) in PHOTON_FILES.items():
    if not (photon_dir / f"{element}.arrow").exists():
        print(f"converting ENDF/B-VIII.1 {element} photon data")
        yamc.convert_photon(ENDF / photoatomic, photon_dir,
                            relaxation_path=ENDF / relaxation, library="endf-b8.1")

print(f"done: Arrow data in {neutron_dir} and {photon_dir}")
