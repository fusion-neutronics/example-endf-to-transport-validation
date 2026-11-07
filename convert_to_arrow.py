"""Step 1: convert the Al27 ENDF evaluations to Arrow.

* Neutrons: TENDL-2025 n-Al027, processed by NJOY to 294 K and written by
  yamc.convert_neutron_transport (cross sections, secondary distributions and
  photon production).
* Photons: TENDL has no photoatomic data, so the gammas are transported with
  the ENDF/B-VIII.1 Al photoatomic and atomic relaxation evaluations.

Archives are downloaded to downloads/ and checked against their sha256.
"""

import hashlib
import pathlib
import shutil
import tarfile
import urllib.request

import yamc

HERE = pathlib.Path(__file__).parent
DOWNLOADS = HERE / "downloads"

ARCHIVES = {  # url, sha256
    "neutron": ("https://tendl.imperial.ac.uk/tendl_2025/tar_files/TENDL-n.tgz",
                "e547527688506cbe09813364dcefa2aed11f474139bfa129d7cd4ca24fae21fa"),
    "photoat": ("https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/photoat/"
                "photoat-version.VIII.1.tar.gz",
                "d39dfc8708b62dd8647bd7ff0af36741f629260d70724f71be13c949c62aa79f"),
    "relax": ("https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/atomic_relax/"
              "atomic_relax-version.VIII.1.tar.gz",
              "95845646a25982aa84d7e4c3459bd51332ee35f345f8d84a9d62181a448ccd23"),
}
# WHAT TO CHANGE FOR ANOTHER MATERIAL: the evaluations to convert.
# * neutron: one TENDL file per nuclide, named n-<Symbol><AAA>.tendl
#   (n-Co059.tendl, n-Fe056.tendl). A natural element needs every isotope.
# * photoat / relax: one ENDF/B-VIII.1 file each per element, named with the
#   atomic number, photoat-<ZZZ>_<Symbol>_000.endf and atom-<ZZZ>_<Symbol>_000.endf.
# The nuclides converted here must match COMPOSITION in run_simulation.py,
# and the elements must match ELEMENTS there.
MEMBERS = {"neutron": "n-Al027.tendl", "photoat": "photoat-013_Al_000.endf",
           "relax": "atom-013_Al_000.endf"}


def fetch(name):
    """Download one archive (once), check it, and extract the Al member."""
    url, digest = ARCHIVES[name]
    archive = DOWNLOADS / name / url.rsplit("/", 1)[1]
    if not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, archive)
    stamp = archive.with_name(archive.name + ".sha256")
    if not stamp.exists():
        h = hashlib.sha256()
        with open(archive, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        if h.hexdigest() != digest:
            raise SystemExit(f"{archive} failed its sha256 check; delete it and rerun")
        stamp.write_text(digest)
    out = DOWNLOADS / "endf" / MEMBERS[name]
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "r|*") as tar:  # one streaming pass
            for member in tar:
                if pathlib.PurePosixPath(member.name).name == MEMBERS[name]:
                    with tar.extractfile(member) as src, open(out, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                    break
    return out


neutron_dir = HERE / "data" / "tendl-2025" / "neutron"
photon_dir = HERE / "data" / "endf-b8.1" / "photon"

# The Arrow names below (Al27, Al) follow from MEMBERS: change them with it.
if not (neutron_dir / "Al27.arrow").exists():
    print("converting TENDL-2025 Al27 (NJOY)")
    yamc.convert_neutron_transport(str(fetch("neutron")), str(neutron_dir),
                                   njoy_exec=shutil.which("njoy") or "njoy",
                                   temperatures=[294.0], library="tendl-2025")
if not (photon_dir / "Al.arrow").exists():
    print("converting ENDF/B-VIII.1 Al photon data")
    yamc.convert_photon(fetch("photoat"), photon_dir, relaxation_path=fetch("relax"),
                        library="endf-b8.1")
print(f"done: {neutron_dir / 'Al27.arrow'}, {photon_dir / 'Al.arrow'}")
