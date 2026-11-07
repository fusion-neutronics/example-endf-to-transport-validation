"""Step 1: convert the Al27 ENDF evaluations to Arrow.

Three evaluations are needed:

* n-Al027.tendl: the TENDL-2025 neutron evaluation of Al27 (CC-BY-4.0,
  tendl.imperial.ac.uk). NJOY processes it to 294 K and
  yamc.convert_neutron_transport writes the cross sections, secondary
  distributions and photon production.
* photoat-013_Al_000.endf and atom-013_Al_000.endf: the ENDF/B-VIII.1 Al
  photoatomic and atomic relaxation evaluations (NNDC). TENDL has no
  photoatomic data, so the gammas are transported with these.

Each file is looked for in the local ENDF directory, nuclear_data/ unless
--endf-dir names another, so a copy you already have saves the download. A
file that is not there is downloaded into nuclear_data/ first. Every file,
found or downloaded, must match its pinned sha256, so the reference case is
guaranteed to be exactly TENDL-2025 and ENDF/B-VIII.1. To test a different
evaluation, pass it to run_simulation.py with --neutron-file and --label
instead.

* TENDL serves each evaluation on its own (3.9 MB for Al27).
* The NNDC serves only whole sublibrary archives, so the archive is streamed
  and the download stops as soon as the wanted file has been read: 1.8 MB of
  the 34 MB photoatomic archive and 20 kB of the atomic relaxation one, since
  the archives are ordered by atomic number and Al is Z = 13.

    python convert_to_arrow.py                          # download what is missing
    python convert_to_arrow.py --endf-dir ~/nuclear_data/tendl-2025

WHAT TO CHANGE FOR ANOTHER MATERIAL: add its files to NEUTRON_FILES and
PHOTON_FILES below.
* Neutrons: one TENDL file per nuclide, n-<Symbol><AAA>.tendl (n-Co059.tendl).
  A natural element needs every isotope. The download URL follows the pattern
  neutron_files/<Symbol>/<Symbol><AAA>/lib/endf/<file>.
* Photons: one photoatomic and one atomic relaxation file per element, named
  with the atomic number, photoat-<ZZZ>_<Symbol>_000.endf and
  atom-<ZZZ>_<Symbol>_000.endf.
Set each sha256 to None until you have checked a first download, then pin it.
The nuclides must match COMPOSITION in run_simulation.py and the elements
must match ELEMENTS there.
"""

import argparse
import hashlib
import pathlib
import shutil
import tarfile
import urllib.request

import yamc

HERE = pathlib.Path(__file__).parent

TENDL = "https://tendl.imperial.ac.uk/tendl_2025/neutron_files"
NNDC = "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1"
PHOTOAT_ARCHIVE = f"{NNDC}/photoat/photoat-version.VIII.1.tar.gz"
RELAX_ARCHIVE = f"{NNDC}/atomic_relax/atomic_relax-version.VIII.1.tar.gz"

# nuclide: (file, sha256)
NEUTRON_FILES = {
    "Al27": ("n-Al027.tendl",
             "61afaf01ac4fca5b5d43bfcedaffaa66b5f36dddb9c1c657b4b0a21ea63d2ffd"),
}
# element: ((photoatomic file, sha256), (atomic relaxation file, sha256))
PHOTON_FILES = {
    "Al": (("photoat-013_Al_000.endf",
            "de8205f2315aa9b5863fe93912bbc0e8cb2cc53a0995a490196b15f114c765c9"),
           ("atom-013_Al_000.endf",
            "a0e4c6c6320949d8387df60ecf5f8c706cb9c70b115bca654fd36a620455765d")),
}


def _open(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}))


def tendl_url(nuclide, filename):
    symbol = nuclide.rstrip("0123456789")
    return f"{TENDL}/{symbol}/{filename[2:-6]}/lib/endf/{filename}"  # n-Al027.tendl -> Al027


def download(filename, dest, url=None, archive=None):
    """Fetch one file, either directly or by streaming it out of an archive."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    if url:
        print(f"  downloading {url}")
        with _open(url) as response, open(partial, "wb") as out:
            shutil.copyfileobj(response, out)
    else:
        print(f"  streaming {filename} out of {archive}")
        with _open(archive) as response, tarfile.open(fileobj=response, mode="r|gz") as tar:
            for member in tar:  # stops reading the network once the file is found
                if pathlib.PurePosixPath(member.name).name == filename:
                    with tar.extractfile(member) as src, open(partial, "wb") as out:
                        shutil.copyfileobj(src, out)
                    break
            else:
                raise SystemExit(f"{archive} has no {filename}")
    partial.rename(dest)


def local_file(filename, sha256, endf_dir, **source):
    """The path of an evaluation, downloading it into endf_dir if it is missing.

    The file, found or downloaded, must match its pinned sha256: a different
    release under the same name, or a truncated download, is refused.
    """
    path = next(iter(sorted(endf_dir.rglob(filename))), None)
    downloaded = path is None
    if downloaded:
        path = endf_dir / filename
        download(filename, path, **source)
    found = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha256 and found != sha256:
        if downloaded:
            path.unlink()
        raise SystemExit(
            f"{path} has sha256 {found}, expected {sha256}: it is not the pinned "
            "evaluation. To test a different evaluation, pass it to run_simulation.py "
            "with --neutron-file and --label.")
    if not sha256:
        print(f"  {filename} has sha256 {found}: check it, then pin it in this script")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--endf-dir", type=pathlib.Path, default=HERE / "nuclear_data",
                        help="directory searched (recursively) for the evaluations, "
                             "and where missing ones are downloaded (default: nuclear_data/)")
    args = parser.parse_args()
    endf_dir = args.endf_dir.expanduser()

    neutron_dir = HERE / "data" / "tendl-2025" / "neutron"
    photon_dir = HERE / "data" / "endf-b8.1" / "photon"

    for nuclide, (filename, sha256) in NEUTRON_FILES.items():
        if (neutron_dir / f"{nuclide}.arrow").exists():
            continue
        endf = local_file(filename, sha256, endf_dir, url=tendl_url(nuclide, filename))
        print(f"converting {endf} (NJOY)")
        yamc.convert_neutron_transport(str(endf), str(neutron_dir),
                                       njoy_exec=shutil.which("njoy") or "njoy",
                                       temperatures=[294.0], library="tendl-2025")

    for element, ((photoat, photoat_sha), (relax, relax_sha)) in PHOTON_FILES.items():
        if (photon_dir / f"{element}.arrow").exists():
            continue
        photoat_path = local_file(photoat, photoat_sha, endf_dir, archive=PHOTOAT_ARCHIVE)
        relax_path = local_file(relax, relax_sha, endf_dir, archive=RELAX_ARCHIVE)
        print(f"converting {photoat_path} and {relax_path}")
        yamc.convert_photon(photoat_path, photon_dir, relaxation_path=relax_path,
                            library="endf-b8.1")

    print(f"done: Arrow data in {neutron_dir} and {photon_dir}")


if __name__ == "__main__":
    main()
