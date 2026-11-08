"""Step 1: convert the TENDL-2025 ENDF evaluations to Arrow.

The TENDL-2025 neutron evaluations (CC-BY-4.0, tendl.imperial.ac.uk) of the
sphere samples, Al27 and Mn55. NJOY processes each to 294 K
and yamc.convert_neutron_transport writes the cross sections, secondary
distributions and photon production yamc transports with. About 20 seconds
per nuclide, once.

Only the sample's evaluations are converted here. Everything else the model
needs (the photon interaction data, and the neutron data of the steel vessel)
is not what the benchmark tests, so run_simulation.py loads it already
processed, by yamc library keyword.

Each file is looked for in nuclear_data/ (or --endf-dir, searched
recursively), and downloaded there if it is missing (TENDL serves each
evaluation on its own, about 4 MB each). Found or downloaded, it must match
its pinned sha256, so the reference case is guaranteed to be exactly
TENDL-2025. To test a different evaluation, pass it to run_simulation.py with
--neutron-file and --label instead.

    python convert_to_arrow.py
    python convert_to_arrow.py --endf-dir ~/nuclear_data/tendl-2025

WHAT TO CHANGE FOR ANOTHER MATERIAL: add its files to NEUTRON_FILES. One
TENDL file per nuclide, n-<Symbol><AAA>.tendl (n-Co059.tendl); a natural
element needs every isotope. Set the sha256 to None until you have checked a
first download, then pin it.
"""

import argparse
import hashlib
import pathlib
import shutil
import urllib.request

import yamc

HERE = pathlib.Path(__file__).parent
TENDL = "https://tendl.imperial.ac.uk/tendl_2025/neutron_files"

# nuclide: (file, sha256)
NEUTRON_FILES = {
    "Al27": ("n-Al027.tendl",
             "61afaf01ac4fca5b5d43bfcedaffaa66b5f36dddb9c1c657b4b0a21ea63d2ffd"),
    "Mn55": ("n-Mn055.tendl",
             "364707f2d85df022559021802eaa583337018017cc964c6dd7e1a66fd0d861a4"),
}


def tendl_url(nuclide, filename):
    symbol = nuclide.rstrip("0123456789")
    return f"{TENDL}/{symbol}/{filename[2:-6]}/lib/endf/{filename}"  # n-Al027.tendl -> Al027


def local_file(nuclide, filename, sha256, endf_dir):
    """The path of an evaluation, downloading it into endf_dir if it is missing.

    The file, found or downloaded, must match its pinned sha256: a different
    release under the same name, or a truncated download, is refused.
    """
    path = next(iter(sorted(endf_dir.rglob(filename))), None)
    downloaded = path is None
    if downloaded:
        path = endf_dir / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        url = tendl_url(nuclide, filename)
        print(f"  downloading {url}")
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request) as response, open(path, "wb") as out:
            shutil.copyfileobj(response, out)
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

    for nuclide, (filename, sha256) in NEUTRON_FILES.items():
        if (neutron_dir / f"{nuclide}.arrow").exists():
            continue
        endf = local_file(nuclide, filename, sha256, endf_dir)
        print(f"converting {endf} (NJOY)")
        yamc.convert_neutron_transport(str(endf), str(neutron_dir),
                                       njoy_exec=shutil.which("njoy") or "njoy",
                                       temperatures=[294.0], library="tendl-2025")
    print(f"done: Arrow data in {neutron_dir}")


if __name__ == "__main__":
    main()
