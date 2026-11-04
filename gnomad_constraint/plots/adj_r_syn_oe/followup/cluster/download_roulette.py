"""Copy Roulette's per-chromosome rate VCFs (v5.2) from the BWH download site to GCS.

Run as a Dataproc job: each file is streamed from the site straight into the bucket
(curl into ``gsutil cp -``), several chromosomes at a time, and its size is checked
against the site's. Files already in the bucket with the right size are skipped.
"""

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor

URL = "https://genetics.bwh.harvard.edu/downloads/Vova/Roulette"
NAME = "{c}_rate_v5.2_TFBS_correction_all.vcf.bgz"


def remote_size(url: str) -> int:
    """Get a file's size on the download site."""
    out = subprocess.run(
        ["curl", "-sSfIL", url], capture_output=True, text=True, check=True
    ).stdout
    sizes = [
        int(line.split(":")[1])
        for line in out.splitlines()
        if line.lower().startswith("content-length")
    ]
    return sizes[-1]


def bucket_size(path: str) -> int:
    """Get a file's size in the bucket (-1 if it isn't there)."""
    out = subprocess.run(["gsutil", "du", path], capture_output=True, text=True)
    fields = out.stdout.split()
    return int(fields[0]) if out.returncode == 0 and fields else -1


def copy(name: str, dest: str) -> str:
    """Stream one file into the bucket, retrying, and check its size."""
    url, path = f"{URL}/{name}", f"{dest}/{name}"
    size = remote_size(url)
    for attempt in range(4):
        if bucket_size(path) == size:
            return f"{name}: {size} bytes"
        curl = subprocess.Popen(["curl", "-sSfL", url], stdout=subprocess.PIPE)
        subprocess.run(["gsutil", "-q", "cp", "-", path], stdin=curl.stdout)
        curl.wait()
    raise RuntimeError(f"{name}: size still differs from {size} after 4 tries!")


def main(args: argparse.Namespace) -> None:
    """Copy every file."""
    names = [NAME.format(c=c) for c in args.chromosomes]
    names += [f"{n}.csi" for n in names]
    with ThreadPoolExecutor(args.n_parallel) as pool:
        for msg in pool.map(lambda n: copy(n, args.dest), names):
            print(msg, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", required=True, help="GCS folder for the files.")
    parser.add_argument(
        "--chromosomes", nargs="+", default=[str(c) for c in range(1, 23)]
    )
    parser.add_argument("--n-parallel", type=int, default=8)
    main(parser.parse_args())
