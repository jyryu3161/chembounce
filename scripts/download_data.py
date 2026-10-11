#!/usr/bin/env python3
"""Download the versioned Zenodo library, check its published digest, extract safely."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import urllib.request
import zipfile

URL = "https://zenodo.org/records/16741967/files/data.zip?download=1"
MD5 = "f4e0439e44da841d822099c10766d099"  # Published by Zenodo record 16741967.


def install(destination):
    destination = Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f"Destination already exists: {destination}; choose an empty path")
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".chembounce-data-") as temp:
        temp = Path(temp)
        archive = temp / "data.zip"
        md5, sha = hashlib.md5(), hashlib.sha256()
        with urllib.request.urlopen(URL, timeout=60) as response, archive.open("wb") as stream:
            for block in iter(lambda: response.read(1024 * 1024), b""):
                stream.write(block)
                md5.update(block)
                sha.update(block)
        if md5.hexdigest() != MD5:
            raise ValueError("Downloaded archive does not match the published Zenodo checksum")
        extracted = temp / "extracted"
        extracted.mkdir()
        with zipfile.ZipFile(archive) as zipped:
            for member in zipped.infolist():
                path = (extracted / member.filename).resolve()
                if (
                    not path.is_relative_to(extracted)
                    or ((member.external_attr >> 16) & 0o170000) == 0o120000
                ):
                    raise ValueError("Unsafe archive member")
            zipped.extractall(extracted)
        root = extracted / "data" if (extracted / "data").is_dir() else extracted
        (root / "download_manifest.json").write_text(
            json.dumps({"url": URL, "md5": md5.hexdigest(), "sha256": sha.hexdigest()}, indent=2) + "\n"
        )
        shutil.move(str(root), destination)
    print(f"Installed reference data: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(Path(__file__).resolve().parents[1] / "data"))
    args = parser.parse_args()
    install(args.output)
