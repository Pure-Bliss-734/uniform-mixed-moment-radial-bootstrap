"""Verify the manuscript distribution and reproduce its finite certificate records."""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import tempfile
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
PREFIX = "mixed-moment-certificates-v2.0.0"
PRIVATE_LOCATION = re.compile(rb"(?:/(?:Users|home|Volumes)/|[A-Z]:[\\/]Users[\\/])[^\s/\\'\"]+", re.I)
CREDENTIAL = re.compile(rb"(?:-----BEGIN (?:[A-Z]+ )*PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})")


def check_archive(path):
    sums = (ROOT / "SHA256SUMS.txt").read_bytes()
    expected = {"SHA256SUMS.txt": hashlib.sha256(sums).hexdigest()}
    for line in sums.decode("utf-8").splitlines():
        digest, name = line.split("  ", 1)
        if name in expected:
            raise ValueError("Duplicate checksum entry")
        expected[name] = digest
    seen = set()
    names = set()
    with ZipFile(path, "r") as archive:
        entries = archive.infolist()
        if len(entries) > 2000 or sum(item.file_size for item in entries) > 300_000_000:
            raise ValueError("Unexpected archive size")
        for item in entries:
            member = PurePosixPath(item.filename)
            if (not member.parts or member.parts[0] != PREFIX or ".." in member.parts
                    or chr(92) in item.filename or str(member) != item.filename.rstrip("/")
                    or member in names or stat.S_ISLNK(item.external_attr >> 16)
                    or any(part in {".git", ".DS_Store", "__MACOSX"} or part.startswith("._") for part in member.parts)):
                raise ValueError("Unsafe archive member")
            names.add(member)
            if item.is_dir():
                continue
            name = member.relative_to(PREFIX).as_posix()
            if name in seen or name not in expected:
                raise ValueError("Unexpected or duplicate member")
            local = ROOT / name
            if local.is_symlink() or not local.resolve().is_relative_to(ROOT):
                raise ValueError("Unsafe extracted file")
            content = archive.read(item)
            if PRIVATE_LOCATION.search(content) or CREDENTIAL.search(content):
                raise ValueError("Local path or credential detected in " + name)
            if hashlib.sha256(content).hexdigest() != expected[name]:
                raise ValueError("ZIP member checksum mismatch: " + name)
            if hashlib.sha256(local.read_bytes()).hexdigest() != expected[name]:
                raise ValueError("Extracted checksum mismatch: " + name)
            seen.add(name)
    if seen != expected.keys():
        raise ValueError("Archive member missing")
    local_files = {file.relative_to(ROOT).as_posix() for file in ROOT.rglob("*") if file.is_file() and ".git" not in file.relative_to(ROOT).parts}
    if local_files != seen:
        raise ValueError("Extracted file inventory differs from the ZIP")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    check_archive(args.archive)
    subprocess.run([sys.executable, "-B", "-S", "reproduce.py"], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
