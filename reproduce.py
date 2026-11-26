"""Rebuild the finite certificates using only the distributed inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
from paper_values import check_paper_values

from compare_records import check_records

def run(stage, script, *args):
    subprocess.run([sys.executable, "-B", "-S", script, *map(str, args)], cwd=stage,
                   env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, check=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, help="Write the verified JSON records to a new directory")
    args = parser.parse_args()
    if args.results and args.results.exists():
        parser.error("Choose a new results directory")
    plan = json.loads((ROOT / "reference.json").read_text())
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        stages = work / "stages"
        shutil.copytree(ROOT / "stages", stages)
        shutil.copyfile(ROOT / "reference.json", work / "reference.json")
        checked = {}
        generated_inputs = json.loads((ROOT / "generated_inputs.json").read_text())
        shutil.copyfile(ROOT / "compare_records.py", work / "compare_records.py")
        for label, expected in sorted(plan.items(), key=lambda row: int(row[0]) if row[0].isdigit() else -1):
            if not label.isdigit():
                continue
            n = int(label)
            stage = stages / ("step" + label)
            out = stage / "results/verified"
            out.mkdir(parents=True, exist_ok=True)
            if n == 38:
                run(stage, "energy_model.py", "--source-dir", checked["37"], "--output", stage / "inputs/energy_pencil.json")
            elif n == 39:
                run(stage, "radial_model.py", "--source", stages / "step38/inputs/energy_pencil.json", "--out", stage / "inputs")
            elif n == 40:
                run(stage, "augment_model.py", "--baseline", stages / "step39", "--output-dir", stage / "inputs")
            for filename, expected_hash in generated_inputs.get(label, {}).items():
                if hashlib.sha256((stage / "inputs" / filename).read_bytes()).hexdigest() != expected_hash:
                    raise ValueError("Reconstructed coefficient data differ: " + filename)
            run(stage, f"verify_step{n}.py", "--output-dir", out)
            checked[label] = out
            count = check_records(out, expected)
            print(f"Step {n}: {count} records matched.", flush=True)
        check_paper_values(checked)
        if args.results:
            args.results.mkdir(parents=True)
            for label, directory in checked.items():
                shutil.copytree(directory, args.results / label)
            for label, files in generated_inputs.items():
                destination = args.results / label / "inputs"
                destination.mkdir()
                for name in files:
                    shutil.copyfile(stages / ("step" + label) / "inputs" / name, destination / name)
        print("All finite reference comparisons and paper-value checks passed.")

if __name__ == "__main__":
    main()
