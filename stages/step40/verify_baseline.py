"""Check all nine previously reconstructed Step 39 reference records."""
from pathlib import Path
import json
import sys

def reproduce(baseline, work):
    root = baseline.parents[1]
    sys.path.insert(0, str(root))
    from compare_records import check_records
    expected = json.loads((root / "reference.json").read_text())["39"]
    count = check_records(baseline / "results/verified", expected)
    return {"core_json_files_matched": count}
