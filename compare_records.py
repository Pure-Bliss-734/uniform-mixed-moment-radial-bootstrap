import hashlib
import json

def check_records(directory, reference):
    actual = {p.stem for p in directory.glob("*.json")}
    if actual != set(reference):
        raise ValueError("Reference and computed record sets differ in " + str(directory))
    for name, expected in reference.items():
        obj = json.loads((directory / (name + ".json")).read_bytes())
        payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError("Computed record differs: " + str(directory / name))
    return len(reference)

