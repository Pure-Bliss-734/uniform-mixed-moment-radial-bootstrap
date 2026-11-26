from fractions import Fraction
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_bytes())

def check_paper_values(paths):
    table = read(ROOT / "paper_values.json")
    for stage, field in [("38","limiting_radius_support"), ("39","limiting_support"), ("40","necessary_limiting_support")]:
        values = read(paths[stage] / "radial_consequences.json")[field]
        if len(values) != 2 or any(Fraction(printed) != Fraction(values[index])
                                   for index, pair in enumerate(table[stage]) for printed in pair):
            raise ValueError("Radius endpoints differ from the paper: Step " + stage)
    variance = read(paths["40"] / "radial_consequences.json")["limsup_variance_upper"]
    if Fraction(variance) != Fraction(table["variance"]):
        raise ValueError("Variance bound differs from the paper")
    print("Matched the six strip endpoints and the variance bound.")
