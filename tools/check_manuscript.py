"""Check manuscript citations, cross-references and displayed finite certificates."""
from fractions import Fraction
import json
from pathlib import Path
import re


def check(source, energy, radial, centered):
    text = source.read_text()
    labels = re.findall(r'\\label\{([^}]+)\}', text)
    refs = re.findall(r'\\(?:eqref|ref|pageref)\{([^}]+)\}', text)
    if len(labels) != len(set(labels)) or set(refs) - set(labels):
        raise ValueError('Unresolved or duplicated manuscript label')
    bib = re.findall(r'\\bibitem\{([^}]+)\}', text)
    cites = {key.strip() for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',text) for key in group.split(',')}
    if len(bib) != len(set(bib)) or cites != set(bib):
        raise ValueError('Bibliography and citations differ')

    def equation(label):
        return text.split(r'\label{'+label+'}',1)[1].split(r'\end{equation}',1)[0]

    for label, directory, field in (
            ('eq:energycertificateradii', energy, 'limiting_radius_support'),
            ('eq:radialcertificateradii', radial, 'limiting_support'),
            ('eq:finalsupportendpoints', centered, 'necessary_limiting_support')):
        facts=json.loads((directory/'radial_consequences.json').read_text())
        endpoints=list(map(Fraction,facts[field]))
        printed=re.findall(r'\\frac\{(\d+)\}\{(\d+)\}\s*=\s*(0\.\d+)',equation(label))
        if [Fraction(int(a),int(b)) for a,b,c in printed] != endpoints or [Fraction(c) for a,b,c in printed]!=endpoints:
            raise ValueError('Radius endpoint transcription mismatch in '+label)
    variance=json.loads((centered/'radial_consequences.json').read_text())['limsup_variance_upper']
    m=re.search(r'\\frac(\d)\{(\d+)\}=([\d.]+)\\times10\^\{(-?\d+)\}',equation('eq:certifiedvariance'))
    if not m or Fraction(int(m[1]),int(m[2]))!=Fraction(variance) or Fraction(m[3])*Fraction(10)**int(m[4])!=Fraction(variance):
        raise ValueError('Variance transcription mismatch')
    print(f"Checked {len(labels)} labels, {len(bib)} references, six endpoints and the variance bound.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manuscript", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    check(args.manuscript, args.results / "38", args.results / "39", args.results / "40")
