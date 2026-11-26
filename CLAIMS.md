# Finite checks and paper locations

- II §§5–7, App D: word quotient, trace identities and covariance forms. Output directory `37/`.
- II §7, Apps D–E: energy pencil and two energy-strip certificates. Output directory `38/`.
- II §7, Apps D–E: radial substitutions and two continuation certificates. Output directory `39/`.
- II §7, Apps D–E: color-centering, two additional certificates and limiting variance. Output directory `40/`.

Appendix D identifies `37/quotient11.json`, `37/gram_blocks.json` and `37/covariance_structure.json` as inputs to the energy reduction. This program generates them from the word rules. The six Appendix E certificates are the lower/upper energy, refined radial and color-centered continuation inputs. The verification reports are `38/{lower,upper}_certificate_checks.json`, `39/{lower,upper}_uniform_certificate.json` and `40/{lower,upper}_new_uniform.json`. The final support and variance are in `40/radial_consequences.json`. This distribution reconstructs the mathematical objects used by Appendix E.3 in local directories. It does not reconstruct the historical distribution bytes. The failed floating support is provided with its selected atoms and objective for the exact rank negative control.

## Appendix material index

The following objects are named in the paper. Paths beginning with a stage number are relative to the output of `--results`. Regenerated mathematical records are checked against frozen references before any downstream stage uses them.

- Reconstructed record: `37/quotient11.json`.
- Reconstructed record: `37/gram_blocks.json`.
- Reconstructed record: `37/covariance_structure.json`.
- Distributed input or implementation: `stages/step38/energy_model.py`.
- Reconstructed coefficient data: `38/inputs/energy_pencil.json`.
- Distributed input or implementation: `stages/step38/inputs/lower_certificate.json`.
- Reconstructed record: `38/lower_certificate_checks.json`.
- Distributed input or implementation: `stages/step38/inputs/upper_certificate.json`.
- Reconstructed record: `38/upper_certificate_checks.json`.
- Reconstructed coefficient data: `39/inputs/radial_pencil.json`.
- Distributed input or implementation: `stages/step39/inputs/lower_refined_continuation.json`.
- Reconstructed record: `39/lower_uniform_certificate.json`.
- Distributed input or implementation: `stages/step39/inputs/upper_refined_continuation.json`.
- Reconstructed record: `39/upper_uniform_certificate.json`.
- Reconstructed coefficient data: `40/inputs/augmented_pencil.json`.
- Distributed input or implementation: `stages/step40/inputs/lower_continuation.json`.
- Reconstructed record: `40/lower_new_uniform.json`.
- Distributed input or implementation: `stages/step40/inputs/upper_continuation.json`.
- Reconstructed record: `40/upper_new_uniform.json`.

With a separately supplied manuscript LaTeX source, run the optional citation and cross-reference checks after reproduction:

```sh
python3 -B -S tools/check_manuscript.py manuscript.tex reproduced
```
