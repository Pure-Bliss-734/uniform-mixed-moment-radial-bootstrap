# Uniform Mixed-Moment Bounds and an Unfactorized Radial Bootstrap for Massive Matrix Quantum Mechanics

How much can the squared radius of a massive two-matrix quantum ground state fluctuate as the matrix size grows? The paper proves uniform bounds on coordinate and momentum moments. At unit mass squared and coupling, it shows that every limiting distribution of the normalized squared radius is supported in **[0.38895, 0.38910]**, without assuming that products of traces factorize. This repository reconstructs the six finite certificates behind that interval and checks the resulting variance bound of **5.625 × 10⁻⁹**. The analytical proofs are in the [paper](https://doi.org/10.5281/zenodo.24453945).

**Author disclosure.** The author is not a professional physics researcher and initiated this project out of personal interest. Large language models (LLMs) performed all the research, including the mathematical derivations and computations, and wrote this manuscript. This work has not undergone peer review.
 
## Reproduce

From the extracted directory, use Python 3.11 or newer with assertions enabled. The finite computations use only the standard library and temporary disk space. To retain the reconstructed scientific records for inspection, supply a new output directory.

```sh
python3 -B -S reproduce.py --results reproduced
```

For an archive check and a full computational run, first extract to a clean directory, then run:

<!-- RELEASE_COMMAND:verify -->
```sh
python3 -B -S tools/verify_public.py --archive <ARCHIVE_PATH>
```
<!-- END_RELEASE_COMMAND:verify -->

## Records and comparisons

See [the paper-to-record guide](CLAIMS.md). The four stages check 9, 11, 9 and 15 reference records, respectively, including the six strip endpoints in Appendices D and E. The input coefficients are integers and exact rational numbers. The numerical certificate checks use exact arithmetic or outward-rounded intervals. The recorded reference digests use sorted, compact UTF-8 JSON. The `--results` directory contains the regenerated certificate records grouped by stage.

## License

Project-authored software is licensed under MIT. Contributor-owned research materials are licensed under CC BY 4.0. Upstream terms apply to third-party materials. See [LICENSE](LICENSE).
