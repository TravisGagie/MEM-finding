# MEM-finding

Code and paper source for *Three Examples of MEM-Finding* (Travis Gagie and Gonzalo Navarro, in preparation). The paper describes three experiments in finding maximal exact matches (MEMs) with indexes based on the run-length compressed Burrows–Wheeler transform. Each experiment has a different bottleneck:

| Paper section | Experiment | Bottleneck | Code |
|---|---|---|---|
| 3 | ONT reads of the ZymoBIOMICS mock community, classified to species (as in SPUMONI 2) | finding the right matches | `katka/` |
| 4 | Simulated 16S rRNA reads, classified to genera with SILVA (as in Cliffy) | answering each MEM | `katka/` |
| 5 | Two chromosome 19 haplotypes against 16–64 HPRC individuals | crossing long matches | `chr19/` |

## Layout

- `paper/`: LaTeX source (LIPIcs). Build with `pdflatex main; bibtex main; pdflatex main; pdflatex main`.
- `katka/`: KATKA, the classifier of Sections 3 and 4 and Appendices A–C: `rz-index/` (the `rz-*` programs), `silva/` (SILVA experiment scripts), `tools/` (scorers) and the `demo_*.sh` scripts. Appendix A of the paper says how to reproduce each result. This is the code previously at `github.com/TravisGagie/KATKA` (commit `65f236e`).
- `katka-readable/`: a cleaned-up rewrite of KATKA's main pipeline (`katka-*` programs), verified to give byte-identical outputs to the original on test data and on the full Zymo reads. Its README refers to sections of the standalone KATKA paper (arXiv 2610.10500), not of this paper.
- `chr19/`: shortcuts, two-level indexing and block jumps (Section 5 and Appendix D): `shortcuts/` (RLBWT and shortcuts, with libsais) and `twolevel/` (the two-level index and the experiment scripts `run_final.sh`, `run_extras_final.sh`, `run_movi_final.sh`, `run_continuation.sh`). This is the code previously at `github.com/TravisGagie/more-rlbwt-tricks` (commit `833de7f`); its README refers to the standalone paper *More RLBWT Tricks* (arXiv 2609.06596).

The two codebases have separate build systems and are independent of each other. The full histories are in the original repositories.

## License

GPL-3.0 (see `LICENSE` and `katka/LICENSE`).
