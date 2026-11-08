# Results that were only on the desktop, for checking against the papers (2026-11-05)

Copied as they were, from the runs of October 2026; only large files were left out (over 300 KB: the read fragments;
the `.lvl` indexes). Nothing here was rerun.

## `shortcuts_chr19/`: Reportable Shortcuts on chr19

From `~/research/rlbwt-tricks/shortcuts/` (`run_chr19_scaling.sh`, its log) and `~/research/rlbwt-tricks/shortcuts/chr19_scaling/hprc_chr19/`:
- `summary_scaling.tsv`: all sizes together; `n16/`, `n32/`, `n64/summary.tsv`: each size;
- `n*/mems/t<tau>_ell<ell>_{bml,fb}.csv`: the MEM CSVs;
- `manifest.tsv`: the sequences.

There is no `results.csv` or `summary.txt` for this run on the desktop or on the external drive. `summary_scaling.tsv`
and `summary.tsv` are its summaries. (`scripts/collect_results.py` in `reportable-shortcuts-code.tar.gz` may be what
produced `results.csv` elsewhere.)

## `twolevel_old/`: the old two-level code (More RLBWT Tricks)

- **The final runs:** `final.log`, `extras_final.log`, `movi_*.log`, `ropebwt3_n16.log`, `morning.log`,
  `overnight.log`, `final/summary*.txt`, and `final/n*/{bml,fb}.{csv,txt}`, from `~/twolevel_runs/twolevel/`
  (on the external drive).
- **The most frequent successor:** `continuation.log` and `continuation/` (build log; `q_{smallest,frequent}_{bml,fb}`),
  from `~/research/hybrid-fm-csa/twolevel/`, where `run_continuation.sh` was run on Oct 28 at n16 (it was not in
  `~/twolevel_runs`).
- **The simulated noisy reads:** `reads.log`, `reads/summary.txt`, `reads/n{16,64}/e{1,5}/*.txt`, `sim_reads.py` and
  `run_reads.sh`. The fragment files `frag_L*.txt` (1–20 MB each) stay on the external drive, in
  `twolevel_runs/twolevel/reads/`.
