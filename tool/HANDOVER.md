# Handover from the cloud code session (2026-11-05)

For the chat (paper) session and the desktop code session, who take over from here. Everything below is on
`origin/tool` at the commit that adds this file. The cloud session's container is ephemeral; nothing of its own
remains outside the repository.

## Who does what from now

The cloud session had three roles; they need new owners.
- **Code changes to `tool/`.** This was the cloud session's alone. From now on the desktop code session owns them.
  It should keep the cloud session's checks (below) and push as it goes.
- **Running experiments.** The desktop, as before.
- **The papers' text, and checking numbers against results.** The chat session, using `tool/results/SILVA_DIGEST.md`
  and the results READMEs as the source of truth. A number not traceable to a file and line there should not go into
  a paper.

To avoid the confusion of three sessions on one paper:
- One session writes each file. The chat session writes paper text only. The desktop writes everything under `tool/`,
  including `NOTES.md` and the results READMEs.
- Requests between sessions go through Travis, or as a dated "For the desktop" / "For the chat session" item at the top
  of `tool/NOTES.md`. That is the only place both sessions are sure to read.

## Standing rules (from Travis)

- Don't change `katka/`, `katka-readable/` or `chr19/`. They are the code of published papers. (`paper/`, the abandoned
  *Three Examples* source, was removed in `3647b68`.)
- Don't work on compressing the phrase dictionary beyond packing it; that idea is someone else's.
- Exact intervals or nothing: every configuration must give the same MEMs and intervals as plain search. No
  approximate shortcut goes in, even as an option. If something can't be made exact, write down why in `NOTES.md`.
- Push to the `tool` branch as you go. `main` was brought up to date with `tool` once (`3647b68`). The papers link to the
  `tool` branch by name, so keep it.
- Desktop limits: 11 GB memory cap per step (`systemd-run`), single-threaded timed searches pinned to CPU 9.

## Checks to run after any code change

- `tool/selftest.sh` (about 9 minutes): MEMs and intervals against brute force for every search mode, both continuation
  rules, tags, digests, STPD, locating. It must print `self-test passed`.
- `python3 tool/tests/test_score_schemes.py EXP TAX D2T T2L`: 44 checks of the credit schemes.
  - Covered: every scheme against an independent implementation, the endpoints (β = 0 and 1, EM with 0 iterations,
    the priors), hand-worked cases, and `bayesauto`'s fit.
  - Arguments: EXP is cliffy-experiments' `src`; TAX is `tax_slv_ssu_138.1.txt`; D2T and T2L are any reference's
    `doc_to_traversal.txt` and `trav_to_length.txt` with at least 12 genera.
- `bash tool/tests/dial_toy.sh MICROBEMIXER EXP D2T WORK` (about 1 minute): the whole dial pipeline on a toy SILVA,
  the endpoints, and `find --classify bayes` against the scorer, call for call.
- For block changes: `katka-test mems PREFIX queries --twolevel --random 20` on the self-test's chr19-like collection,
  with each `--continuation` and `--max-block`. `experiments/block_scaling.sh` also checks that every variant gives
  identical MEMs.

All four passed on 2026-11-05, the last two after the desktop's `4f68295` change to `score_silva.py`.

## The papers

### *More RLBWT Tricks* v2 (arXiv:2609.40299; code `tool/`, reproduction `tool/CHR19.md`)
- **Posted.** The cloud session checked every number against `tool/results/chr19*`. Most match. The mismatches are
  listed in the conversation's report, with files and lines.
  - **Peaks in mixed units:** six peaks were converted as KiB/10⁶ instead of KiB × 1024/10⁹, so they are about 2.4% low.
    Table 1's 16-individual build: 4.6 GB, should be 4.7. Table 2's STPD over the text: 4.95, should be 5.06. Two-level
    with locating: 2.22, should be 2.28. The published STPD-index: 0.96, should be 0.98. 32-individual STPD text: 8.87,
    should be 9.08. "25 to 31% less memory" should be 24 to 29%.
  - **Rounding and arithmetic:** 70.2 s should be 70.1; 770 MB should be 771; 1,982 MB should be 1,983.
  - **Table 1's build peaks aren't like for like:** 16 individuals used the old two-level build, 32 and 64 the new one.
  - **Not checked:** §4's "3.3 times longer, 5–6% slower" for the most frequent successor, and the Discussion's
    simulated-read figures. Their raw logs are now in `tool/results/desktop_archive/` (`353ff51`).
  - **Owner:** the chat session, for the corrections; the desktop or chat session, for the unchecked figures.
- **Follow-up work on blocks** (`NOTES.md`, "Block length as genomes are added"; `results/block_scaling/`):
  - **Block length:** under the smallest successor, blocks shrink as genomes are added; under the most frequent they
    stay long (synthetic, 16 to 1,024 haplotypes; chr19, 16 to 64 individuals, about 30 against 7–10 phrases).
  - **Search time:** the most frequent successor doesn't make search faster. It takes more parse steps, all at match
    ends.
  - **Capped blocks** (`index --max-block B`) are exact and checked, but no better.
  - **Not implemented, analysed in NOTES:** a range start in Δ isn't exact as proposed, because of the search's
    right-end tail. One-step match ends need phrase texts, Ψ on Δ and a co-lex wavelet tree.

### *Reportable Shortcuts with RLBWTs* v3 (arXiv:2609.06596; code `chr19/shortcuts`)
- **Posted.** Only internal consistency was checked: text against tables.
  - "1.3 to 1.7" in §7.2 should be 1.3 to 1.6 (Table 2).
  - "8.0 times" in §7.3 is 8.1 (348.7 / 43.1).
  - "10.0 MiB" in §7.1 is 9.9 in Table 1, possibly rounding.
- **Not checked against raw results.** The chr19 results are now in `tool/results/desktop_archive/shortcuts_chr19/`.
  The cere and bacterial results are on the laptop and not yet committed. **Owner:** the chat session, once the laptop's
  results are in.

### *KATKA* (arXiv:2610.10500; code `katka/`, `katka-readable/`)
Posted. Superseded by `tool/`'s `zymo` and `silva` presets. Each repository's README says so.

### The SILVA paper against Cliffy (in progress)
- **Source of truth:** `tool/results/SILVA_DIGEST.md`, every number with its file and line and protocol, current up to
  the desktop's `6a6bd22`. Also `tool/results/silva/README.md`.
- **Not yet in the digest** (the desktop's `4f68295`, `4907378`; `results/silva/README.md`, from line 179):
  - **Timing:** 5 repeats give a mean per-read time of **45.57 µs**, range 45.46–45.87, and **167,224 reads/s** with
    12 threads. These replace the queue's single 52.3 µs and 145,500 reads/s; the slower run's cause is unknown.
  - **Smaller ε, and the ε → 0 limit, for naive Bayes with a flat prior:** on the aquatic reads, 94.70% (ε = 10⁻⁶) down
    to 94.08% (the limit). That is above proportional (93.77%) and normalised (92.76%) at every ε.
  - **Accuracy excluding reads from genera with ambiguous names:** the aquatic reads have none. The simulated tables
    fall by 0.1–1.5 points, e.g. uniform proportional 88.09 → 87.28, normalised 91.75 → 91.20.
  - **`bayesauto` on the grid to 0.99:** ε = 0.5 in 9 of 12 runs and 0.6 in 3, none at the edge, so it is a real fit
    now. The grid is coarse there (0.3, 0.5, 0.6, 0.7). Accuracy 91.42, 94.26, 95.19 at α = 1, 0, −1.

  - **And `aa8294f`, which landed while this was written:** `bayesauto` on the aquatic reads, naive Bayes on the
    uniform reads, timing with one thread and all threads, and what each SILVA index size counts. All are in
    `results/silva/README.md`.

  **Owner:** the desktop, to update the digest, or the chat session, to read these from the README directly.
- **Main findings so far** (details and sources in the digest):
  - **Proportional credit wins on Cliffy-style reads** (93.77% strict, aquatic), but much of that is the prior: it
    falls to 88.09% under uniform abundances, while normalised stays at 91.75%.
  - **The aquatic reads' genera are well represented in SILVA** (14.7% of its sequences, median rank 620 of 9,118), but
    their read counts don't follow N_g (α ≈ 0.06).
  - **Naive Bayes over MEMs** with a flat prior is the strongest scheme on the dial and on the aquatic reads (about
    94–96%). Its ε matters: the held-out choice at α = 0 (10⁻⁴) falls to 84.49% at α = 1.
  - **The best genus lacks about 5–6 of a read pair's MEMs** on SILVA (from `bayesauto`'s likelihood), probably
    because of MicrobeMixer's 2% mutation and ART's errors.
  - **Space:** the default index is 1.87 GB with its RMQ (0.31 GB, unused by counting), 1.56 GB without; the paper says
    1.44.
  - **Accuracy definition:** strict, ties wrong. Reads of genera with ambiguous names count as correct (Cliffy's rule),
    which matters only for simulated reads.
- **Open:**
  - Gonzalo's run of Cliffy's approximate listing on our reads (`scripts/score_cliffy_listing.sh`, one command).
  - Whether Cliffy's `tab:others` numbers are on 5,000 or 200,000 pairs per region (`results/silva/README.md:14–25`).
  - Cliffy's index sizes on the server (`scripts/cliffy_index_size.sh`).
  - `find --classify bayes` (`edb28a0`) has no clade filter, unlike the scorer since `4f68295`. On SILVA the two can
    differ on reads whose top genus is "uncultured". They are still equal on the toy, which has none. `find` also
    classifies each mate alone, while the scorer combines mates.
  - Methods text: the dial section of `NOTES.md` has an exact description (simulation, every scheme as a formula, EM,
    `bayesauto`, the accuracy definition). It doesn't yet mention the ε → 0 limit or the scored-clades restriction.

### The KeBaB revision (in progress)
- **Experiments:** `tool/results/kebab_revision/`, `tool/results/zymo_kebab/`.
- **Zymo reruns without homopolymer compression:** prepared and checked end to end on synthetic data
  (`tool/results/zymo/RERUN_without_hpc.md`, `78b9c4d`); not yet run.
  - **Memory:** the RLCSA may grow to about 7.6 GB, so compare KATKA's `zp.csa` with `zy.csa` before starting.
  - **Comparison scripts:** `same_longest_mems.py`, `compare_longest_t.py` and `zymo_kebab/topt_lengths.py` measure
    compressed lengths only, and need a plain-length mode before the KeBaB-rule comparisons are rerun.

## Known pitfalls

- **Python 3.14** defaults to the forkserver start method. Scripts with a `Pool` need
  `multiprocessing.get_context("fork")` (done in `score_silva.py`, `simulate_dial.py`).
- **Big-BWT's multithreaded parser** (`-t`) failed an internal assertion on the synthetic pangenomes; use it
  single-threaded.
- **bigrepair:** use `-m 0` (plain RePair). Its memory-bounded RePair stalled on 64 M tag runs.
- **Timing noise:** one configuration timed at 52.3, 45.9 and 43.0 µs per read in single runs. Use repeats
  (`timing_repeats.txt`) for anything finer than about 20%.
- **c129 on the external drive:** `~/tool/chr19/c129.S` and `c129.bwt` are symlinks to it, so the drive must be
  plugged in to rebuild anything for 64 individuals.
- **`silva_queue.sh collect`** now copies only small text files (`fb7c06c`); committed results stay a few MB.
- **The cloud container can't reach arb-silva.de;** anything needing the SILVA FASTA runs on the desktop.

## Where things are

- **Tool:** `tool/README.md` (commands and presets), `tool/NOTES.md` (design, desktop commands, every experiment's
  write-up, newest sections near the top of each topic).
- **Reproduction:** `tool/CHR19.md` (*More RLBWT Tricks* v2).
- **Runners:** `tool/experiments/silva_queue.sh`, `block_scaling.sh`, `grammar_tags.sh`, `uniform_microbemixer.sh`.
- **Scorer and tables:** `tool/scripts/score_silva.py`, `silva_tables.py`, `dial_tables.py`, `silva_alpha_fit.py`.
- **Simulators:** `tool/scripts/simulate_dial.py`, `pangenome_sim.py`.
- **Results:** `tool/results/*/README.md`, and `tool/results/SILVA_DIGEST.md`.
