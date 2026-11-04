Subject: KeBaB: memory leak in process_sequences with more than one thread

Hi,

While rerunning the KeBaB paper's Section 5.1 setting (the SPUMONI 2 Zymo mock community), I found a memory leak in
`kebab build` and `kebab scan` when they run with more than one thread. The version is commit 4bed76f.

What happens. In `process_sequences` (src/kebab.cpp, about lines 46–80), when `threads > 1`, every sequence read is
copied with `strdup` (both `seq->seq.s` and `seq->name.s`). The copies are freed only once, after each thread's
`while` loop ends. So every sequence except the last one a thread handles is never freed, and memory grows with the
total size of the input. With one thread there is no copy (kseq's own buffer is used), so there is no leak.

Effect. `kebab build -k 20 -e 0.1 -f 1 -t 12` on the 27.2 GB uncompressed Zymo reference reached 11.5 GB resident
43% of the way through "Estimating Cardinality", and was killed by an 11 GB memory limit. The resident size tracked
the bytes read almost exactly. `kebab scan` calls the same function, so with several threads on the 14.8 Gbases of
positive reads it would hold all of them in memory too. On a machine with plenty of RAM this goes unnoticed. The
paper's machine had 270 GB.

Workaround. `-t 1`. `kebab build -t 1` on the same input ran in about 10 MB, and the filter is the same: the same
k-mers are inserted, and HyperLogLog gives the same estimate.

Suggested fix. Free the copies after each sequence is processed, inside the loop:

```cpp
            process_func(seq_info);
            if (threads > 1) {
                free(seq_info.seq_content);
                free(seq_info.seq_name);
            }
        }   // end while
```

and drop the `free` calls after the loop. Otherwise the last sequence would be freed twice, or a stale pointer
freed when the final `kseq_read` returns -1. Equivalently, set both pointers to `nullptr` after freeing them inside
the loop and keep the final cleanup.

Best,
Travis
