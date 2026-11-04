| configuration | time (s) | peak (GB) | files (MB) |
|---|---|---|---|
| plain backward search (RLCSA) | 90.79 / 93.88 | 0.89 | 642 |
| two-level with blocks | 13.05 / 18.01 | 1.19 | 780 |
|   with locating | 17.98 / 30.57 | 1.73 | 1320 |
| STPD over the text | 7.03 / 12.82 | 1.16 | 903 |
|   with locating | 7.01 / 13.68 | 1.16 | 903 |
| STPD over the text, read through the parse | 45.83 / 90.04 | 1.01 | 753 |
| STPD over the parse, plain parse | 11.77 / 20.48 | 1.84 | 1022 |
|   with locating | 11.41 / 22.40 | 1.84 | 1315 |
| STPD over the parse, RLZ parse | 15.07 / 26.74 | 1.80 | 979 |
| published STPD-index (RLZ text) | (set STPD to run it) | | |

1 individuals. Search time for L = 1000 / L = 100 (excluding loading), the peak memory of the search at L = 1000, and the size of what the search reads. Every configuration gave the same MEMs as plain search.
