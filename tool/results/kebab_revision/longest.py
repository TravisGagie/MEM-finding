#!/usr/bin/env python3
# longest.py mode ...: per-read t longest MEMs from ropebwt3 mem output (name, start, end, count; 0-based, end
# exclusive), with fragment names read:START-END (1-based inclusive) mapped back to read coordinates.
#   longest.py fallback mems lambda.tsv t > reads_to_rerun.txt      (b: reads with T+1 fragments and fewer than t
#                                                                       MEMs at least as long as the (T+1)st)
#   longest.py compare run.mem[,fallback.mem] ref.mem reads.fq t   (reads whose t longest lengths / spans differ;
#                                                                       streams run.mem and ref.mem in read order)
# .zst inputs are read through zstdcat. MEMs are deduplicated by read span (a MEM inside two overlapping fragments).
import sys, subprocess, collections
def open_any(f):
    if f.endswith('.zst'): return subprocess.Popen(['zstdcat', f], stdout=subprocess.PIPE, text=True).stdout
    return open(f)
def mems(files):
    d = collections.defaultdict(set)
    for f in files:
        for line in open_any(f):
            x = line.split('\t')
            if len(x) < 3: continue
            name, s, e = x[0], int(x[1]), int(x[2])
            if ':' in name and '-' in name.rsplit(':', 1)[1]:
                name, rng = name.rsplit(':', 1); a = int(rng.split('-')[0]) - 1; s += a; e += a
            d[name].add((s, e))
    return d
def topt(spans, t):     # the t longest with ties at the t-th place, as (sorted lengths of the first t, set of spans)
    L = sorted((e - s for s, e in spans), reverse=True)
    if not L: return (), frozenset()
    cut = L[min(t, len(L)) - 1]
    return tuple(L[:t]), frozenset(x for x in spans if x[1] - x[0] >= cut)
mode = sys.argv[1]
if mode == 'fallback':
    # streaming: lambda.tsv and the MEMs of the searched fragments are both in read order
    t = int(sys.argv[4]); n = rerun = few = 0
    g = None
    def grp(f):
        cur, sp = None, []
        for line in open_any(f):
            x = line.split('\t')
            if len(x) < 3: continue
            name = x[0].rsplit(':', 1)[0] if ':' in x[0] else x[0]
            if name != cur:
                if cur is not None: yield cur, sp
                cur, sp = name, []
            sp.append(int(x[2]) - int(x[1]))
        if cur is not None: yield cur, sp
    g = grp(sys.argv[2]); nxt = next(g, None)
    for line in open(sys.argv[3]):
        name, cnt, lt = line.split('\t'); cnt, lt = int(cnt), int(lt); n += 1
        lens = []
        if nxt is not None and nxt[0] == name: lens = nxt[1]; nxt = next(g, None)
        if lt == 0: few += 1; continue                       # at most T fragments: every MEM >= L was searched
        if sum(1 for l in lens if l >= lt) < t: print(name); rerun += 1
    if nxt is not None: sys.exit('MEMs not in the fragments order')
    print('reads %d: %d with at most T fragments (no check needed), %d fall back (%.3f%%)' % (n, few, rerun, 100 * rerun / n), file=sys.stderr)
else:
    # streaming: ropebwt3 (one thread) and kebab scan write reads in input order, so each file is read once, a read
    # at a time; the fallback file (b) is small and is held in memory
    t = int(sys.argv[5])
    def groups(f):     # (read name, set of read spans) for consecutive lines of the same read
        cur, sp = None, set()
        for line in open_any(f):
            x = line.split('\t')
            if len(x) < 3: continue
            name, s_, e_ = x[0], int(x[1]), int(x[2])
            if ':' in name and '-' in name.rsplit(':', 1)[1]:
                name, rng = name.rsplit(':', 1); a = int(rng.split('-')[0]) - 1; s_ += a; e_ += a
            if name != cur:
                if cur is not None: yield cur, sp
                cur, sp = name, set()
            sp.add((s_, e_))
        if cur is not None: yield cur, sp
    class cursor:
        def __init__(self, f): self.g = groups(f); self.nxt = next(self.g, None)
        def take(self, name):
            if self.nxt is not None and self.nxt[0] == name:
                v = self.nxt[1]; self.nxt = next(self.g, None); return v
            return set()
    runs = sys.argv[2].split(',')
    main = cursor(runs[0]); extra = mems(runs[1:]) if len(runs) > 1 else {}
    ref = cursor(sys.argv[3])
    n = ld = sd = none_ref = none_run = 0
    with open(sys.argv[4]) as F:
        for i, line in enumerate(F):
            if i % 4: continue
            nm = line[1:].split()[0]; n += 1
            a, b = topt(main.take(nm) | extra.get(nm, set()), t), topt(ref.take(nm), t)
            if not b[0]: none_ref += 1
            if not a[0] and b[0]: none_run += 1
            if a[0] != b[0]: ld += 1
            if a[1] != b[1]: sd += 1
    if main.nxt is not None or ref.nxt is not None: sys.exit('listings not in read order (or names not in reads.fq)')
    print('reads %d, t = %d: the t longest lengths differ in %d (%.3f%%), spans (with ties) differ in %d (%.3f%%); '
          'reads with no MEM >= L: %d in the reference, %d more in this run' % (n, t, ld, 100 * ld / n, sd, 100 * sd / n, none_ref, none_run))
