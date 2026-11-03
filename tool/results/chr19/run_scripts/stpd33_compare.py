import sys
def mems(f, occ_col=None):
    d = {}
    for line in open(f):
        x = line.rstrip('\n').split('\t')
        d[(x[0], int(x[1]), int(x[2]))] = (int(x[occ_col]) if occ_col is not None and x[occ_col] not in ('-',) else None, x[-1] == 'junction')
    return d
ref = mems('ours.plain.tsv')
counts = {k: None for k in ref}
for line in open('ours.plain.tsv'):
    x = line.split('\t'); counts[(x[0], int(x[1]), int(x[2]))] = int(x[4]) - int(x[3])
print('comparison (reference: katka plain BML, %d MEMs):' % len(ref))
for name, occ in [('ours.twolevel-blocks', None), ('ours.stpd-nointervals', None), ('ours.stpd', None), ('ours.stpd-twolevel', None),
                  ('ours.stpd-locate', 5), ('theirs.bml', None), ('theirs.bml-locate', 3)]:
    d = mems(name + '.tsv', occ)
    only_t = [k for k in d if k not in ref]; only_r = [k for k in ref if k not in d]
    junc = sum(1 for k in only_t if d[k][1])
    line = '  %-24s %7d MEMs, %d only here (%d with the first occurrence across a junction), %d missing' % (name, len(d), len(only_t), junc, len(only_r))
    if occ is not None:
        diff = sum(1 for k in d if k in counts and d[k][0] != counts[k])
        line += '; occurrence counts different from ep - sp: %d' % diff
    print(line)
fbr = mems('ours.fb-twolevel-blocks.tsv'); fbt = mems('theirs.fb.tsv')
print('  forward-backward: ours %d MEMs, theirs %d, %d in both' % (len(fbr), len(fbt), len(set(fbr) & set(fbt))))
