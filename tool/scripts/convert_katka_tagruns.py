# KATKA tag runs (n u64, r u64, start u64[r], tag u16[r]) -> the new tool's (KTRUNS01, n, r, start u64[r], r, tag u32[r])
import sys, os, numpy as np
src, dst = sys.argv[1], sys.argv[2]
with open(src, 'rb') as f: n, r = np.frombuffer(f.read(16), dtype='<u8')
n, r = int(n), int(r)
w = (os.path.getsize(src) - 16) // r - 2 if r else 8      # KATKA's start width: 4 bytes if n < 2^32, else 8
if w not in (4, 8) or os.path.getsize(src) != 16 + (w + 2) * r: sys.exit(f"{src}: size {os.path.getsize(src)} fits no layout (r = {r})")
st = np.memmap(src, dtype='<u%d' % w, mode='r', offset=16, shape=(r,))
tg = np.memmap(src, dtype='<u2', mode='r', offset=16 + w * r, shape=(r,))
C = 1 << 26
with open(dst + '.part', 'wb') as o:
    o.write(b'KTRUNS01'); o.write(np.array([n, r], dtype='<u8').tobytes())
    prev = -1
    for i in range(0, r, C):
        a = np.asarray(st[i:i + C])
        if a[0] <= prev or (len(a) > 1 and not (np.diff(a.astype(np.int64)) > 0).all()): sys.exit("starts not increasing")
        prev = int(a[-1]); o.write(a.astype('<u8').tobytes())
    if r and (st[0] != 0 or prev >= n): sys.exit("starts out of range")
    o.write(np.array([r], dtype='<u8').tobytes())
    for i in range(0, r, C):
        b = np.asarray(tg[i:i + C])
        if len(b) > 1 and not (b[1:] != b[:-1]).all(): sys.exit("two adjacent runs with the same tag")
        o.write(b.astype('<u4').tobytes())
os.rename(dst + '.part', dst)
print(f"n = {n}, runs = {r}, start width {w} bytes, tags used = {len(np.unique(tg))}, wrote {os.path.getsize(dst)} bytes")
