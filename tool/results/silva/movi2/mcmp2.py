import glob,os,re
os.chdir(os.path.expanduser('~/research/katka/work/pc/movi')); tot=bad=0; cb=ct=0
for f in sorted(glob.glob('ours.*.tsv')):
    b=f[5:-4]
    O=set(tuple(l.rstrip('\n').split('\t')[:5]) for l in open(f))
    M=set()
    for l in open(f'r4/iv.{b}'):
        p=l.rstrip('\n').split('\t'); n,o=p[0].rsplit(':',1); o=int(o); M.add((n,str(o+int(p[1])),str(o+int(p[2])),p[3],p[4]))
    tot+=len(O); bad+=len(O^M); print(b,'intervals',len(O),'different',len(O^M))
    # tag answers: katka's -c listing against the bridge's, MEM by MEM
    A={}
    name=None
    for l in open(f'ourc.{b}.out'):
        if l[0]=='>': name=l[1:].strip()
        else:
            for a,e,t in re.findall(r'\[(\d+),(\d+)\] \{([^}]*)\}',l): A[(name,int(a),int(e))]=t
    B={}
    for l in open(f'r4/tags.{b}'):
        if l[0]=='>': n,o=l[1:].strip().rsplit(':',1); o=int(o)
        else:
            for a,e,t in re.findall(r'\[(\d+),(\d+)\] \{([^}]*)\}',l): B[(n,o+int(a),o+int(e))]=t
    print('   tag answers',len(A),len(B),'different',sum(1 for k in A if A[k]!=B.get(k))+len(set(B)-set(A))); ct+=len(A); cb+=sum(1 for k in A if A[k]!=B.get(k))+len(set(B)-set(A))
print('intervals',tot,'different',bad,'; tag answers',ct,'different',cb)
