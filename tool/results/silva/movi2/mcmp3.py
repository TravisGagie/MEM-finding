import re,sys
ref,got,name,f=sys.argv[1:5]
A={}
for l in open(ref):
    if l[0]=='>': n=l[1:].strip()
    else:
        for a,e,t in re.findall(r'\[(\d+),(\d+)\] \{([^}]*)\}',l): A[(n,int(a),int(e))]=t
d=0;c=0
for l in open(got):
    if l[0]=='>': n,o=l[1:].strip().rsplit(':',1); o=int(o)
    else:
        for a,e,t in re.findall(r'\[(\d+),(\d+)\] \{([^}]*)\}',l):
            c+=1
            if A.get((n,o+int(a),o+int(e)))!=t: d+=1
print('%s %s: %d of %d MEMs; different %d'%(name,f,c,len(A),d+abs(len(A)-c)) + ('' if d==0 and c==len(A) else ' DIFFERS'),flush=True)
