# Generates the TikZ for the de Bruijn graph figure (draft; to be redrawn in Ipe).
from search import graph, kmers
F,U,V,Vp,G='GTA','AACGCTC','GTCAAGG','GTCCAGG','TGG'
S=[F+U+V+G]*5+[F+U+Vp+G]+[F+V+U+G]*5+[F+Vp+U+G]
Gr=graph(S)
mult={}
for s in S:
    km=kmers(s)
    for a,b in zip(km,km[1:]): mult[(a,b)]=mult.get((a,b),0)+1
br={'GTA':(0,0),'AAC':(3,1.6),'CTC':(7,1.6),'GTC':(3,-1.6),'AGG':(7,-1.6),'TGG':(10,0)}
# chains between branching k-mers: (start, internal..., end, perpendicular offset)
chains=[('GTA',['TAA','AAA'],'AAC',0),('GTA',['TAG','AGT'],'GTC',0),
        ('AAC',['ACG','CGC','GCT'],'CTC',0),('CTC',['TCG','CGT'],'GTC',0),
        ('CTC',['TCT','CTG'],'TGG',0),('GTC',['TCA','CAA','AAG'],'AGG',0),
        ('GTC',['TCC','CCA','CAG'],'AGG',-1.0),('AGG',['GGA','GAA'],'AAC',0),
        ('AGG',['GGT','GTG'],'TGG',0)]
pos=dict(br)
for a,mid,b,off in chains:
    (x0,y0),(x1,y1)=br[a],br[b]
    dx,dy=x1-x0,y1-y0; L=(dx*dx+dy*dy)**.5; px,py=-dy/L,dx/L
    for i,m in enumerate(mid):
        t=(i+1)/(len(mid)+1)
        pos[m]=(x0+t*dx+off*px*(1 if 0<t<1 else 0), y0+t*dy+off*py)
assert set(pos)==set(Gr.nodes()), set(Gr.nodes())^set(pos)
out=[]
out.append(r'\begin{tikzpicture}[x=0.62cm,y=0.62cm,>={Stealth[length=1.4mm]},v/.style={draw,rounded corners=1pt,inner sep=1pt,font=\tiny\ttfamily}]')
for k,(x,y) in pos.items(): out.append(rf'\node[v{",fill=gray!20" if k in br else ""}] ({k}) at ({x:.2f},{y:.2f}) {{{k}}};')
for (a,b) in Gr.edges(): out.append(rf'\draw[->] ({a}) -- ({b});')
out.append(r'\end{tikzpicture}')
open('dbg_full.tex','w').write('\n'.join(out)+'\n')
# compacted
cn={'GTA':(0,0),'TAAA':(1.5,0.8),'TAGT':(1.5,-0.8),'AACGCTC':(5,1.6),'TCGT':(5.6,0.45),'GGAA':(4.4,0.45),
    'GTC':(3,-1.6),'TCAAG':(5,-1.6),'TCCAG':(5,-2.6),'AGG':(7,-1.6),'TCTG':(8.5,0.8),'GGTG':(8.5,-0.8),'TGG':(10,0)}
ce=[('GTA','TAAA'),('GTA','TAGT'),('TAAA','AACGCTC'),('TAGT','GTC'),('AACGCTC','TCGT'),('AACGCTC','TCTG'),
    ('TCGT','GTC'),('GTC','TCAAG'),('GTC','TCCAG'),('TCAAG','AGG'),('TCCAG','AGG'),('AGG','GGAA'),('AGG','GGTG'),
    ('GGAA','AACGCTC'),('TCTG','TGG'),('GGTG','TGG')]
w={('GTC','TCAAG'):10,('TCAAG','AGG'):10,('GTC','TCCAG'):2,('TCCAG','AGG'):2}
out=[r'\begin{tikzpicture}[x=0.62cm,y=0.62cm,>={Stealth[length=1.4mm]},v/.style={draw,rounded corners=1pt,inner sep=1.5pt,font=\scriptsize\ttfamily}]']
for k,(x,y) in cn.items():
    style='v,fill=gray!20' if k in ('AACGCTC',) else ('v,fill=gray!8' if k in ('GTC','TCAAG','TCCAG','AGG') else 'v')
    out.append(rf'\node[{style}] ({k}) at ({x},{y}) {{{k}}};')
for a,b in ce:
    lab=w.get((a,b),'')
    out.append(rf'\draw[->] ({a}) -- node[font=\tiny,fill=white,inner sep=0.5pt] {{{lab}}} ({b});' if lab else rf'\draw[->] ({a}) -- ({b});')
out.append(r'\end{tikzpicture}')
open('dbg_compact.tex','w').write('\n'.join(out)+'\n')
print('ok')
