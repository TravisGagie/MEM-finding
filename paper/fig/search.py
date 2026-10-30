import random, networkx as nx
random.seed(7)
k=3
def kmers(s): return [s[i:i+k] for i in range(len(s)-k+1)]
def graph(strings):
    G=nx.DiGraph()
    for s in strings:
        km=kmers(s)
        for a in km: G.add_node(a)
        for a,b in zip(km,km[1:]): G.add_edge(a,b)
    return G
def compact(G):
    internal=lambda v: G.in_degree(v)==1 and G.out_degree(v)==1
    out=[]
    for v in G:
        if internal(v): continue
        if G.out_degree(v)==0: out.append([v])
        for w in G.successors(v):
            p=[v,w]
            while internal(w):
                w=next(G.successors(w)); p.append(w)
            out.append(p)
    return out
