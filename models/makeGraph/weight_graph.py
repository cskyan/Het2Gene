import networkx as nx
import pandas as pd
import pickle
G = nx.DiGraph(nx.read_gpickle("train100_weight.gpickle"))
G1 = nx.Graph(nx.read_gpickle("train100.gpickle"))
temp_dict={}
edge_dict={}
all_weight=[]

for term1,term2,term3 in list(G.edges.data()) :
    all_weight.append(term3['weight'])
    temp_dict.update({(term1,term2):int(term3['weight'])})

for t1,t2 in list(G1.edges) :
    w = temp_dict[(t1,t2)]
    edge_dict.update({(t1,t2):int(w)})

print(len(edge_dict))

#########################################################################3

hpoIndex=0
geneIndex=0
hpoMap={}
geneMap={}
hp_gene_edge =[[],[]]
hp_hp_edge = [[],[]]

for k in edge_dict :
    term1,term2=k
    w = edge_dict[k]
    if "HP" in term2 :          # hp-hp
        if term1 not in hpoMap:
            hp1 = hpoIndex
            hpoMap.update({term1: hpoIndex})
            hpoIndex += 1
        else:
            hp1 = hpoMap[term1]
        if term2 not in hpoMap:
            hp2 = hpoIndex
            hpoMap.update({term2: hpoIndex})
            hpoIndex += 1
        else:
            hp2 = hpoMap[term2]

        hp_hp_edge[0].append(hp1)
        hp_hp_edge[1].append(hp2)

    else:                       # hp-gene
        if term1 not in hpoMap:
            hi = hpoIndex
            hpoMap.update({term1: hpoIndex})
            hpoIndex += 1
        else:
            hi = hpoMap[term1]
        if term2 not in geneMap:
            gi = geneIndex
            geneMap.update({term2: geneIndex})
            geneIndex += 1
        else:
            gi = geneMap[term2]

        for i in range(w):
            hp_gene_edge[0].append(hi)
            hp_gene_edge[1].append(gi)

print(len(hp_gene_edge[0]))
print(len(hp_hp_edge[0]))

import pickle
with open("hp_gene_edge.pickle","wb+") as f :
    f.write(pickle.dumps(hp_gene_edge))

with open("hpoMap.pickle","wb+") as f :
    f.write(pickle.dumps(hpoMap))

with open("geneMap.pickle","wb+") as f :
    f.write(pickle.dumps(geneMap))

with open("hp_hp_edge.pickle","wb+") as f :
    f.write(pickle.dumps(hp_hp_edge))
