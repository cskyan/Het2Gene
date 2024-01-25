import pickle

import torch_geometric
import random
from scipy import spatial
import torch
import torch_geometric.nn as nn
import numpy as np

seed = 5
# fix random seed
def same_seeds(seed):
    torch.manual_seed(seed)  # （CPU）
    if torch.cuda.is_available():  # （GPU)
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.deterministic = True
    print("set")
same_seeds(seed)

def TransPose(ls:list):
    ls_temp=ls.copy()
    temp=ls_temp[0].copy()
    ls_temp[0]=ls_temp[1].copy()
    ls_temp[1]=temp
    return ls_temp


def get_hpo_gene_edges() :
    with open("hp_gene_edge.pickle", "rb+") as f:
        hpo_gene_edges = list(pickle.load(f))
    return hpo_gene_edges
import pickle
def get_hpo_hpo_edge():
    with open("hp_hp_edge.pickle","rb+") as f :
        hpo_hpo_edges = list(pickle.load(f))
    return hpo_hpo_edges


hpo_hpo_edges=get_hpo_hpo_edge()   #  cadaNet
hp_gene_edge=get_hpo_gene_edges()
import random
random.seed(5)

dict_edge_index={("hpo","resultFrom","gene"):torch.LongTensor(hp_gene_edge),
                 ("gene","resultTo","hpo"):torch.LongTensor(TransPose(hp_gene_edge)),
                 ("hpo","hhc","hpo"):torch.LongTensor(hpo_hpo_edges),
                 ("hpo","hhc","hpo"):torch.LongTensor(TransPose(hpo_hpo_edges))}
metapath=[

("gene","resultTo","hpo"),
("hpo","hhc","hpo"),
("hpo","hhc","hpo"),
("hpo","resultFrom","gene"),


]
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
model=nn.MetaPath2Vec(dict_edge_index,embedding_dim=300,metapath=metapath,walk_length=200,context_size=3,sparse=False,walks_per_node=15,num_negative_samples=3).to(device)
loader = model.loader(batch_size=256, shuffle=True)
optimizer = torch.optim.Adam(list(model.parameters()), lr=0.1)
StepLR = torch.optim.lr_scheduler.StepLR(optimizer, step_size=1, gamma=0.9)
print("训练开始")
def getEembedding(epoch, log_steps=100 ,eval_steps=20) :
    model.train()
    total_loss = 0
    for i, (pos_rw, neg_rw) in enumerate(loader):
        optimizer.zero_grad()
        loss = model.loss(pos_rw.to(device),neg_rw.to(device))
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
        if (i + 1) % log_steps == 0:
            print((f'Epoch: {epoch}, Step: {i + 1:05d}/{len(loader)}, '
                   f'Loss: {total_loss / log_steps:.4f}'),end="\t")
            #print("")
            print('scheduler: {}'.format(StepLR.get_lr()[0]))
            StepLR.step()
            total_loss = 0




for i in range(1,70) :
    try:
        getEembedding(i,10)
    except KeyboardInterrupt:
        print("break")
        break

embedding_hpo = model.forward("hpo").tolist()
embedding_gene = model.forward("gene").tolist()
with open("embedding_hpo.pickle", "wb+") as f:
    f.write(pickle.dumps(embedding_hpo))


with open("embedding_gene.pickle", "wb+") as f:
    f.write(pickle.dumps(embedding_gene))


print("saved")