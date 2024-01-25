# Het2Gene
a phenotype-driven model for gene prioritization.The heterogeneous graph embedding algorithm is used to learn the embedding representation of heterogeneous graph nodes, and the score of candidate causal genes is calculated according to the embedding, so as to prioritize.

# Folder description
The data used in the training model are placed in `./data/` .The edge relation, weight and trained embedding information have been encapsulated by pickle module. In addition, the test data used are also included. In `./models/`, it includes the method of constructing weighted and unweighted graphs, node coding method, training method, test code, etc

# Usage
In folder `./models/prioritize/Het2Gene/` ,run the following command(The same goes for weighted graphs) can use Het2Gene:
```
python het2gene.py --hps [parameter] --out_dir [-options][parameter] --topn [-options][parameter]
```
**hps**: Abnormal phenotype set from patients，required

**out_dir**: Directory for outputting results，default `RankResult`

**topn**: Output candidate genes in the top n,default `1000`

An example command:
```
python het2gene.py --hps HP:0000573,HP:0001102,HP:0003115,HP:0001681,HP:0008067,HP:0004417 --out_dir myResult --topn 5
```

Output Example：
```
Rank	EntrezID	Symbol	Score
1	368	ABCC6	11.56
2	64132	XYLT2	1.983
3	5167	ENPP1	-0.815
4	64131	XYLT1	-0.936
5	4000	LMNA	-3.271

```
# Web Resources for Comparison Methods
CADA:[https://cada.gene-talk.de/webservice](https://cada.gene-talk.de/webservice)

AMELIE:[https://amelie.stanford.edu](https://amelie.stanford.edu)

Phen2Gene:[https://phen2gene.wglab.org](https://phen2gene.wglab.org)

GADO:[https://genenetwork.nl/gado](https://genenetwork.nl/gado)

# License
See the LICENSE file for license rights and limitations



