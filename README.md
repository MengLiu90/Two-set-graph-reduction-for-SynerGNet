# Two-set-graph-reduction-for-SynerGNet
Two set graph reduction is a tool to create reduced graphs from the full-sized graphs.
## Dependencies
1. Networkx 2.7.1
2. Pandas 1.1.3
## Input graph
The input to this algorithm is a full-sized graph constructed by mapping gene expression, copy number variation, mutation, drug-protein association score and gene ontology terms onto the PPI network. The proteins are represented as nodes in the graph and the interactions between them are denoted as edges.

Each full-sized graph is represented by two parts: a node table and an edge table. 

The format of a node table:
| Nodes         | Feature_1 | Feature_2 | ... |Feature_m |
| ------------- | ----------|---------- | ----|----------|
| UNiprot ID 1  | x1_1      |x1_2       | ... |x1_m      |
| UNiprot ID 2  | x2_1      |x2_2       | ... |x2_m      |
| ...           | ...       |...        | ... |...       |
| UNiprot ID n  | xn_1      |xn_2       | ... |xn_m      |

The format of an edge table:
| Node1 ID   | Node2 ID  | Edge score | 
| ---------- | ----------|----------  | 
| UNiprot 1  | UNiprot 2 |xxx         | 
| UNiprot 1  | UNiprot 3 |xxx         | 
| ...        | ...       |...         | 
| UNiprot n  | UNiprot k |xxx         |

```./Dataset/full_sized_graph/NodeTable_22RV1.csv``` shows the node table of cell line 22RV1 as an example. ```./Dataset/ppi_maxSubG.csv``` provides the PPI base graph, which serves as the edge table for the full-sized graphs.
## Graph reduction
To perform graph redcution, run ```python two_set_graph_reduction.py Dataset/instances_list.csv```. 
## Output graph
The node tables, edge tables and statistics of the reduced graph can be found in ```./Dataset/ReducedGraphs/``` directory. Each reduced graph is represented by its node table and edge table. 




