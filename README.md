# Two-set-graph-reduction-for-SynerGNet
Two set graph reduction is a tool to create reduced graphs from the full-sized graphs.
## Input graph
The input to this algorithm is a full-sized graph constructed by mapping gene expression, copy number variation, mutation, drug-protein association score and gene ontology terms onto the PPI network. The proteins are represented as nodes in the graph and the interactions between them are denoted as edges.

Each full-sized graph is represented by two parts: a node table and an edge table. 

The format of a node table:
| Nodes         | Feature_1 |
| ------------- | ----------|
| Protein ID 1  | x1_1      |
| Protein ID 2  | x2_1      |
```./Dataset/full_sized_graph/NodeTable_22RV1.csv``` shows the node table of cell line 22RV1 as an example. ```./Dataset/ppi_maxSubG.csv``` provides the PPI base graph.
## Graph reduction
To perform graph redcution, run ```python two_set_graph_reduction.py instances_list.csv```. The node tables, edge tables and statistics of the reduced graph can be found in ```./Dataset/ReducedGraphs/``` directory.
## Output graph




