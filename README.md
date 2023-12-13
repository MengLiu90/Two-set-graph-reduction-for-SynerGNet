# Two-set-graph-reduction-for-SynerGNet
Two set graph reduction is a tool to create reduced graphs from the full-sized graphs.
## Input graph
The input to this algorithm is a full-sized graph constructed by mapping gene expression, copy number variation, mutation and drug-protein association score onto the PPI network. 

This creates a node table for each cell line. The PPI serves as the edge table for all the full-sized graphs. ```./Dataset/full_sized_graph/NodeTable_22RV1.csv``` shows the node table of cell line 22RV1 as an example.
## Graph reduction
To perform graph redcution, run ```python two_set_graph_reduction.py instances_list.csv```. The node tables, edge tables and statistics of the reduced graph will be saved in ```./Dataset/ReducedGraphs/``` directory.




