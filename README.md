# Two-set-graph-reduction-for-SynerGNet
Two set graph reduction is a tool to create reduced graphs from the full-sized graphs.
## Input graph
The input to this algorithm is a full-sized graph constructed by mapping gene expression, copy number variation, mutation and drug-protein association score onto the PPI network. (Gene ontology (GO) term is saved at this step as it was not used during the reduction process, and the dimension of this feature is 200, which would increase the processing burden. GO term was assigned to the reduced graph for training purpose.) This creates a node table for each cell line. The PPI serves as the edge table for all the full-sized graphs.

```NodeTable_22RV1.csv``` shows the node table of cell line 22RV1 as an example.


