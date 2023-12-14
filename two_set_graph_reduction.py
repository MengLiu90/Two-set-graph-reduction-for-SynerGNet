import networkx as nx
import pandas as pd
import time
import argparse
import os, sys

def GraphReduction(instance_value):
    df_association = pd.read_csv('Dataset/drug_protein_association.csv')
    drug_nodes = df_association.protein.unique().tolist()

    df_ppi_nodes = pd.read_csv('Dataset/uni_protein.csv')
    ppi_nodes = df_ppi_nodes.Uni_protein.tolist()
    set_ppi_nds = set(ppi_nodes)
    intersection = set_ppi_nds.intersection(drug_nodes)  # intersection of nodes in ppi and in association
    druggable = list(intersection)
    non_druggable = list(set_ppi_nds - intersection)

    df_edge = pd.read_csv('Dataset/ppi_maxSubG.csv')  # base graph
    graph_stats = {'CellLines': [], 'Drug1': [], 'Drug2': [], 'Num_of_nodes': [], 'Num_of_edges': [],
                   'Avg_node_degree': [], 'Graph_density': [], 'Graph_diameter': []}

    start = time.time()
    # read the instance
    sample = instance_value
    string = sample.replace('\r', '')
    cell_drug_list = string.split('_')
    print('instance being processed', cell_drug_list)
    cell = cell_drug_list[0]
    drug1 = cell_drug_list[1]
    drug2 = cell_drug_list[2]

    NT_1 = pd.read_csv(f'Dataset/full_sized_graph/NodeTable_{cell}_part1.csv')
    NT_2 = pd.read_csv(f'Dataset/full_sized_graph/NodeTable_{cell}_part2.csv')
    NT_3 = pd.read_csv(f'Dataset/full_sized_graph/NodeTable_{cell}_part3.csv')
    df_nodetable = pd.concat([NT_1, NT_2, NT_3])
    df_nodetable.loc[df_nodetable.Nodes.isin(druggable), "druggable"] = 1
    df_nodetable.loc[df_nodetable.Nodes.isin(non_druggable), "druggable"] = 0

    table_col = df_nodetable.columns.tolist()
    col_len = len(table_col)

    assoc_drug1_target = df_association.loc[df_association.chemical == drug1]
    assoc_drug2_target = df_association.loc[df_association.chemical == drug2]

    assoc_target_part = df_association.loc[df_association.chemical.isin([drug1, drug2])]
    target_nodes = assoc_target_part.protein.unique().tolist()  # target nodes in this instance

    EdgeTable = df_edge.copy()
    EdgeTable['gex1'] = EdgeTable['Prot1'].map(df_nodetable.set_index('Nodes')['gex'])
    EdgeTable['gex2'] = EdgeTable['Prot2'].map(df_nodetable.set_index('Nodes')['gex'])
    EdgeTable['druggable1'] = EdgeTable['Prot1'].map(df_nodetable.set_index('Nodes')['druggable'])
    EdgeTable['druggable2'] = EdgeTable['Prot2'].map(df_nodetable.set_index('Nodes')['druggable'])

    # keepers isolate the target nodes and separate drugable/nondrugable, different gex
    keepers = EdgeTable[(EdgeTable['Prot1'].isin(target_nodes)) |
                        (EdgeTable['Prot2'].isin(target_nodes)) |
                        (EdgeTable['gex1'] != EdgeTable['gex2']) |
                        (EdgeTable['druggable1'] != EdgeTable['druggable2'])]


    temp_G = nx.from_pandas_edgelist(EdgeTable, 'Prot1', 'Prot2', 'Score')
    keepers_tuple = [(x, y) for x, y in zip(keepers['Prot1'].tolist(), keepers['Prot2'].tolist())]
    temp_G.remove_edges_from(keepers_tuple)
    Induced_subG = [temp_G.subgraph(c).copy() for c in nx.connected_components(temp_G)]
    connected_subG = [x for x in Induced_subG if len(x.nodes) > 1]

    NodeTable = df_nodetable.copy()
    d = {i: i for i in NodeTable['Nodes'].tolist()}
    ids = 0
    for G in connected_subG:
        nodes = list(G.nodes)
        node_attrs = NodeTable[NodeTable.Nodes.isin(nodes)]
        feature_sum = node_attrs.sum(axis=0)
        ### merge nodes
        for i in range(col_len - 2):
            col = table_col[i + 2]
            NodeTable.loc[NodeTable.Nodes.isin(nodes), col] = feature_sum[col]
        NodeTable.loc[NodeTable.Nodes.isin(nodes), 'Nodes'] = f'v_{ids}'

        for n in nodes:
            d[n] = f'v_{ids}'

        EdgeTable.loc[EdgeTable.Prot1.isin(nodes), 'Prot1'] = f'v_{ids}'
        EdgeTable.loc[EdgeTable.Prot2.isin(nodes), 'Prot2'] = f'v_{ids}'

        ids += 1

    EdgeTable = EdgeTable[EdgeTable.Prot1 != EdgeTable.Prot2]
    NodeTable.drop_duplicates('Nodes', inplace=True)

    ### merge edges
    grouped = EdgeTable.groupby(['Prot1', 'Prot2'], as_index=False)
    EdgeTable_1 = grouped.agg({'Score': 'median'})

    ### remove A-B/ B-A duplication in the edges
    EdgeTable_1['check_string'] = EdgeTable_1.apply(
        lambda row: ''.join(sorted([str(row['Prot1']), str(row['Prot2'])])),
        axis=1)
    EdgeTable_1.drop_duplicates('check_string', inplace=True)
    EdgeTable_1.drop('check_string', inplace=True, axis=1)

    ### check if all nodes are connected after reduction
    uni_node_nodeTable = NodeTable.Nodes.tolist()
    protein1 = EdgeTable_1.Prot1
    protein2 = EdgeTable_1.Prot2
    protein = pd.concat([protein1, protein2])
    uni_node_edgeTable = protein.unique().tolist()

    # recalculate mutation and cnv of nodes obtained by merging
    for j in range(0, 16):
        col = table_col[j + 2]
        NodeTable[col] = NodeTable[col].apply(lambda x: 1 if x > 0 else 0)

    # # assign affinity score to target nodes
    NodeTable['Affinity1'] = NodeTable['Nodes'].map(assoc_drug1_target.set_index('protein')['combined_score'])
    NodeTable['Affinity1'] = NodeTable['Affinity1'].fillna(0)
    NodeTable['Affinity2'] = NodeTable['Nodes'].map(assoc_drug2_target.set_index('protein')['combined_score'])
    NodeTable['Affinity2'] = NodeTable['Affinity2'].fillna(0)
    ## get the max(Affinity_1, Affinity_2) as the Affinity score of each node
    NodeTable['Affinity'] = NodeTable[['Affinity1', 'Affinity2']].max(axis=1)
    NodeTable.drop(columns=['Affinity1', 'Affinity2'], inplace=True)
    NodeTable.drop(columns=['druggable'], inplace=True)

    nodetable_dir = 'Dataset/ReducedGraphs/NodeTables/'
    if not os.path.isdir(nodetable_dir):
        os.mkdir(nodetable_dir)
    edgetable_dir = 'Dataset/ReducedGraphs/EdgeTables/'
    if not os.path.isdir(edgetable_dir):
        os.mkdir(edgetable_dir)
    noderename_dir = 'Dataset/ReducedGraphs/NodeRename/'
    if not os.path.isdir(noderename_dir):
        os.mkdir(noderename_dir)
    nodedegree_dir = 'Dataset/ReducedGraphs/NodeDegree/'
    if not os.path.isdir(nodedegree_dir):
        os.mkdir(nodedegree_dir)
    statistics_dir = 'Dataset/ReducedGraphs/statistics/'
    if not os.path.isdir(statistics_dir):
        os.mkdir(statistics_dir)

    # Get the column names
    columns = NodeTable.columns.tolist()
    # Move the last column to the 13th position
    last_column = columns.pop(-1)  # Remove the last column
    columns.insert(18, last_column)  # Insert the last column at the 13th position
    # Reorder the DataFrame columns
    NodeTable = NodeTable[columns]

    NodeTable.to_csv(nodetable_dir + f'NodeTable_{cell}_{drug1}_{drug2}.csv', index=False)
    EdgeTable_1.to_csv(edgetable_dir + f'EdgeTable_{cell}_{drug1}_{drug2}.csv', index=False)
    elapsed_time_fl = (time.time() - start)

    ## record the name-change on the merged nodes
    d1 = {'original_node': [], 'new_node': []}
    for k, v in d.items():
        d1['original_node'].append(k)
        d1['new_node'].append(v)
    df_node_rename = pd.DataFrame(d1)
    df_node_rename.to_csv(noderename_dir + f'NodesIn_{cell}_{drug1}_{drug2}.csv', index=False)

    ## graph statistics for each instance
    ## since the graph reduction based on target nodes of drugs not druggable nodes
    ## therefore, for each synergy instance, the topology is different
    # Create new graph (reduced graph)
    ReduG = nx.from_pandas_edgelist(EdgeTable_1, 'Prot1', 'Prot2', 'Score')
    # Graph statistics
    num_nodes = len(ReduG.nodes)
    num_edges = len(ReduG.edges)

    Nodes_degree = []
    nodes_proc = []
    for node in ReduG.nodes:
        degree = ReduG.degree[node]
        nodes_proc.append(node)
        Nodes_degree.append(degree)
    degree_tuples = list(zip(nodes_proc, Nodes_degree))
    degreeCheck = pd.DataFrame(degree_tuples, columns=['Nodes', 'Degree'])
    degreeCheck.to_csv(nodedegree_dir + f'NodeDegree_{cell}_{drug1}_{drug2}.csv', index=False)

    total_degree = sum(dict(ReduG.degree).values())
    avg_degree = total_degree / num_nodes

    density = 2 * num_edges / (num_nodes * (num_nodes - 1))

    ## diameter of graph is the maximum distance among all pairs of vertices
    diamtr = nx.diameter(ReduG)

    ## keep a record of the graph statistics
    graph_stats['CellLines'].append(cell)
    graph_stats['Drug1'].append(drug1)
    graph_stats['Drug2'].append(drug2)
    graph_stats['Num_of_nodes'].append(num_nodes)
    graph_stats['Num_of_edges'].append(num_edges)
    graph_stats['Avg_node_degree'].append(avg_degree)
    graph_stats['Graph_density'].append(density)
    graph_stats['Graph_diameter'].append(diamtr)

    df_graph_stats = pd.DataFrame.from_dict(graph_stats)
    df_graph_stats.to_csv(statistics_dir + f'GraphStatistics_{cell}_{drug1}_{drug2}.csv', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python two_set_graph_reduction_original_data.py input_file")
    else:
        input_file_path = sys.argv[1]

        with open(input_file_path, 'r') as file:
            for line in file:
                instance_value = line.strip()
                GraphReduction(instance_value)









