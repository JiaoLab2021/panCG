import os
import pandas as pd
import networkx as nx
import pickle
import copy
import io


class BulitPhySynNet:
    def __init__(self, gene_index_workDir, graph_file, speciesList):
        self.gene_index_workDir = gene_index_workDir
        self.graph_file = graph_file
        self.speciesList = speciesList

    @staticmethod
    def read_net(file_):
        if not os.path.exists(file_):
            # 文件不存在，返回空图
            return nx.Graph()
        sequenceIO = io.StringIO()
        n = 0
        with open(file_, "r") as f:
            for line in f:
                if line.startswith("#"):
                    n += 1
                else:
                    line_li = line.strip().split()
                    line_li.append(str(n))
                    sequenceIO.write("\t".join(line_li) + "\n")
        sequenceIO.seek(0)
        try:
            df = pd.read_csv(sequenceIO, sep="\t", header=None, names=["source", "target", "synteny_weight", "block"])
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=["source", "target", "synteny_weight", "block"])
        # df = pd.read_csv(sequenceIO, sep="\t", header=None, names=["source", "target", "synteny_weight", "block"])
        sequenceIO.close()
        df['synteny_weight'] = df['synteny_weight'].astype('str')
        df["synteny_weight"] = df["synteny_weight"].apply(lambda x: x[:-1] if x.endswith('L') else x)
        df["synteny_weight"] = df["synteny_weight"].astype(int)
        df = df.loc[df.groupby(['source', 'target'])['synteny_weight'].idxmax()]
        net_G = nx.from_pandas_edgelist(
            df,  # Data source DataFrame
            'source',  # Source Node Column
            'target',  # Target Node Column
            ['synteny_weight', 'block']  # Edge Weight
        )
        return net_G

    @staticmethod
    def overlap_net(anchors_G, lifted_anchors_G):
        """
        Extract the nodes that do not exist in anchors_G but exist in lifted_anchors_G,
        and then extract the edges including the nodes in lifted_anchors_G and add them to anchors_G
        Returns:
            out_graph: Deep copy anchors_G and add the egde network
        """
        out_graph = copy.deepcopy(anchors_G)
        other_nodes = set(lifted_anchors_G.nodes) - set(anchors_G.nodes)
        edges_to_add = []
        for node in other_nodes:
            block_max_edges = {}
            adjacent_nodes = lifted_anchors_G.adj[node]
            for neighbor, attributes in adjacent_nodes.items():
                block = attributes['block']
                if block not in block_max_edges:
                    block_max_edges[block] = {
                        'node1': node,
                        'node2': neighbor,
                        'synteny_weight': attributes['synteny_weight']
                    }
                else:
                    if attributes['synteny_weight'] > block_max_edges[block]['synteny_weight']:
                        block_max_edges[block] = {
                            'node1': node,
                            'node2': neighbor,
                            'synteny_weight': attributes['synteny_weight']
                        }
            for block, info in block_max_edges.items():
                edges_to_add.append((info["node1"], info["node2"], {"synteny_weight": info["synteny_weight"], "block": block}))
        out_graph.add_edges_from(edges_to_add)
        return out_graph

    @staticmethod
    def find_subgraph(G, nodes: list):
        """
        Extract the subnet of a given node in the network
        """
        subgraph = G.subgraph(nodes).copy()
        return subgraph

    @staticmethod
    def G2dataframe(graph, data=False):
        """
        Output the network as a data frame
        """
        if data:
            edges_df = pd.DataFrame(list(graph.edges(data=True)), columns=["source", "target", "attributes"])
        else:
            edges_df = pd.DataFrame(list(graph.edges()), columns=["source", "target"])
        return edges_df

    @staticmethod
    def get_gene_from_index(df, speciesList):
        result = df[speciesList].values.tolist()
        li = [item for sublist in result for j in sublist if j != "." for item in j.split(",")]
        return li

    def read_Synteny_Network(self):
        if os.path.exists(self.graph_file):
            with open(self.graph_file, "rb") as f:
                G = pickle.load(f)
        else:
            G = nx.Graph()
            JCVIDir = os.path.join(self.gene_index_workDir, "JCVIDir")
            for i in range(len(self.speciesList)):
                for j in range(i, len(self.speciesList)):
                    que, ref = self.speciesList[i], self.speciesList[j]
                    anchor_file = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.anchors")
                    anchor_G = self.read_net(anchor_file)
                    lifted_anchor_file = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.lifted.anchors")
                    lifted_anchor_G = self.read_net(lifted_anchor_file)
                    out_graph = self.overlap_net(anchor_G, lifted_anchor_G)
                    G.add_edges_from(out_graph.edges(data=True))  # The edges attribute must be added
            with open(self.graph_file, "wb") as fo:
                pickle.dump(G, fo)
        return G
