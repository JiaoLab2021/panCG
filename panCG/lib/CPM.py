
import networkx as nx
import numpy as np
import logging
import pandas as pd


class CPM:
    # Since the cmp algorithm consumes a lot of resources for communities with a huge number of cliques,
    # a maximum threshold is set here. Communities exceeding the threshold are not processed and are subsequently
    # classified according to their evolutionary relationships.
    k = 3
    max_cliques_num = 100000  # Due to excessive resource consumption, CPM is not applicable if the number of cliques exceeds max_cliques_num
    large_cliques_num = 5000
    chunksize = 100  # If there are too many cliques (more than large_cliques_num), chunksize cliques are used as the initial

    def __init__(self, G):
        self.G = G

    @classmethod
    def set_max_cliques_num(cls, value_):
        cls.max_cliques_num = value_

    @classmethod
    def set_large_cliques_num(cls, value_):
        cls.large_cliques_num = value_

    @classmethod
    def set_k(cls, value_):
        cls.k = value_

    @staticmethod
    def get_dtype(max_value):
        if max_value <= pow(2, 8) - 1:  # 255
            dtype = np.uint8
        elif max_value <= pow(2, 16) - 1:  # 65535
            dtype = np.uint16
        elif max_value <= pow(2, 32) - 1:  # 4294967295
            dtype = np.uint32
        else:
            dtype = np.uint64
        return dtype

    @staticmethod
    def set_to_vector(all_node_li, set_):
        return np.array([1 if elem in set_ else 0 for elem in all_node_li])

    def nodes2matrix(self, all_nodes: set, set_li: list):
        node_li = sorted(list(all_nodes))
        # Replace set with a vector, 1 indicates existence, 0 indicates non-existence,
        # and the position indicates the position of the element in the sorted all_nodes
        vector_list = [self.set_to_vector(node_li, set_) for set_ in set_li]
        matrix = np.array(vector_list, dtype=self.get_dtype(len(all_nodes)))
        result_matrix = np.dot(matrix, matrix.T)
        return result_matrix

    def find_cliques_percolation(self):
        """
        Implements Clique Percolation Method (CPM) for community detection.

        Args:
        G : NetworkX graph
            The input graph.

        Returns:
        list of sets
            A list of communities, where each community is a set of nodes.
        """

        # # Step 1: Find all k-cliques
        # cliques = nx.find_cliques(self.G)  # Calculate the maximum cliques
        # cliques = [set(clique) for clique in cliques if len(clique) >= self.k]
        # # if len(cliques) > self.large_cliques_num:
        # #     # logging.warning(f"The community in k-cliques exceeds the threshold {self.large_cliques_num}. len(cliques) = {len(cliques)}. len(nodes) = {len(set(self.G.nodes()))}")
        # #     return [set(self.G.nodes())]
        #
        # # Step 2: Build the Clique Graph
        # clique_graph = nx.Graph()
        # clique_graph.add_nodes_from(range(len(cliques)))
        # for (i, set1), (j, set2) in combinations(enumerate(cliques), 2):
        #     if len(set1.intersection(set2)) >= self.k - 1:
        #         clique_graph.add_edge(i, j)
        #
        # # Step 3: Find connected components in the Clique Graph
        # communities = []
        # for component in nx.connected_components(clique_graph):
        #     community = set()
        #     for clique_index in component:
        #         community.update(cliques[clique_index])
        #     communities.append(community)
        # return communities

        # Step 1: Find all k-cliques
        cliques = nx.find_cliques(self.G)  # Calculate the maximum cliques
        cliques = [set(clique) for clique in cliques if len(clique) >= self.k]

        if len(cliques) >= self.max_cliques_num:
            set_ = set(self.G.nodes())
            logging.warning("The community in k-cliques exceeds the threshold {}. "
                            "len(cliques) = {}. "
                            "len(nodes) = {}. "
                            "The list is {}".format(self.max_cliques_num,
                                                    len(cliques),
                                                    len(set(self.G.nodes())),
                                                    ",".join(set_)
                                                    )
                            )
            return [set_]

        if len(cliques) <= 1:
            communities = []
            for i in cliques:
                communities.append(i)
            return communities

        # Step 2: Build the Clique Graph
        all_nodes = set(self.G.nodes())
        result_matrix = self.nodes2matrix(all_nodes, cliques)
        if len(cliques) >= self.large_cliques_num:
            chunk_m = result_matrix[0:self.chunksize, :]
            row_indices, col_indices = np.where(chunk_m >= self.k - 1)
            values = result_matrix[row_indices, col_indices]
            df = pd.DataFrame({
                'row_index': row_indices,
                'col_index': col_indices,
                'value': values
            })
            clique_graph = nx.from_pandas_edgelist(df, 'row_index', 'col_index')
            # Add the edge of its own node
            clique_graph.add_edges_from([(i, i) for i in range(result_matrix.shape[0])])
            # Convert the connected graph into a loop graph to simplify calculations
            huan_G = nx.Graph()
            for component in list(nx.connected_components(clique_graph)):
                component_li = list(component)
                edges = [(component_li[i], component_li[(i + 1) % len(component_li)]) for i in range(len(component_li))]
                huan_G.add_edges_from(edges)
            clique_graph, huan_G = huan_G, clique_graph
            del huan_G
            len_ = 0
            while True:
                connected_components = list(nx.connected_components(clique_graph))
                if len(connected_components) == len_:  # Indicates that no more edges can be added
                    break
                if len(connected_components) == 1:  # It means there is only one
                    break
                len_ = len(connected_components)
                for i in range(len(connected_components)):
                    for j in range(i + 1, len(connected_components)):
                        chunk_m = result_matrix[list(connected_components[i]), :][:, list(connected_components[j])]
                        if np.max(chunk_m) >= self.k - 1:
                            # Randomly pick one of the two connected graphs and add an edge
                            clique_graph.add_edge(list(connected_components[i])[0], list(connected_components[j])[0])
        else:
            row_indices, col_indices = np.where(result_matrix >= self.k - 1)
            values = result_matrix[row_indices, col_indices]
            df = pd.DataFrame({
                'row_index': row_indices,
                'col_index': col_indices,
                'value': values
            })
            clique_graph = nx.from_pandas_edgelist(df, 'row_index', 'col_index')

        # Step 3: Find connected components in the Clique Graph
        communities = []
        for component in nx.connected_components(clique_graph):
            community = set()
            for clique_index in component:
                community.update(cliques[clique_index])
            communities.append(community)
        return communities
