import pandas as pd
import copy
import networkx as nx
import heapq

from panCG.lib.base import DetermineReference
from panCG.lib.CPM import CPM


def multi_IndexAssign(df, speciesList, JCVI_Dict):
    """
    This function uses the results of jcvi to assign indexes to the clustered cns
    Args:
        df: cns cluster
        speciesList: species list
        JCVI_Dict: jcvi dict
    Returns:
        assigned cns
    """
    data_li = []
    for _, row in df.iterrows():
        new_Group_Dict = {}
        row_dict = row.to_dict()
        columnsLi = ["Group", "Index"]
        columnsLi.extend(speciesList)
        check_n = 0
        while True:
            check_n += 1
            if check_n > len(speciesList) + 1:
                raise Exception("The number of loops exceeds the threshold")
            new_Group_Dict_copy = copy.deepcopy(new_Group_Dict)
            Maps_Dict = {key: [] for key in speciesList}
            for k, v in new_Group_Dict.items():
                for specie in speciesList:
                    Maps_Dict[specie].extend(v.get(specie, []))
            UnMapDict = {key: [] for key in speciesList}
            for specie in speciesList:
                mapCNSsList = Maps_Dict[specie]
                allCNSsList = row_dict[specie].split(",") if row_dict[specie] != "." else []
                UnMapDict[specie] = list(set(allCNSsList) - set(mapCNSsList))
            if all([i == [] for i in UnMapDict.values()]):
                del new_Group_Dict_copy
                break
            else:
                reference = DetermineReference(UnMapDict, speciesList)
                Start = len(new_Group_Dict) + 1
                for count, refCNS in enumerate(UnMapDict[reference], start=Start):
                    index_cns = "{}.{}".format(row_dict["Group"], count)
                    new_Group_Dict_copy[index_cns] = {key: [] for key in columnsLi}
                    new_Group_Dict_copy[index_cns]["Group"].append(row_dict["Group"])
                    new_Group_Dict_copy[index_cns]["Index"].append(index_cns)
                    new_Group_Dict_copy[index_cns][reference].append(refCNS)
                    for specie in speciesList:
                        if specie == reference:
                            continue
                        allCNSsList = row_dict[specie].split(",") if row_dict[specie] != "." else []
                        hitJCVICNSsList = JCVI_Dict[f"{reference}.{specie}"].get(refCNS, [])
                        new_Group_Dict_copy[index_cns][specie].extend(
                            list(set(hitJCVICNSsList).intersection(set(allCNSsList))))
                new_Group_Dict, new_Group_Dict_copy = new_Group_Dict_copy, new_Group_Dict
                del new_Group_Dict_copy

        for v_dict in new_Group_Dict.values():
            newDict = {k1: ",".join(list(set(v1))) if v1 else "." for k1, v1 in v_dict.items()}
            data_li.append(pd.DataFrame([newDict]))
    return pd.concat(data_li, ignore_index=True, axis=0)


class GeneIndex:
    single_abstract_cutoff = 5
    top_max_socre_num = 3

    def __init__(self, diamond_blastp_db, Orthogroups_data, G, k):
        """
        Args:
            diamond_blastp_db: diamond_blastp_db
            Orthogroups_data: Orthogroups_data
            G: Phylogenomic Synteny Network
            k: K value in cpm algorithm
        """
        self.diamond_blastp_db = diamond_blastp_db
        self.Orthogroups_data = Orthogroups_data
        self.G = G
        self.k = k

    @classmethod
    def set_single_abstract_cutoff(cls, value_):
        cls.single_abstract_cutoff = value_

    @staticmethod
    def classify_scattered_nodes(component_sub_G, communities_dict, non_cluster_nodes):
        """
        Determine which communities each given node should belong to based on the synteny_weight and blastp_weight
        between the given node and communities_dict
        Args:
            component_sub_G:
            communities_dict:
            non_cluster_nodes:
        Returns:
        """
        non_cluster_node_dict = {}  # Stores which cluster non_cluster_node belongs to
        for node in non_cluster_nodes:
            dict_ = {}  # Store the best edge weight of the non_cluster_node for each cluster
            for idx, cluster in communities_dict.items():
                tmp_dict = {}  # Store the best edge weight of this node for each node in the cluster
                for i in cluster:
                    if component_sub_G.has_edge(node, i):
                        attr_dict = component_sub_G.get_edge_data(node, i)
                        tmp_dict[i] = {"synteny_weight": attr_dict.get("synteny_weight", 0),
                                       "blastp_weight": attr_dict.get("blastp_weight", 0)}
                if len(tmp_dict) > 0:
                    max_key = max(tmp_dict, key=lambda k: (tmp_dict[k]["synteny_weight"], tmp_dict[k]["blastp_weight"]))
                    dict_[idx] = tmp_dict[max_key]
            if len(dict_) > 0:
                cluster_max_key = max(dict_, key=lambda k: (dict_[k]["synteny_weight"], dict_[k]["blastp_weight"]))
                all_max_keys = [
                    k for k, v in dict_.items()
                    if v["synteny_weight"] == dict_[cluster_max_key]["synteny_weight"] and v["blastp_weight"] == dict_[cluster_max_key]["blastp_weight"]
                ]
                for i in all_max_keys:
                    non_cluster_node_dict.setdefault(node, []).append(i)
        return non_cluster_node_dict

    @staticmethod
    def merge_df(df):
        """
        This function is used to merge rows that only have genes in the same column
        Args:
            df:
        Returns:
        """
        data_li = [df[df.apply(lambda row: sum(row != '.') > 1, axis=1)]]
        tmp_df = df[df.apply(lambda row: sum(row != '.') == 1, axis=1)]
        for col in tmp_df.columns:
            df_ = tmp_df[tmp_df[col] != '.']
            if len(df_) == 0:
                continue
            series_tmp = df_.apply(lambda x: ','.join(x[x != '.']) if (x != '.').any() else '.', axis=0)
            df_tmp = pd.DataFrame([series_tmp], columns=df_.columns)
            data_li.append(df_tmp)
        result = pd.concat(data_li, axis=0, ignore_index=True)
        return result

    def parse_row(self, row, group_index_):
        node_df = (row.apply(lambda x: pd.Series(x.split(',')))
                   .stack()
                   .reset_index(level=1, drop=True)
                   .to_frame(name='gene_id')
                   .reset_index()
                   .rename(columns={'index': 'species'}))
        node_df = node_df[node_df["gene_id"] != "."]
        subgraph = self.G.subgraph(node_df["gene_id"].tolist()).copy()

        # step 1. Extract the Phylogenomic Synteny Network in the group and add a separate node
        gene_dict = dict(zip(node_df['gene_id'], node_df['species']))
        sub_G = nx.Graph(subgraph)  # Deep copy of the graph
        scattered_nodes = set(node_df["gene_id"].tolist()) - set(subgraph.nodes())
        for i in scattered_nodes:
            sub_G.add_node(i)

        # step 2. Extract the connected subnetworks in the Phylogenomic Synteny Network, and abstract each subnetwork into a node.
        abstract_nodes = list(nx.connected_components(sub_G))
        # Dictionary storing index abstract_nodes, key is used to build the network
        abstract_nodes_dict = {index_: abstract_nodes for index_, abstract_nodes in enumerate(abstract_nodes)}
        single_abstract_nodes_keys = []
        normal_abstract_nodes_keys = []
        for index_, abstract_nodes in abstract_nodes_dict.items():
            spe_li = [gene_dict[i] for i in abstract_nodes]
            if len(set(spe_li)) <= self.single_abstract_cutoff:
                single_abstract_nodes_keys.append(index_)
            else:
                normal_abstract_nodes_keys.append(index_)

        # step 3. For each single net, add an edge to other networks based on the subnet’s best hit
        for single_abstract_nodes_key in single_abstract_nodes_keys:
            single_abstract_nodes_set = abstract_nodes_dict[single_abstract_nodes_key]
            result = []  # [(que1, ref1, score), (que1, ref2, score))]
            for i in single_abstract_nodes_set:
                top_li = heapq.nlargest(self.top_max_socre_num, self.diamond_blastp_db.query_data(i),
                                        key=lambda x: x[1])
                result.extend([(i, *j) for j in top_li])
            if len(result) == 0:
                continue
            all_top_li = sorted(result, key=lambda x: x[2], reverse=True)
            for tuple_i in all_top_li:
                if tuple_i[1] in gene_dict:
                    sub_G.add_edge(tuple_i[0], tuple_i[1], blastp_weight=tuple_i[2])
                    break
                else:
                    pass

        # step 4. Get the connected subnet of Synteny and blastp net and traverse it
        df_li = []
        for component in nx.connected_components(sub_G):
            component_sub_G = sub_G.subgraph(component).copy()
            CPM.set_k(self.k)
            CPMer = CPM(component_sub_G)
            communities = CPMer.find_cliques_percolation()
            communities_dict = {index_: communities for index_, communities in enumerate(communities)}
            if len(communities_dict) > 1:
                input_dict = copy.deepcopy(communities_dict)
                while True:
                    all_cluster_nodes = set().union(*communities_dict.values())
                    non_cluster_nodes = set(component_sub_G.nodes()) - all_cluster_nodes
                    if len(non_cluster_nodes) == 0:
                        break
                    non_cluster_node_dict = self.classify_scattered_nodes(component_sub_G, input_dict, non_cluster_nodes)
                    input_dict = {}
                    for node, idx_li in non_cluster_node_dict.items():
                        for idx_ in idx_li:
                            communities_dict[idx_].add(node)
                            input_dict.setdefault(idx_, set()).add(node)
            else:  # There is a connected graph with 1 or 0 clusters
                all_cluster_nodes = set().union(*communities_dict.values())
                non_cluster_nodes = set(component_sub_G.nodes()) - all_cluster_nodes
                for i in non_cluster_nodes:
                    communities_dict.setdefault(0, set()).add(i)
            tmp_li = []
            for gene_set in communities_dict.values():
                tmp_dict = {}
                for i in gene_set:
                    tmp_dict.setdefault(gene_dict[i], set()).add(i)
                tmp_li.append({spe: ",".join(genes) for spe, genes in tmp_dict.items()})
            df = pd.DataFrame(tmp_li)
            df_li.append(df)

        result = pd.concat(df_li, axis=0, ignore_index=True)
        result = result.fillna('.')
        result = self.merge_df(result)
        df = pd.DataFrame({
            'Group': [group_index_] * len(result),
            'Index': [f'{group_index_}.{i + 1}' for i in range(len(result))]
        })
        result = pd.concat([df, result], axis=1)
        result["Index"] = result.apply(lambda row: str(row['Index']) + ".Un" if sum([1 for cell in row[2:] if cell != '.']) == 1 else row['Index'], axis=1)
        return result

    def run(self):
        df_li = []
        for index, row in self.Orthogroups_data.iterrows():
            count_non_dot = (row != '.').sum()
            if count_non_dot == 1:
                df = pd.DataFrame([row])
                df.insert(0, "Group", index)
                df.insert(1, "Index", f"{index}.1")
                df_li.append(df)
            else:
                df = self.parse_row(row, index)
                df_li.append(df)
        result = pd.concat(df_li, axis=0, ignore_index=True)
        return result



