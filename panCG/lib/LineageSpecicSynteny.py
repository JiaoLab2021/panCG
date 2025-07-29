import pandas as pd
import networkx as nx
import os
import pickle
from panCG.lib.CPM import CPM
import logging
import uuid


class GeneLineageSpecicSynteny:
    k = 3  # k value in CPM algorithm
    miss_spe_num = 0  # Maximum number of species allowed to lack synteny in Lineage_species
    max_gene_num = 500

    def __init__(self, G, Lineage_species, Orthogroups_data):
        self.G = G
        self.Lineage_species = Lineage_species
        self.Orthogroups_data = Orthogroups_data
        self.species = list(Orthogroups_data.columns)

    @classmethod
    def set_k(cls, value_):
        cls.k = value_

    @classmethod
    def set_miss_spe_num(cls, value_):
        cls.miss_spe_num = value_

    def get_gene2spe_dict(self, row):
        gene2spe_dict = {}
        for spe in self.species:
            if row[spe] != ".":
                for i in row[spe].split(","):
                    gene2spe_dict[i] = spe
        return gene2spe_dict

    def parse_row(self, gene2spe_dict, gene_set, group_index_):
        subgraph = self.G.subgraph(gene_set).copy()

        CPM.set_k(self.k)
        CPMer = CPM(subgraph)
        communities = CPMer.find_cliques_percolation()  # [{node1, node2}, {node3, node4, node5}]

        data_li = []
        all_data_li = []
        for index_, community in enumerate(communities):
            if len(community) > self.max_gene_num:
                logging.warning("No processing for too many Gene")
                continue

            spes_set = set([gene2spe_dict[i] for i in community])
            sub_G = subgraph.subgraph(community).copy()
            edges_df = pd.DataFrame(sub_G.edges, columns=["Source", "Target"])
            edges_df["Source_specie"] = edges_df["Source"].apply(lambda x: gene2spe_dict[x])
            edges_df["Target_specie"] = edges_df["Target"].apply(lambda x: gene2spe_dict[x])
            edges_df["group"] = f"{group_index_}.{index_}"
            all_data_li.append(edges_df)
            if spes_set.issubset(set(self.Lineage_species)) and len(spes_set) >= len(set(self.Lineage_species)) - self.miss_spe_num:
                data_li.append(edges_df)
            else:
                pass

        if len(data_li) != 0:
            data = pd.concat(data_li, axis=0, ignore_index=True)
        else:
            data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        if len(all_data_li) != 0:
            all_data = pd.concat(all_data_li, axis=0, ignore_index=True)
        else:
            all_data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        return data, all_data

    def run(self):
        df_li = []
        all_df_li = []
        for index, row in self.Orthogroups_data.iterrows():
            node_df = (row.apply(lambda x: pd.Series(x.split(',')))
                       .stack()
                       .reset_index(level=1, drop=True)
                       .to_frame(name='gene_id')
                       .reset_index()
                       .rename(columns={'index': 'species'}))
            node_df = node_df[node_df["gene_id"] != "."]
            gene_set = set(node_df["gene_id"].tolist())

            gene2spe_dict = self.get_gene2spe_dict(row)
            df, all_df = self.parse_row(gene2spe_dict, gene_set, index)
            all_df_li.append(all_df)

            not_dot_names = set(row[row != "."].index.tolist()) & set(self.Lineage_species)
            if len(not_dot_names) >= len(set(self.Lineage_species)) - self.miss_spe_num:
                df_li.append(df)
            else:
                pass

        if len(df_li) != 0:
            df = pd.concat(df_li, axis=0, ignore_index=True)
        else:
            df = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        if len(all_df_li) != 0:
            all_df = pd.concat(all_df_li, axis=0, ignore_index=True)
        else:
            all_df = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        return df, all_df


class GeneLineageSpeSynNoOrth:
    k = 3  # k value in CPM algorithm
    miss_spe_num = 0  # Maximum number of species allowed to lack synteny in Lineage_species
    max_gene_num = 500

    def __init__(self, G, Lineage_species):
        self.G = G
        self.Lineage_species = Lineage_species

    @classmethod
    def set_k(cls, value_):
        cls.k = value_

    @classmethod
    def set_miss_spe_num(cls, value_):
        cls.miss_spe_num = value_

    def run(self, gene2spe_dict):
        cluster_df_li = []
        data_li = []
        connected_components = list(nx.connected_components(self.G))
        for n, gene_set in enumerate(connected_components):
            if len(gene_set) > self.max_gene_num:
                logging.warning("No processing for too many genes")
                continue
            sub_G = self.G.subgraph(gene_set).copy()
            CPM.set_k(self.k)
            CPMer = CPM(sub_G)
            communities = CPMer.find_cliques_percolation()
            for community in communities:
                # Store the edges corresponding to the community
                community_G = sub_G.subgraph(community).copy()
                community_G_edges_df = pd.DataFrame(community_G.edges, columns=["Source", "Target"])
                community_G_edges_df["Source_specie"] = community_G_edges_df["Source"].apply(lambda x: gene2spe_dict[x])
                community_G_edges_df["Target_specie"] = community_G_edges_df["Target"].apply(lambda x: gene2spe_dict[x])
                group_ = uuid.uuid4().hex  # Generate global random numbers
                community_G_edges_df["group"] = group_
                cluster_df_li.append(community_G_edges_df)

                # The output includes the cluster of the specified Lineage Species Synteny
                spes_set = set([gene2spe_dict[i] for i in community])
                if (spes_set.issubset(set(self.Lineage_species)) and
                        len(spes_set) >= len(set(self.Lineage_species)) - self.miss_spe_num):
                    data_li.append(community_G_edges_df)
                else:
                    pass

        if len(cluster_df_li) != 0:
            all_cluster_data = pd.concat(cluster_df_li, axis=0, ignore_index=True)
        else:
            all_cluster_data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        if len(data_li) != 0:
            Lineage_Species_data = pd.concat(data_li, axis=0, ignore_index=True)
        else:
            Lineage_Species_data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        return all_cluster_data, Lineage_Species_data


class CnsPhySynNet:
    def __init__(self, cns_index_workDir, graph_file, speciesList):
        self.cns_index_workDir = cns_index_workDir
        self.graph_file = graph_file
        self.speciesList = speciesList

    def build_net(self):
        if os.path.exists(self.graph_file):
            with open(self.graph_file, "rb") as f:
                G = pickle.load(f)
        else:
            G = nx.Graph()
            for i in range(len(self.speciesList)):
                for j in range(i + 1, len(self.speciesList)):
                    que, ref = self.speciesList[i], self.speciesList[j]
                    syn_file = os.path.join(self.cns_index_workDir, "JCVIDir", f"{que}.{ref}", f"{que}.{ref}.halLiftoverFilter.rescue.lifted.anchors")
                    anchorData = pd.read_csv(syn_file, sep="\t", header=None, names=[que, ref])
                    anchorData[que] = anchorData[que].apply(lambda x: f"{que}::{x}")
                    anchorData[ref] = anchorData[ref].apply(lambda x: f"{ref}::{x}")
                    graph = nx.from_pandas_edgelist(anchorData, source=que, target=ref)
                    G.update(graph)
            with open(self.graph_file, "wb") as fo:
                pickle.dump(G, fo)
        return G


class CnsLineageSpecicSynteny:
    k = 3
    miss_spe_num = 0  # Maximum number of species allowed to lack collinearity in Lineage_species
    max_CNS_num = 300

    def __init__(self, G, Lineage_species, connected_component_dict):
        self.G = G
        self.Lineage_species = Lineage_species
        self.connected_component_dict = connected_component_dict

    @classmethod
    def set_k(cls, value_):
        cls.k = value_

    @classmethod
    def set_miss_spe_num(cls, value_):
        cls.miss_spe_num = value_

    def parse_row(self, components, group_index_):
        subgraph = self.G.subgraph(components).copy()

        CPM.set_k(self.k)
        CPM.set_max_cliques_num(10000)
        CPMer = CPM(subgraph)
        communities = CPMer.find_cliques_percolation()  # [{node1, node2}, {node3, node4, node5}]
        data_li = []
        all_data_li = []
        for index_, community in enumerate(communities):
            spes_set = set([i.split("::")[0] for i in community])
            sub_G = subgraph.subgraph(community).copy()
            edges_df = pd.DataFrame(sub_G.edges, columns=["Source", "Target"])
            edges_df["Source_specie"] = edges_df["Source"].apply(lambda x: x.split("::")[0])
            edges_df["Target_specie"] = edges_df["Target"].apply(lambda x: x.split("::")[0])
            edges_df["group"] = f"{group_index_}.{index_}"
            all_data_li.append(edges_df)
            if spes_set.issubset(set(self.Lineage_species)) and len(spes_set) >= len(set(self.Lineage_species)) - self.miss_spe_num:
                data_li.append(edges_df)
            else:
                pass
        if len(data_li) != 0:
            data = pd.concat(data_li, axis=0, ignore_index=True)
        else:
            data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        if len(all_data_li) != 0:
            all_data = pd.concat(all_data_li, axis=0, ignore_index=True)
        else:
            all_data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        return data, all_data

    def run(self):
        df_li = []
        all_df_li = []

        for components_index, components in self.connected_component_dict.items():
            if len(components) > self.max_CNS_num:
                logging.warning("No processing for too many CNS")
                continue
            # else:
            #     print("num of cns_node: ", len(components))

            df, all_df = self.parse_row(components, components_index)
            df_li.append(df)
            all_df_li.append(all_df)

        if len(df_li) != 0:
            data = pd.concat(df_li, axis=0, ignore_index=True)
        else:
            data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        if len(all_df_li) != 0:
            all_data = pd.concat(all_df_li, axis=0, ignore_index=True)
        else:
            all_data = pd.DataFrame(columns=["Source", "Target", "Source_specie", "Target_specie", "group"])

        return data, all_data

