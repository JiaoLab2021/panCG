import os
import networkx as nx
import multiprocessing
import pickle

from panCG.lib.CnsGeneLink import *
from panCG import SCIPTS_DIR


def run_cns_gene_link(logger, filter_cns_index, gene_index_file, link_yaml, threads, out_dir):
    BasePan.pymkdir(out_dir)
    link_dict = BasePan.read_yaml(link_yaml)
    species_li = list(link_dict["species"].keys())

    CnsIndexReader = CnsIndexRead(filter_cns_index)
    GeneIndexReader = GeneIndexRead(gene_index_file)

    anno_R_file = os.path.join(SCIPTS_DIR, "CNS.anno.R")
    pool = multiprocessing.Pool(processes=threads)
    for spe in species_li:
        df = CnsIndexReader.ToBed("Index", spe)
        cns_bed_file = os.path.join(out_dir, f"{spe}.bed")
        df.to_csv(cns_bed_file, sep="\t", index=False, header=False)
        anno_file = os.path.join(out_dir, f"{spe}.anno.tsv")
        log_file = os.path.join(out_dir, f"{spe}.anno.log")
        cmd = "Rscript {} -b {} -g {} -o {} > {} 2>&1".format(anno_R_file,
                                                              cns_bed_file,
                                                              link_dict["species"][spe]["longest_pep_gff"],
                                                              anno_file,
                                                              log_file)
        pool.apply_async(BasePan.cmd_linux, (cmd,))
    pool.close()
    pool.join()

    G = nx.DiGraph()  # Network graph of CNS index and Gene index
    anno_data_li = []
    for spe in species_li:
        logger.info(f"parse {spe} ...")
        gene2index_dict = GeneIndexReader.ToDict("Index", spe)
        anno_file = os.path.join(out_dir, f"{spe}.anno.tsv")
        LinkCnsGeneIndexer = LinkCnsGeneIndex(anno_file)
        data = LinkCnsGeneIndexer.parse_anno(gene2index_dict)

        directed_edge_data = data[["seqnames", "start", "end", "transcriptId", "annotation"]].copy()
        directed_edge_data["cns_id"] = directed_edge_data.apply(lambda row: "{}::{}:{}-{}".format(
            spe, row["seqnames"], row["start"], row["end"]), axis=1)
        directed_edge_data = directed_edge_data[["cns_id", "transcriptId", "annotation"]]
        anno_data_li.append(directed_edge_data)

        edge_data = data[["cns_index", "gene_index"]].copy()
        edge_data['gene_index'] = edge_data['gene_index'].str.split(',')
        edge_data = edge_data.explode('gene_index', ignore_index=True)
        edge_data = edge_data.drop_duplicates()

        for _, row in edge_data.iterrows():
            cns_index, gene_index = row['cns_index'], row['gene_index']
            if G.has_edge(cns_index, gene_index):
                G[cns_index][gene_index]['weight'] += 1
            else:
                G.add_edge(cns_index, gene_index, weight=1)

    logger.info(f"Constructing a directed network graph of CNS and genes ...")
    DiG_file = os.path.join(out_dir, f"CnsGeneLink.net.graph.pkl")
    directed_edge_file = os.path.join(out_dir, f"CnsGeneLink.net.graph.tsv")

    directed_edge_df = pd.concat(anno_data_li, axis=0, ignore_index=True)
    directed_edge_df.to_csv(directed_edge_file, sep="\t", index=False)

    DiG = nx.from_pandas_edgelist(
        directed_edge_df,
        source='cns_id',
        target='transcriptId',
        edge_attr=["annotation"],
        create_using=nx.DiGraph()
    )
    with open(DiG_file, "wb") as fo:
        pickle.dump(DiG, fo)

    # Traverse each cns_index node and take the gene_index with the largest number as its target
    row_dict_li = []
    for cns_index_ in CnsIndexReader.cns_index_data["Index"].tolist():
        downstream_edges = G.out_edges(cns_index_, data=True)
        if len(downstream_edges) == 0:
            continue
        max_weight = -float('inf')
        max_node = set()
        # Traverse all downstream edges and find the edge and node with the largest weight
        for u, v, weight in downstream_edges:
            if weight['weight'] > max_weight:
                max_weight = weight['weight']
                max_node = set()
                max_node.add(v)
            elif weight['weight'] == max_weight:
                max_node.add(v)
        row_dict = {"cns_index": cns_index_, "gene_index": ",".join(max_node), "weight": max_weight}
        row_dict_li.append(row_dict)
    result_data = pd.DataFrame(row_dict_li)
    result_data['gene_index'] = result_data['gene_index'].str.split(',')
    result_data = result_data.explode('gene_index', ignore_index=True)

    valid_cns_count_df = CnsIndexReader.get_valid_cns_count(species_li)
    gene_count_df = GeneIndexReader.get_gene_count(species_li)

    result_data = pd.merge(result_data, valid_cns_count_df, on="cns_index", how="left")
    result_data = pd.merge(result_data, gene_count_df, on="gene_index", how="left")

    out_file = os.path.join(out_dir, "CnsGeneLink.tsv")
    result_data.to_csv(out_file, sep="\t", index=False, header=True)

