#!/usr/bin/env python3
import sys
import os
from collections import defaultdict
from collections import OrderedDict
import pandas as pd
import numpy as np

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
DICT_RANK_TO_INDEX = dict(zip(RANKS, list(range(len(RANKS)))))
DICT_RANK_TO_INDEX["kingdom"] = 0
RANKS_LOW2HIGH = list(reversed(RANKS))


def load_merged(names_file_path):
    taxid_to_taxid = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            taxid_to_taxid[line[0]] = line[1]
    return taxid_to_taxid


def load_names(taxid_to_rank, names_file_path):
    taxid_to_name = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            taxid = line[0]

            if line[3] == "scientific name":
                taxid_to_name[taxid] = line[1]
            else:
                continue

            if taxid_to_rank[taxid] == "species" or taxid_to_rank[taxid] == "strain":
                names = taxid_to_name[taxid].split(" ")
                if len(names) > 2:
                    if taxid_to_rank[taxid] == "strain":
                        taxid_to_name[taxid] = "{} {} strain".format(names[0], names[1])
                    else:
                        taxid_to_name[taxid] = "{} {}".format(names[0], names[1])
    taxid_to_name[""] = ""
    return taxid_to_name


def check_parent_is_species(taxid_to_parent, taxid_to_rank, taxid):
    if taxid_to_parent[taxid] in taxid_to_parent:
        if taxid_to_rank[taxid_to_parent[taxid]] == "species":
            return True
        elif taxid_to_parent[taxid] != '1' and taxid_to_rank[taxid_to_parent[taxid]] not in RANKS:
            return check_parent_is_species(taxid_to_parent, taxid_to_rank, taxid_to_parent[taxid])
    return False


def load_tax_info(ncbi_nodes_file):
    taxid_to_parent = {}
    taxid_to_rank = {}
    with open(ncbi_nodes_file) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            taxid = line[0]
            taxid_to_parent[taxid] = line[1]
            taxid_to_rank[taxid] = line[2]

    for taxid, rank in taxid_to_rank.items():
        if taxid_to_rank[taxid] == "no rank" and check_parent_is_species(taxid_to_parent, taxid_to_rank, taxid):
            taxid_to_rank[taxid] = "strain"

    return taxid_to_parent, taxid_to_rank


def get_id_path(taxid, taxid_to_parent, taxid_to_rank):
    if taxid not in taxid_to_rank:
        # TODO report this in a log file
        return None

    while taxid_to_rank[taxid] not in RANKS:
        taxid = taxid_to_parent[taxid]
        if taxid == '1':
            return None

    index = DICT_RANK_TO_INDEX[taxid_to_rank[taxid]]
    path = [''] * (index + 1)
    path[index] = taxid

    id = taxid
    while id in taxid_to_parent:
        id = taxid_to_parent[id]
        if id == '1':
            break
        if taxid_to_rank[id] not in RANKS:
            continue
        index = DICT_RANK_TO_INDEX[taxid_to_rank[id]]
        path[index] = id
        if taxid_to_rank[id] == "superkingdom":
            break
    return path

if __name__ == "__main__":
    taxid_to_parent, taxid_to_rank = load_tax_info(f"{sys.argv[3]}/nodes.dmp")
    taxid_to_name = load_names(taxid_to_rank, f"{sys.argv[3]}/names.dmp")

    df_all = pd.read_csv(sys.argv[2], sep="\t")
    df_all["TAXONOMY_ID"] =  df_all["TAXONOMY_ID"].astype(str)

    th = float(sys.argv[1])
    df_species = df_all.loc[df_all["TAXONOMY_LEVEL"] == "species"]
    df_species = df_species.loc[df_species["FRACTION_TOTAL"].astype(float) > th]

    species_set = set(df_species["TAXONOMY_ID"]) - {'0'}
    parent_taxa = set()
    
    for taxid in species_set:
        parent_taxa.update(get_id_path(taxid, taxid_to_parent, taxid_to_rank))
    df_all = df_all.loc[df_all["TAXONOMY_ID"].isin(parent_taxa) | (df_all["FRACTION_TOTAL"].astype(float) > th)]
    # df_all = df_all.loc[df_all["TAXONOMY_ID"].isin(parent_taxa)]

    df_genome_info = pd.read_csv(sys.argv[4], sep="\t")
    df_genome_info = df_genome_info[["species_taxid", "genome_size"]]
    gsize = defaultdict(list)
    for index, row in df_genome_info.iterrows():
        path = get_id_path(str(row["species_taxid"]), taxid_to_parent, taxid_to_rank)
        if path is not None:
            for taxid in path:
                gsize[taxid].append(row["genome_size"])
    gsize = {key: sum(val)/len(val) for key, val in gsize.items()}
    gsize["0"] = df_genome_info["genome_size"].mean()

    df_all["AVG_GENOME_SIZE"] = df_all["TAXONOMY_ID"].apply(lambda x: gsize.get(x,  gsize["0"]))
    df_all["FRACTION_TOTAL"] = df_all["FRACTION_TOTAL"] / df_all["AVG_GENOME_SIZE"]
    df_all["FRACTION_TOTAL"] = df_all["FRACTION_TOTAL"] / df_all["FRACTION_TOTAL"].sum()

    rank_to_taxid_to_score = defaultdict(dict)
    rank_to_taxid_to_id_path = defaultdict(dict)

    for ix, row in df_all.iterrows():
        taxid = row["TAXONOMY_ID"]
        curr_rank = row["TAXONOMY_LEVEL"]
        rank_index = DICT_RANK_TO_INDEX[curr_rank]
        rank_to_taxid_to_score[rank_index][taxid] = row["FRACTION_TOTAL"]
        rank_to_taxid_to_id_path[rank_index][taxid] = get_id_path(taxid, taxid_to_parent, taxid_to_rank)


    sample_id = sys.argv[5]
    with open(sys.argv[2] + ".profile", "w") as f:
        f.write("@SampleID:{}\n@Version:0.9.1\n@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n".format(sample_id))
        for rank in RANKS:
            rank_index = DICT_RANK_TO_INDEX[rank]
            for taxid in rank_to_taxid_to_score[rank_index]:
                if rank_to_taxid_to_score[rank_index][taxid] == .0:
                    continue
                if taxid != "0":
                    id_path = rank_to_taxid_to_id_path[rank_index][taxid]
                    name_path = []
                    for id in id_path:
                        name_path.append(taxid_to_name[id])
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(taxid, taxid_to_rank[taxid], "|".join(id_path), "|".join(name_path), rank_to_taxid_to_score[rank_index][taxid] * 100))
