#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import csv
import allel
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse


def make_header(file_name):
    """
    saves the vep file header
        :param file_name: vep file name
    :return header_list[0]: header list
    """

    header_list = []
    with open(file_name, "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            if line.split()[0] == "#Uploaded_variation":  # сохраним заголовок
                header_list.append(line.split())
                return header_list[0]

def omim_list_func(omim_list_file):
    """
    create a list of genes
        :param omim_list_file: file containing raw list of genes
    :return omim_list_2: prepared list of genes
    """


    omim_list = []
    with open(omim_list_file, "r") as file: #создаем список генов
        while True:
            line = file.readline()
            if not line:
                break
            omim_list.append(line)
            omim_list_1 = [i.replace('\n', '') for i in omim_list ]
            omim_list_2 = [i for i in omim_list_1 if 'NA' not in i]
    return omim_list_2


def make_lists(df):
    """
    an iterator that creates the necessary lists
        :param df: dataframe
    :yield ref_amino: reference amino acid
    :yield alt_amino: atlernative amino acid
    :yield ref_codon: reference codon
    :yield alt_codon: atlernative codon
    """


    for i, row in df.iterrows():
        ref_amino = row[6][:1]
        alt_amino = row[6][2:]
        ref_codon = row[7][:3]
        alt_codon = row[7][4:]
        yield ref_amino, alt_amino, ref_codon, alt_codon


def compare_triple(triple_a, triple_b):
    """
    сompare two triples, if they differ by one letter
        :param triple_a: string of first triples
        :param triple_b: string of second triples
    :return tuple of triple_a: string of first triples
        :triple_b: string of second triples
        :list_dif_nucl[0][0]: reference nucleotide
        :list_dif_nucl[0][1]: atlernative nucleotide 
        :list_dif_nucl[0][2]: position
    """


    list_dif_nucl = []
    pos = 0
    for x, y in zip(triple_a, triple_b):
        pos += 1
        if x != y:
            list_dif_nucl.append([x, y, pos])
    if len(list_dif_nucl) == 1:
        return (
            triple_a,
            triple_b,
            list_dif_nucl[0][0],
            list_dif_nucl[0][1],
            list_dif_nucl[0][2],
        )


def compare_ref_with_codons(ref_codon, alt_codons, alt_codon):
    """
    compares the reference codon with alternative codons
        :param ref_codon: reference codon string
        :param alt_codons: list of alternative codons
        :param alt_codon: alternative codon string 
    :return list_compare_triple: a list containing the string of first triples
        string of second triples
        nucleotide reference
        atlernative nucleotide
        position
    """


    list_compare_triple = []
    ref_pos = 1
    alt_pos = 1
    for i in alt_codons:
        list_a = []
        dif = 0
        if compare_triple(ref_codon, i) and i == alt_codon:
            ref_pos = compare_triple(ref_codon, i)[4]
    for i in alt_codons:
        list_a = []
        if compare_triple(ref_codon, i) and i != alt_codon:
            alt_pos = compare_triple(ref_codon, i)[4]
            list_a.append(compare_triple(ref_codon, i)[0:4])
            dif = alt_pos - ref_pos
            list_a.append(dif)
            list_compare_triple.append(list_a)
    return list_compare_triple


def compare_lists(ref_codon_list, alt_amino_list, alt_codon_list):
    """
    compares the list of reference codons with the list of alternative codons
        :param ref_codon_list: reference codons list
        :param alt_amino_list: list of alternative amino acids
        :param alt_codon_list: alternative codons list 
    :return additional_amino_list: a list containing the lists needed for the dataframe
    """


    additional_amino_list = []
    for ref_codon, alt_amino, alt_codon in zip(
        ref_codon_list, alt_amino_list, alt_codon_list
    ):
        try:
            additional_amino_list.append(
                compare_ref_with_codons(
                    ref_codon.upper(), codon_table[alt_amino], alt_codon.upper()
                )
            )  # сравниваем референсный кодон с кортежем альтернативных кодонов
        except KeyError:
            additional_amino_list.append("#")  # Если ошибка ключа помечаем "#"
            pass
    return additional_amino_list


def make_df_list(df):
    """
    creates a list to create a dataframe 
        :param df: df_add
    :return df_list: list containing information for vcf formate
    """


    df_list = []
    for row in df.itertuples():
        for i in range(9, len(row)):
            if row[i] and row[i] != "#":
                str_list = []
                pos = int(row[2]) + int(row[i][1])
                ID = str(row[3]) + "_" + str(i - 8)
                REF = row[i][0][2]
                ALT_1 = row[i][0][3]
                REF_codon = row[i][0][0]
                ALT_codon = row[i][0][1]
                str_list.append(row[1])
                str_list.append(pos)
                str_list.append(ID)
                str_list.append(REF)
                str_list.append(ALT_1)
                str_list.append(".")  # 'QUAL'
                str_list.append(".")  # 'FILTER'
                info_str = ""
                for k in range(6, 8):
                    info_str += str(row[k]) + ";"
                info_str += REF_codon + ";" + ALT_codon
                str_list.append(info_str)
                df_list.append(str_list)
    return df_list


def save_VCF(df, name):
    """
    save dataframe to vcf file 
        :param df: dataframe in the vcf formate
    """


    header = """##fileformat=
##fileDate=
##source=
##reference=
#CHROM POS ID REF ALT QUAL FILTER INFO
"""
    with open(name, "w") as vcf:
        vcf.write(header)
    df.to_csv(name, sep="\t", mode="a", header=False)


codon_table = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "I": ("ATT", "ATC", "ATA"),
    "H": ("CAT", "CAC"),
    "K": ("AAA", "AAG"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    "*": ("TAA", "TAG", "TGA"),
}


def processing(vcf_file, vep_file):
    """
    creates four vcf files: benign.vcf, likely_benign.vcf, pathogenic.vcf, likely_pathogenic.vcf
        :vcf_file: clinvar.vcf
        :vep_file: vep file
    """

    header_list = make_header(vep_file)
    df_clinvar = allel.vcf_to_dataframe(input=vcf_file, fields="*")
    clinvar_ann_df = pd.read_table(
        vep_file, header=None, comment='#'
    )
    clinvar_ann_df.columns = header_list
    omim_list_file = 'omim_list.txt' # https://www.omim.org
    omim_list = omim_list_func(omim_list_file)
    omim_pd = pd.DataFrame(omim_list, columns=['Feature'])
    search_omim_df = clinvar_ann_df.merge(omim_pd, on='Feature', how='right')
    omim_missence_canonical_df = search_omim_df[(search_omim_df['Consequence'] == 'missense_variant') & (search_omim_df['Extra'].str.contains("CANONICAL=YES")==True)]
    omim_missence_canonical_df['#Uploaded_variation']=omim_missence_canonical_df['#Uploaded_variation'].astype(int)
    df_1 = clinvar_ann_df[
        (clinvar_ann_df["Consequence"] == "missense_variant")
        & (clinvar_ann_df["Extra"].str.contains("CANONICAL=YES") == True)
    ]
    df_1["#Uploaded_variation"] = df_1["#Uploaded_variation"].astype(int)
    df_1 = omim_missence_canonical_df.reset_index(drop=True)
    df_2 = df_1.iloc[:, [0, 2, 3, 7, 8, 9, 10, 11]]
    generator = make_lists(df_2)
    all_list = []
    ref_amino_list = []
    alt_amino_list = []
    ref_codon_list = []
    alt_codon_list = []
    i = 0
    while True:
        try:
            all_list.append(next(generator))
            ref_amino_list.append(all_list[i][0])
            alt_amino_list.append(all_list[i][1])
            ref_codon_list.append(all_list[i][2])
            alt_codon_list.append(all_list[i][3])
            i += 1
        except StopIteration:
            break
    additional_amino_list = compare_lists(
        ref_codon_list, alt_amino_list, alt_codon_list
    )
    additional_amino_df = pd.DataFrame(additional_amino_list)
    df_3 = pd.concat([df_2, additional_amino_df], axis=1)
    ncol = len(df_3.columns)
    col_list = [0]
    for i in range(6, ncol):
        col_list.append(i)
    df_4 = df_3.iloc[:, col_list]
    df_5 = df_4.loc[df_4[0].isna() != True]
    df_5.rename(columns={"#Uploaded_variation": "ID"}, inplace=True)
    df_clinvar_2 = df_clinvar.iloc[:, [0, 1, 2, 3, 4, 18]]
    df_clinvar_2["ID"] = df_clinvar_2["ID"].astype(int)
    df_add = df_clinvar_2.merge(df_5, on="ID", how="right")
    df_list = make_df_list(df_add)
    df_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    df_6 = pd.DataFrame(df_list, columns=df_columns)
    df_6.set_index("#CHROM", inplace=True)
    benign_df = df_6[
        (df_6["INFO"].str.contains("Benign") == True)
        & (df_6["INFO"].str.contains("Likely_benign") == False)
    ]
    likely_benign_df = df_6[df_6["INFO"].str.contains("Likely_benign") == True]
    pathogenic_df = df_6[
        (df_6["INFO"].str.contains("Pathogenic") == True)
        & (df_6["INFO"].str.contains("Likely_pathogenic") == False)
    ]
    likely_pathogenic_df = df_6[df_6["INFO"].str.contains("Likely_pathogenic") == True]
    save_VCF(benign_df, "benign.vcf")
    print("benign.vcf    SAVED")
    save_VCF(likely_benign_df, "likely_benign.vcf")
    print("likely_benign.vcf    SAVED")
    save_VCF(pathogenic_df, "pathogenic.vcf")
    print("pathogenic.vcf    SAVED")
    save_VCF(likely_pathogenic_df, "likely_pathogenic.vcf")
    print("likely_pathogenic.vcf    SAVED")
    df_1.rename(columns={"#Uploaded_variation": "ID"}, inplace=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CSV file processing")
    parser.add_argument(
        "-a",
        "--vep_file",
        type=str,
        help="vep file",
        required=True,
    )
    parser.add_argument(
        "-b", "--vcf_file", type=str, help="clinvar.vcf", required=True
    )
    args = parser.parse_args()
    vcf_file = args.vcf_file
    vep_file = args.vep_file
    processing(vcf_file, vep_file)
