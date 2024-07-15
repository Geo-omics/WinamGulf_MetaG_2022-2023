import re
import requests
import bs4
from Bio import SeqIO
import pandas as pd
import threading
import os
import time


def get_genbank_number(cluster):
    url = f"https://mibig.secondarymetabolites.org/repository/{cluster}"
    html = requests.get(url).text
    soup = bs4.BeautifulSoup(html, 'html.parser')
    genbank_number = None
    all_text = ""
    div = soup.find('div', {'class': 'description-text'})
    if div:
        all_text = div.get_text()

    match = re.search(r'This entry is originally from NCBI GenBank (\S+)', all_text)
    if match:
        genbank_number = match.group(1)[:-1]
    else:
        print("GenBank number not found")

    # find link to blast.ncbi.nlm.nih.gov
    blast_link = url
    
    return genbank_number, blast_link



def get_additional_cluster_info(cluster):
    url = f"https://mibig.secondarymetabolites.org/repository/{cluster}/{cluster}.gbk"
    with open(f"temp/{cluster}.gbk", "w") as f:
        f.write(requests.get(url).text)
    parsed = SeqIO.parse(f"temp/{cluster}.gbk", "genbank")
    gene_info = {"genes": {}}
    gene_number = 1
    for record in parsed:
        actual_start = int(record.annotations['structured_comment']['antiSMASH-Data']['Orig. start']) + 1
        actual_end = int(record.annotations['structured_comment']['antiSMASH-Data']['Orig. end'])
        for feature in record.features:
            if feature.type == "CDS" and 'gene' in feature.qualifiers:
                gene_start = int(feature.location.start) + actual_start
                gene_end = int(feature.location.end) + actual_start
                gene_name = feature.qualifiers['gene'][0]
                gene_functions = []
                if 'gene_functions' in feature.qualifiers:
                    for function in feature.qualifiers['gene_functions']:
                        gene_functions.append(function)
                gene_kind = "other"
                if 'gene_kind' in feature.qualifiers: 
                    gene_kind = feature.qualifiers['gene_kind'][0]
                protein_id = ""
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                gene_info["genes"][gene_name] = {'start': gene_start, 'end': gene_end, 'functions': gene_functions, 'gene-kind': gene_kind, 'number': gene_number, 'protein_id': protein_id}
                gene_number += 1
            elif feature.type == "CDS" and 'codon_start' in feature.qualifiers:
                gene_start = int(feature.location.start) + actual_start
                gene_end = int(feature.location.end) + actual_start
                if 'locus_tag' in feature.qualifiers:
                    gene_name = feature.qualifiers['locus_tag'][0]
                elif 'protein_id' in feature.qualifiers:
                    gene_name = feature.qualifiers['protein_id'][0]
                gene_functions = []
                if 'gene_functions' in feature.qualifiers:
                    for function in feature.qualifiers['gene_functions']:
                        gene_functions.append(function)
                gene_kind = "Unknown"
                if 'gene_kind' in feature.qualifiers:
                    gene_kind = feature.qualifiers['gene_kind'][0]
                protein_id = ""
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                gene_info["genes"][gene_name] = {'start': gene_start, 'end': gene_end, 'functions': gene_functions, 'gene-kind': gene_kind, 'number': gene_number, 'protein_id': protein_id}
                gene_number += 1

                
        source_feature = record.features[0]
        total_length = len(record.seq)
        gene_info['start'] = actual_start
        gene_info['end'] = actual_end
        gene_info['total_length'] = total_length
                
    gene_info['genbank_number'], gene_info['blast_link'] = get_genbank_number(cluster)

    return gene_info


def get_search_results(name):
    url = f"https://mibig.secondarymetabolites.org/api/v1/search"
    params = {"search_string": name}
    response = requests.post(url, json=params)
    return response.json()


def get_cluster_info(cluster):
    inital_results_info = {
        "accession": cluster["accession"],
        "complete": cluster["complete"],
        "mibig_info": cluster["minimal"],
        "products": cluster["products"],
        "classes": cluster["classes"],
        "organism": cluster["organism"],
    }
    
    return inital_results_info

def create_template_excel():
    columns = [
        "Accession",
        "Complete_cluster",
        "MiBIG_Info",
        "Main_Product",
        "Biosynthetic_Class",
        "Organism",
        "Start",
        "End",
        "Total_Length",
        "NCBI_Genbank",
        "Protein_ID",
        "Gene",
        "Gene_start",
        "Gene_end",
        "#_in_BGC",
        "Type",
        "Rule-based",
        "Functions",
        "NCBi_BlastP",
        "Potential_Homologs?",
        "DB_name"
    ]
    df = pd.DataFrame(columns=columns)
    df.to_excel("template.xlsx", index=False)


def process_cluster(cluster, df, completed_dfs):
    # merge cluster info with gene info
    cluster_info = get_cluster_info(cluster)
    gene_info = get_additional_cluster_info(cluster['accession'])
    for gene in gene_info:
        cluster_info[gene] = gene_info[gene]

    # add rows for each gene
    for gene in gene_info['genes']:
        df.loc[len(df)] = ([
            cluster_info['accession'],
            cluster_info['complete'],
            cluster_info['mibig_info'],
            ', '.join(list(map(lambda x: x['name'], cluster_info['products']))),
            ', '.join(list(map(lambda x: x['name'], cluster_info['classes']))),
            cluster_info['organism'],
            cluster_info['start'],
            cluster_info['end'],
            cluster_info['total_length'],
            cluster_info['genbank_number'],
            gene_info['genes'][gene]['protein_id'],
            gene,
            gene_info['genes'][gene]['start'],
            gene_info['genes'][gene]['end'],
            gene_info['genes'][gene]['number'],
            gene_info['genes'][gene]['gene-kind'],
            "",
            ', '.join(gene_info['genes'][gene]['functions']),
            gene_info['blast_link'],
            "",
            "",
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        ])

    completed_dfs.append(df)

def make_df():
    columns = [
    "Accession",
    "Complete_cluster",
    "MiBIG_Info",
    "Main_Product",
    "Biosynthetic_Class",
    "Organism",
    "Start",
    "End",
    "Total_Length",
    "NCBI_Genbank",
    "Protein_ID",
    "Gene",
    "Gene_start",
    "Gene_end",
    "#_in_BGC",
    "Type",
    "Rule-based",
    "Functions",
    "NCBi_BlastP",
    "Potential_Homologs?",
    "DB_name",
    "Date Added",
]
        
    df = pd.DataFrame(columns=columns)

    return df

def process_search_results(name):
    if not os.path.exists("temp"):
        os.mkdir("temp")
    print("Processing search results...")
    completed_dfs = []
    threads = []
    results = get_search_results(name)
    if 'clusters' not in results or len(results['clusters']) == 0:
        print("No results found")
        return
    clusters = results['clusters']
    for cluster in clusters:
        threads.append(threading.Thread(target=process_cluster, args=(cluster, make_df(), completed_dfs)))
        threads[-1].start()

    for thread in threads:
        thread.join()
    
    df = pd.concat(completed_dfs)
    if os.path.exists("mibig_search_results.xlsx"):
        should_extend = None
        while should_extend not in ['y', 'n']:
            should_extend = input("File already exists, do you want to add to it? Selecting 'n' will OVERWRITE the file. (y/n) ")

        if should_extend == 'y':
            old_df = pd.read_excel("mibig_search_results.xlsx")
            df = pd.concat([old_df, df])

    df = df.sort_values(by=['Accession', '#_in_BGC'])
    df.to_excel(f"mibig_search_results.xlsx", index=False)
    
    # clean up all .gbk files using the os module
    for file in os.listdir("temp"):
        if file.endswith(".gbk"):
            os.remove(os.path.join("temp", file))

    return df


if __name__ == "__main__":
    name = "microcystis"
    process_search_results(input("Enter search term: "))