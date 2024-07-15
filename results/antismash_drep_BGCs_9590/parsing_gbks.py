from Bio import SeqIO
import csv

# Function to extract information from a feature
def extract_feature_info(feature):
    start = feature.location.start
    end = feature.location.end
    length = len(feature)
    product = feature.qualifiers.get('product', [''])[0]
    locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
    note = feature.qualifiers.get('note', [''])[0]
    gene_functions_list = feature.qualifiers.get('gene_functions', [])  # Store gene functions in a list
    description = feature.qualifiers.get('description', [''])[0]
    gene_ontologies = feature.qualifiers.get('gene_ontologies', [''])[0]
    gene_kind = feature.qualifiers.get('gene_kind', [''])[0]  # Collect gene_kind
    
    # Join multiple gene functions into a single string
    gene_functions = ';'.join(gene_functions_list)

    return [start, end, length, product, locus_tag, note, gene_functions, description, gene_ontologies, gene_kind]

# Open the output CSV file
with open('drep95_antismash.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Write header row
    csvwriter.writerow(['BGC Locus', 'BGC Definition', 'Start', 'End', 'Length', 'Product', 'Locus Tag', 'Note', 'Gene Functions', 'Description', 'Gene Ontologies', 'Gene Kind'])
    
    # Parse the .gbk file
    for record in SeqIO.parse("/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/data/quality_bins/drep_99/antismash_bins99/antismash_drep_BGCs/extractedbgcs_95_90_rep_seq/extractedbgcs_95_90_rep_seq.gbk", "genbank"):
        bgc_locus = record.name
        bgc_definition = record.description
        
        # Keep track of previously encountered locus tags and their descriptions
        locus_tag_description = {}
        
        # Iterate over features in the record
        for feature in record.features:
            if feature.type == 'CDS':  # Only process CDS features
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                description = feature.qualifiers.get('description', [''])[0]
                
                # Check if locus_tag has been encountered before
                if locus_tag in locus_tag_description:
                    # If a description exists for this locus_tag, add it to the second locus column
                    if locus_tag_description[locus_tag]:
                        feature_info = extract_feature_info(feature)
                        csvwriter.writerow([bgc_locus, bgc_definition] + feature_info + [''] + [locus_tag_description[locus_tag]])
                        locus_tag_description[locus_tag] = ''  # Clear the description to avoid repeating it
                    else:
                        feature_info = extract_feature_info(feature)
                        csvwriter.writerow([bgc_locus, bgc_definition] + feature_info)
                else:
                    # Store the description for this locus_tag
                    locus_tag_description[locus_tag] = description
                    feature_info = extract_feature_info(feature)
                    csvwriter.writerow([bgc_locus, bgc_definition] + feature_info)



# from Bio import SeqIO
# import csv

# # Function to extract information from a feature
# def extract_feature_info(feature):
#     start = feature.location.start
#     end = feature.location.end
#     length = len(feature)
#     product = feature.qualifiers.get('product', [''])[0]
#     locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
#     note = feature.qualifiers.get('note', [''])[0]
#     gene_functions = feature.qualifiers.get('gene_functions', [''])[0]
#     description = feature.qualifiers.get('description', [''])[0]
#     gene_ontologies = feature.qualifiers.get('gene_ontologies', [''])[0]

#     return [start, end, length, product, locus_tag, note, gene_functions, description, gene_ontologies]

# # Open the output CSV file
# with open('output.csv', 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
    
#     # Write header row
#     csvwriter.writerow(['BGC Locus', 'BGC Definition', 'Start', 'End', 'Length', 'Product', 'Locus Tag', 'Note', 'Gene Functions', 'Description', 'Gene Ontologies', 'Gene Description Second Locus'])
    
#     # Parse the .gbk file
#     for record in SeqIO.parse("/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/data/quality_bins/drep_99/antismash_bins99/antismash_drep_BGCs/extractedbgcs_95_90_rep_seq/samp_4444_metabat2_75_contig_746_microviridin.region001.gbk", "genbank"):
#         bgc_locus = record.name
#         bgc_definition = record.description
        
#         # Keep track of previously encountered locus tags and their descriptions
#         locus_tag_description = {}
        
#         # Iterate over features in the record
#         for feature in record.features:
#             if feature.type == 'CDS':  # Only process CDS features
#                 locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
#                 description = feature.qualifiers.get('description', [''])[0]
                
#                 # Check if locus_tag has been encountered before
#                 if locus_tag in locus_tag_description:
#                     # If a description exists for this locus_tag, add it to the second locus column
#                     if locus_tag_description[locus_tag]:
#                         feature_info = extract_feature_info(feature)
#                         csvwriter.writerow([bgc_locus, bgc_definition] + feature_info + [''] + [locus_tag_description[locus_tag]])
#                         locus_tag_description[locus_tag] = ''  # Clear the description to avoid repeating it
#                     else:
#                         feature_info = extract_feature_info(feature)
#                         csvwriter.writerow([bgc_locus, bgc_definition] + feature_info)
#                 else:
#                     # Store the description for this locus_tag
#                     locus_tag_description[locus_tag] = description
#                     feature_info = extract_feature_info(feature)
#                     csvwriter.writerow([bgc_locus, bgc_definition] + feature_info)





# from Bio import SeqIO
# import csv

# # Function to extract information from a feature
# def extract_feature_info(feature):
#     start = feature.location.start
#     end = feature.location.end
#     length = len(feature)
#     product = feature.qualifiers.get('product', [''])[0]
#     locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
#     note = feature.qualifiers.get('note', [''])[0]
#     gene_functions = feature.qualifiers.get('gene_functions', [''])[0]
#     description = feature.qualifiers.get('description', [''])[0]
#     gene_ontologies = feature.qualifiers.get('gene_ontologies', [''])[0]

#     return [start, end, length, product, locus_tag, note, gene_functions, description, gene_ontologies]

# # Open the output CSV file
# with open('output.csv', 'w', newline='') as csvfile:
#     csvwriter = csv.writer(csvfile)
    
#     # Write header row
#     csvwriter.writerow(['BGC Locus', 'BGC Definition', 'Start', 'End', 'Length', 'Product', 'Locus Tag', 'Note', 'Gene Functions', 'Description', 'Gene Ontologies'])
    
#     # Parse the .gbk file
#     for record in SeqIO.parse("/geomicro/data21/lnhart/Projects/Kenya_2023/flagship_manuscript/FINAL/data/bgc_mining/data/quality_bins/drep_99/antismash_bins99/antismash_drep_BGCs/extractedbgcs_95_90_rep_seq/samp_4434_metabat2_58_contig_50_T3PKS.region001.gbk", "genbank"):
#         bgc_locus = record.name
#         bgc_definition = record.description
        
#         # Keep track of previously encountered locus tags
#         previous_locus_tags = set()
        
#         # Iterate over features in the record
#         for feature in record.features:
#             if feature.type == 'CDS':  # Only process CDS features
#                 locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                
#                 # Check if locus_tag has been encountered before
#                 if locus_tag not in previous_locus_tags:
#                     feature_info = extract_feature_info(feature)
#                     csvwriter.writerow([bgc_locus, bgc_definition] + feature_info)
#                     previous_locus_tags.add(locus_tag)





