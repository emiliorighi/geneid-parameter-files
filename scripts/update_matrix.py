import os
import pandas as pd
import requests
from lxml import etree
import time


MATRIX_PATH="matrix.tsv"
PARAM_FILES_PATH="parameter_files"
NCBI_RESPONSE_MAPPER = {
    'lineage': 'taxidlist',
    'tax_id':'taxid',
    'organism_name':'species'
}

def get_taxon_from_ENA(taxid):
    taxon = dict()
    try:
        response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{taxid}")
        xml = response.content
        if not xml or len(xml) == 0:
            return taxon
        root = etree.fromstring(xml)
        for child in root:
            child_taxid = int(child.attrib["taxId"])
            if(child_taxid == taxid):
                lineage = []
                for taxon_node in child:
                    if taxon_node.tag == 'lineage':
                        lineage.append(int(child.attrib["taxId"]))
                        lineage.extend([int(node.attrib['taxId']) for node in taxon_node])
                        taxon['lineage'] = lineage
                        taxon['taxid'] = child_taxid
                        taxon['species'] = child.attrib['scientificName']
                        break
        return taxon
    except:
        return taxon
    
def get_taxon_from_NCBI(taxid):
    taxon = dict()
    try:
        response = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v1/taxonomy/taxon/{taxid}")
        if not response.json() or not 'taxonomy_nodes' in response.json():
            return taxon
        taxon_to_insert = response.json()['taxonomy_nodes'][0]['taxonomy']
        for k in NCBI_RESPONSE_MAPPER.keys():
            value = NCBI_RESPONSE_MAPPER[k]
            taxon[value] = taxon_to_insert[k]
        return taxon
    except:
        return taxon

def append_row(df, row):
    return pd.concat([
                df, 
                pd.DataFrame([row], columns=row.index)]
           ).reset_index(drop=True)

param_files = [entry.name for entry in os.scandir(PARAM_FILES_PATH) if '.param' in entry.name]

taxids = [int(name.split('.')[1]) for name in param_files]

table = pd.read_table(MATRIX_PATH)

columns = list(table.columns)

existing_taxids = table.loc[:,"taxid"].tolist()

new_taxids = [int(new_taxid) for new_taxid in taxids if int(new_taxid) not in existing_taxids]

## counter to wait 1 second after 2 requests
counter = 0

for new_taxid in new_taxids:
    if(counter > 2):
        time.sleep(1)
    taxon = get_taxon_from_NCBI(new_taxid)
    if not taxon:
        taxon = get_taxon_from_ENA(new_taxid)
        if not taxon:
            counter += 1
            #skip taxon if does not exists in ncbi
            continue
    taxon['lineage'] = taxon['lineage'].reverse()
    counter += 1
    param_file_index = taxids.index(taxon['taxid'])
    taxon['parameter_file'] = param_files[param_file_index]
    new_row = pd.Series(taxon)
    table = append_row(table, new_row)

table.to_csv(MATRIX_PATH, sep = "\t", index = False)