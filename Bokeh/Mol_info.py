import pubchempy as pcp
import ssl

# Disable SSL certificate verification
ssl._create_default_https_context = ssl._create_unverified_context

def mol_info(smi):
    # Retrieve the compound using SMILES notation
    compound = pcp.get_compounds(smi, 'smiles')
    # Extract PubChem Compound ID (CID) from the compound
    pubchem_cid = compound[0].cid
    # Get properties including identifiers for the compound
    properties = pcp.get_properties(pubchem_cid, 'cid,canonical_smiles,inchi,inchikey')
    return properties

input_smi = input('SMILES: ')
print(mol_info(input_smi))

