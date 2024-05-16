from __future__ import division
import os
import re

def search_pdb_file(uniprot_id, folder_path):
    pattern = re.compile(r"AF-{}-F\d+-model_v4.pdb".format(uniprot_id))
    for subdir, _, files in os.walk(folder_path):
        if pattern.match(os.path.basename(subdir)):
            pdb_files = [f for f in files if f.endswith('.pdb')]
            if pdb_files:
                return os.path.join(subdir, pdb_files[0])
            else:
                return None  # Return None if PDB file not found

pdb_file_path = search_pdb_file(uniprot_id, r"C:\Users\Sabrina\PycharmProjects\intrinsic_disorder\proteome_human")
