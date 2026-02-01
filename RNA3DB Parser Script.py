import os
import json
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import three_to_one

# Path to your extracted RNA3DB mmCIF files
CIF_DIR = "rna3db-mmcifs.v2"
JSON_DIR = "rna3db-jsons"

# Output file
OUTPUT_CSV = "rna3db_processed.csv"

def extract_rna_sequence_and_coords(cif_path):
    """Extracts sequence and 3D coordinates from an mmCIF RNA file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("RNA", cif_path)
    
    sequence = ""
    coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() in ['A', 'U', 'G', 'C']:  # Only canonical RNA bases
                    try:
                        nt = residue.get_resname().strip()
                        sequence += nt
                        # Extracting phosphate atom (P) coordinates
                        if "P" in residue:
                            atom = residue["P"]
                            coords.append(atom.coord)
                    except:
                        continue
            break  # Only one chain per file

    return sequence, np.array(coords)


def load_metadata(json_dir):
    """Optional: Load Rfam + split info from RNA3DB JSONs."""
    all_meta = []
    for fname in os.listdir(json_dir):
        if fname.endswith(".json"):
            with open(os.path.join(json_dir, fname)) as f:
                data = json.load(f)
                all_meta.append(data)
    return all_meta


def main():
    records = []
    for root, dirs, files in os.walk(CIF_DIR):
        for file in files:
            if file.endswith(".cif"):
                file_path = os.path.join(root, file)
                try:
                    seq, coords = extract_rna_sequence_and_coords(file_path)
                    if len(seq) > 10 and coords.shape[0] == len(seq):  # quality check
                        records.append({
                            "filename": file,
                            "sequence": seq,
                            "num_nt": len(seq),
                            "coords": coords.tolist()  # You can also save as npz
                        })
                except Exception as e:
                    print(f"Error parsing {file}: {e}")
                    continue

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\n✅ Processed {len(df)} RNAs into {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
