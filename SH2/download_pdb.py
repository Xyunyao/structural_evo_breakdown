import argparse
import requests

def download_pdb(pdb_id):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)
    
    if response.status_code == 200:
        with open(f"{pdb_id}.pdb", "w") as file:
            file.write(response.text)
        print(f"PDB file '{pdb_id}.pdb' downloaded successfully.")
    else:
        print(f"Failed to download PDB file '{pdb_id}.pdb'. Please check the PDB ID.")

def download_fasta(pdb_id):
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    response = requests.get(fasta_url)
    
    if response.status_code == 200:
        with open(f"{pdb_id}.fasta", "w") as file:
            file.write(response.text)
        print(f"FASTA file '{pdb_id}.fasta' downloaded successfully.")
    else:
        print(f"Failed to download FASTA file '{pdb_id}.fasta'. Please check the PDB ID.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download PDB and FASTA files by PDB ID")
    parser.add_argument("pdb_id", help="PDB ID of the files to download")
    
    args = parser.parse_args()
    download_pdb(args.pdb_id)
    download_fasta(args.pdb_id)
