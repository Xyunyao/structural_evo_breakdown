import argparse
from dms_utils import deep_mutational_scan
from pathlib import Path
import numpy as np

#import esm
#from util import load_structure, extract_coords_from_structure
#import biotite.structure
from collections import defaultdict

# Extract the native sequence from the fasta fiel
# def get_native_seq(pdbfile, chain):
#     structure = load_structure(pdbfile, chain)
#     _ , native_seq = extract_coords_from_structure(structure)
#     return native_seq
def read_fasta(file_path):
    '''Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values.'''
    sequences = {}
    with open(file_path, 'r') as f:
        seq_id = None
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    sequences[seq_id] = seq
                seq_id = line[1:5]  # Remove '>'
                seq = ''
            else:
                seq += line
        if seq_id is not None:
            sequences[seq_id] = seq  # Add the last sequence
    return seq

def write_dms_lib(args):
    '''Writes a deep mutational scanning library, including the native/wildtype (wt) of the 
    indicated target chain in the structure to an output Fasta file'''

    sequence = read_fasta(args.fastafile)
    Path(args.dmspath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.dmspath, 'w') as f:
        f.write('>wt\n')
        f.write(sequence+'\n')
        for pos, wt, mt in deep_mutational_scan(sequence):
            assert(sequence[pos] == wt)
            mut_seq = sequence[:pos] + mt + sequence[(pos + 1):]

            f.write('>' + str(wt) + str(pos+1) + str(mt) + '\n')
            f.write(mut_seq + '\n')


def main():
    parser = argparse.ArgumentParser(
            description='Create a DMS library based on full protein sequence.'
    )
    parser.add_argument(
            '--protein_ID', type=str,
            help='protein name',
    )
    parser.add_argument(
            '--fastafile', type=str,
            help='input filepath, .fasta',
    )
    parser.add_argument(
            '--dmspath', type=str,
            help='output filepath for dms library',
    )
   

    args = parser.parse_args()

    if args.dmspath is None:
        args.dmspath = f'predictions/{args.protein_ID}_dms.fasta'

    write_dms_lib(args)

if __name__ == '__main__':
    main()