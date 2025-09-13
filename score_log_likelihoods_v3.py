import argparse
from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
import torch.nn.functional as F
from tqdm import tqdm
import warnings
import os
from dms_utils import deep_mutational_scan
import esm
import pandas as pd
from multichain_util import extract_coords_from_complex, _concatenate_coords, _concatenate_seqs, score_sequence_in_complex
from util import get_sequence_loss, load_structure, load_coords, score_sequence, extract_coords_from_structure
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def get_native_seq(pdbfile, chain):
        structure = load_structure(pdbfile, chain)
        _ , native_seq = extract_coords_from_structure(structure)
        return native_seq


def align_sequences(seq1, seq2):
    # Perform global alignment using the Needleman-Wunsch algorithm
    alignments = pairwise2.align.globalxx(seq1, seq2)

    # Take the best alignment (highest score)
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, start, end = best_alignment

    # Generate the mask
    mask = []
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b and a != "-":
            mask.append(1)  # Match
        else:
            mask.append(0)  # Mismatch or gap

    return mask
#
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

def updata_cords(coords, mask,start_index=0):
    """
    Update the coordinates based on the mask.
    Args:
        coords: numpy array of shape (N, 3) where N is the number of atoms.
        mask: list of 0s and 1s indicating which atoms to keep.
    Returns:
        Updated coordinates as a numpy array.
    """
    updated_coords= np.full((len(mask), 3, 3), np.inf)

    # Fill the result array with coord values based on the mask
    coord_index = start_index
    for i, value in enumerate(mask):
        if value == 1:
            updated_coords[i] = coords[coord_index]
            coord_index += 1
    return np.array(updated_coords)

def score_singlechain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    coords, coords_seq = load_coords(args.pdbfile, args.chain)

    # code to modify the native sequenca and coords
    native_seq = ''.join(read_fasta(args.seqfile_path).split())
    mask= align_sequences(coords_seq, native_seq)
    coords = updata_cords(coords, mask)
    coords = coords[args.start_index:]
    # make sure the coords and native_seq are aligned
    assert len(coords) == len(native_seq), "Coords and native sequence lengths do not match"
   # convert coords to float32
    coords = coords.astype(np.float32)
    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')
    ll, _ = score_sequence(
            model, alphabet, coords, native_seq)
    print('Native sequence')
    print(f'Log likelihood: {ll:.2f}')
    print(f'Perplexity: {np.exp(-ll):.2f}')
    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqpath)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood\n')
        for header, seq in tqdm(seqs.items()):
            ll, _ = score_sequence(
                    model, alphabet, coords, str(seq))
            fout.write(header + ',' + str(ll) + '\n')
    print(f'Results saved to {args.outpath}') 


def score_multichain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    structure = load_structure(args.pdbfile)
    coords, native_seqs = extract_coords_from_complex(structure)
    target_chain_id = args.chain
    native_seq = native_seqs[target_chain_id]
    order = args.order

    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')

    ll_complex, ll_targetchain = score_sequence_in_complex(
        model,
        alphabet,
        coords,
        native_seqs,
        target_chain_id,
        native_seq,
        order=order,
    )
    print('Native sequence')
    print(f'Log likelihood of complex: {ll_complex:.2f}')
    print(f'Log likelihood of target chain: {ll_targetchain:.2f}')
    print(f'Perplexity: {np.exp(ll_complex):.2f}')

    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqpath)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood, log_likelihood_target\n')
        for header, seq in tqdm(seqs.items()):
            ll_complex, ll_targetchain = score_sequence_in_complex(
                model,
                alphabet,
                coords,
                native_seqs,
                target_chain_id,
                str(seq),
                order=order,
            )
            fout.write(header + ',' + str(ll_complex) + ',' + str(ll_targetchain) + '\n')
    print(f'Results saved to {args.outpath}') 

def get_model_checkpoint_path(filename):
    # Expanding the user's home directory
    return os.path.expanduser(f"~/.cache/torch/hub/checkpoints/{filename}")

def main():
    parser = argparse.ArgumentParser(
            description='Score sequences based on a given structure.'
    )
    parser.add_argument(
            '--pdbfile', type=str,
            help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
            '--seqfile_path', type=str,
            help='input filepath,.fasta file containing the native sequence',
    )
    parser.add_argument(
            '--start_index', type= int, default=0,
            help='int, the start index for coords, default is 0',
    )
    parser.add_argument(
            '--seqpath', type=str,
            help='input filepath for variant sequences in a .fasta file',
    )
    parser.add_argument(
            '--outpath', type=str,
            help='output filepath for scores of variant sequences',
    )
    parser.add_argument(
            '--chain', type=str,
            help='chain id for the chain of interest', default='A',
    )
    parser.set_defaults(multichain_backbone=True)
    parser.add_argument(
            '--multichain-backbone', action='store_true',
            help='use the backbones of all chains in the input for conditioning'
    )
    parser.add_argument(
            '--order', type=str, default=None,
            help='for multichain, specify order'
    )
    parser.add_argument(
            '--singlechain-backbone', dest='multichain_backbone',
            action='store_false',
            help='use the backbone of only target chain in the input for conditioning'
    )

    parser.add_argument(
        "--nogpu", action="store_true", 
        help="Do not use GPU even if available"
    )
    args = parser.parse_args()

    if args.outpath is None:
        args.outpath = f'output/{args.pdbfile[:-4]}-chain{args.chain}_scores.csv'

    model_checkpoint_path = get_model_checkpoint_path('esm_if1_20220410.pt')
    with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            model, alphabet = esm.pretrained.load_model_and_alphabet( \
                model_checkpoint_path \
            )
    model = model.eval()


    if args.multichain_backbone:
        score_multichain_backbone(model, alphabet, args)
    else:
        score_singlechain_backbone(model, alphabet, args)


if __name__ == '__main__':
    main()
