import os
import argparse
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
import re
import pandas as pd


mod_motifs = {
    'm6A': [
        'AAACA',
        'AAACC',
        'AAACT',
        'AGACA',
        'AGACC',
        'AGACT',
        'GAACA',
        'GAACC',
        'GAACT',
        'GGACA',
        'GGACC',
        'GGACT',
        'TAACA',
        'TAACC',
        'TAACT',
        'TGACA',
        'TGACC',
        'TGACT'
    ],
    'psi': [
        'AGTGG',
        'ATTTG',
        'CATAA',
        'CATCC',
        'CCTCC',
        'CTTTA',
        'GATGC',
        'GGTCC',
        'GGTGG',
        'GTTCA',
        'GTTCC',
        'GTTCG',
        'GTTCT',
        'TATAA',
        'TGTAG',
        'TGTGG',
    ],
    # 'Gm': [
    #     "ACGTC",
    #     "ATGAC",
    #     "ATGCT",
    #     "ATGTC",
    #     "ATGTG",
    #     "ATGTT",
    #     "CCGCC",
    #     "CTGCC",
    #     "CTGCG",
    #     "CTGCT",
    #     "CTGTA",
    #     "CTGTC",
    #     "CTGTG",
    #     "GCGCC",
    #     "GTGCA",
    #     "GTGCC",
    #     "GTGTC",
    #     "TCGCC"
    # ],
}


fields = [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'ref5mer'
]

parser = argparse.ArgumentParser(
    prog='Mod sites generator',
    description='Given an input reference and modification types, generate bed file annotating all possible sites',
)
parser.add_argument('--ref_file', required=True)
parser.add_argument('--out_dir', required=True)
parser.add_argument('--chroms', default=None, help='chromosomes to include from reference, comma-separated')
parser.add_argument('--mods', default=None, help='mods to include, {m6A,psi},  comma-separated')
parser.add_argument('--out_prefix', default=None)


def get_reference(ref_file):
    print(f'Loading reference from {ref_file}...')
    out_ref = {}
    for record in tqdm(SeqIO.parse(ref_file, 'fasta')):
        out_ref[record.id] = str(record.seq)
    return out_ref


def main():
    args = parser.parse_args()

    ref = get_reference(args.ref_file)

    os.makedirs(args.out_dir, exist_ok=True)

    if args.chroms is not None:
        all_chroms = args.chroms.split(',')
    else:
        all_chroms = [str(x) for x in ref.keys()]

    if args.mods is not None:
        all_mods = args.mods.split(',')
    else:
        all_mods = list(mod_motifs.keys())

    if args.out_prefix is not None:
        out_prefix = args.out_prefix
    else:
        ref_name = os.path.basename(args.ref_file).rstrip('.fasta')
        out_prefix = '.'.join(all_mods) + '.' + ref_name

    for this_chrom in all_chroms:
        print(f'Parsing chr{this_chrom}')
        seq_forward = ref[this_chrom]
        seq_revcomp = str(Seq(ref[this_chrom]).reverse_complement())
        seq_len = len(seq_forward)
        chromStart_name_strand_ref5mer = []
        for this_mod, this_mod_motifs in mod_motifs.items():
            for this_motif in this_mod_motifs:
                for this_match in re.finditer(this_motif, seq_forward):
                    chromStart_name_strand_ref5mer.append((this_match.start()+2, this_mod, '+', this_motif))
                for this_match in re.finditer(this_motif, seq_revcomp):
                    chromStart_name_strand_ref5mer.append((seq_len-this_match.start()-3, this_mod, '-', this_motif))
        chromStart_name_strand_ref5mer.sort()
        this_chrom_df = pd.DataFrame(chromStart_name_strand_ref5mer, columns=['chromStart', 'name', 'strand', 'ref5mer'])
        this_chrom_df['chrom'] = this_chrom
        this_chrom_df['chromEnd'] = this_chrom_df['chromStart'] + 1
        this_chrom_df['score'] = '.'
        this_chrom_df = this_chrom_df[fields]

        print(f'Writing out {len(this_chrom_df)} sites')
        this_chrom_df.to_csv(os.path.join(args.out_dir, f'{out_prefix}.chr{this_chrom}.bed'),
                             sep='\t', index=False, header=True)


if __name__ == "__main__":
    main()
