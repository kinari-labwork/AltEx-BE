import argparse
from pathlib import Path

def parse_path_from_args(args: argparse.Namespace):
    refflat_path = Path(args.refflat_path)
    fasta_path = Path(args.fasta_path)
    output_directory = Path(args.output_dir)
    gene_file = Path(args.gene_file) if args.gene_file else None
    return refflat_path, fasta_path, output_directory, gene_file
