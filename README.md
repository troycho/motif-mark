# Motif-Mark

## Overview
The goal of this script is to use object-oriented programming in Python to readily visualize predefined motif binding sites in nucleotide sequences from a FASTA file. The input to the program is a FASTA file with up to 10 records and a file with up to 5 motif sequences (which can include ambiguous bases). The script will output one PNG figure per FASTA file that visualizes the locations of all the provided motifs in each gene sequence. The horizontal lengths of all the genetic features in the figure are to scale. 

## Usage
The required argparse options are `-f` for the path to the input FASTA file and `-m` for the path to the motifs file. Each motif sequence in the latter file should be on its own line. The gene sequences should be up to 1000 bases in length. When the script is run, it will output a PNG file with the same prefix as the input FASTA file. 

## Files
Script: [motif-mark-oop.py](./motif-mark-oop.py)  

Example PNG output: [Figure_1.png](./Figure_1.png)
