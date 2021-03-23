import difflib
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import numpy as np


def main():
    with open('names.csv', 'w', newline='') as csvfile:
        fieldnames = ['Sequence ID', 'Reference Nucleotide','Mutation nucleotide','location']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        records = list(SeqIO.parse("names.fasta", "fasta"))
        a=records[0].id
        referenceSequence = list(records[0])
        otherSequences = records[1:]
        for record in otherSequences:
            a=record.id
            a=list(record)
            for i,nucleotide in enumerate((record)):
                if nucleotide != "-" and referenceSequence[i] != nucleotide:
                    writer.writerow({'Sequence ID': record.id, 'Reference Nucleotide': referenceSequence[i],'Mutation nucleotide':nucleotide,'location':i})


if __name__ == "__main__":
    main()
