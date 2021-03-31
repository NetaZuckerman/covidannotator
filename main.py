import difflib
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import numpy as np


def getRegion(i,regionsList):
    regionTitles = []
    for region in regionsList:
        if i + 1 in range(int(region["start"]), int(region["end"])):
            regionTitles.append(region["id"])
    if len(regionTitles) > 1:
        return str(regionTitles[len(regionTitles)-1])
    else: return str(regionTitles[0])

def main():
    # importing the Regions list
    with open('regions.csv', 'r') as f:
        reader = csv.reader(f)
        my_list = []
        for row in reader:
            my_list.append({'segment': row[0], 'id': row[1],
                            'region': row[2], 'start': row[3], 'end': row[4], 'function': row[5]})
        regionsList = my_list[1:]

    # Creating the output csv file
    with open('output.csv', 'w', newline='') as csvfile:
        fieldnames = ['Sequence ID', 'Reference Nucleotide', 'Mutation nucleotide', 'location', 'nuc name','protein']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        # parsing the pasta file
        records = list(SeqIO.parse("FrenchVar_aligned.fasta", "fasta"))
        # the first sequence is the reference
        referenceSequence = list(records[0])
        otherSequences = records[1:]
        # Find mutations = iterating over each sequence and compare to the reference sequence
        for record in otherSequences:
            a = record.id
            a = list(record)
            for i, nucleotide in enumerate((record)):
                if nucleotide != "-" and referenceSequence[i] != nucleotide:
                    regionTitle=getRegion(i,regionsList)
                    # Each mutation is written to the output csv file
                    writer.writerow({'Sequence ID': record.id, 'Reference Nucleotide': referenceSequence[i],
                                     'Mutation nucleotide': nucleotide, 'location': i + 1,
                                     'nuc name': str(i + 1) + " " + referenceSequence[i] + " -> " + nucleotide,'protein':str(regionTitle)})


if __name__ == "__main__":
    main()
