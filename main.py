import difflib
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import numpy as np
import re
from math import ceil

codon_map = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def full_codon_gaps(sequence, start, end, gap='-'):
    """
    avoid partial gaps in codon and convert to whole gaps codon
    exmple: from A-- to ---
    :param sequence: Seq object (Biopython)
    :param start: start of reading frame
    :param end: end of reading frame
    :param gap: gap character. default: '-'
    :return: new sequence, with full codons
    """
    old_seq = str(sequence)
    new_seq = old_seq[:start]
    for i in range(start - 1, end, 3):
        codon = old_seq[i: i + 3]
        if '-' in codon:
            codon = '---'
        new_seq += codon
    new_seq += old_seq[end:]
    return Seq.Seq(new_seq)


# get start and end of region
# don't forget -1 in position
def translate(sequence, start, end, codon_table):
    """
    translate nucleotides sequence in given region to amino acid sequence according to codons in region start-end
    :param sequence: nucleotides sequence as str #
    :param start: position of first nucleotide
    :param end: position of last nucleotide
    :return: translated sequence (aa seq)
    """
    tranlsated = []
    for i in range(start, end, 3):
        codon = sequence[i:i + 3]
        if codon in codon_table:
            aa = codon_table[codon]
        else:
            aa = 'X'  # ignore frameshift
        tranlsated.append(aa)
    return tranlsated


def getRegion(i, regionsList):
    regionTitles = {}
    for region in (regionsList):
        if i + 1 in range(int(region["start"]), int(region["end"])):
            regionTitles[region["id"]] = [int(region["start"]), int(region["end"])]
        if i + 1 < int(region["start"]):
            break
    if len(regionTitles) > 1:
        try:
            regex = re.compile("^(?!ORF)")
            filteredregionTitles = {k: v for k, v in regionTitles.items() if not k.startswith('ORF')}
            region_id = str([*filteredregionTitles][0])
            return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]
        except:
            filteredregionTitles = {k: v for k, v in regionTitles.items()}
            region_id = str([*filteredregionTitles][0])
            return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]
    else:
        filteredregionTitles = {k: v for k, v in regionTitles.items()}
        region_id = str([*filteredregionTitles][0])
        return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]


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
        fieldnames = ['Sequence ID', 'Reference Nucleotide', 'Mutation nucleotide', 'location', 'nuc name', 'protein',
                      'AAMutation', 'varname']
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
                if nucleotide != "-" and nucleotide != "N" and referenceSequence[i] != nucleotide:
                    regionTitle, regionStart, regionEnd = getRegion(i, regionsList)
                    RefAA = translate(''.join(map(str, referenceSequence)), regionStart - 1, regionEnd, codon_map)
                    OtherSeqAA = translate(str(record.seq), regionStart - 1, regionEnd, codon_map)
                    aaMutIndexr = (i + 1 - regionStart) / 3
                    if aaMutIndexr % 1 == 0:
                        aaMutIndex = int(aaMutIndexr + 1)
                    else:
                        aaMutIndex = ceil((i + 1 - regionStart) / 3)
                    refaaseq = RefAA[aaMutIndex - 2:aaMutIndex + 2]
                    otheraaseq = OtherSeqAA[aaMutIndex - 2:aaMutIndex + 2]
                    AAMutToCSv = str(RefAA[aaMutIndex - 1] + str(aaMutIndex) + OtherSeqAA[aaMutIndex - 1])
                    # Each mutation is written to the output csv file
                    writer.writerow({'Sequence ID': record.id, 'Reference Nucleotide': referenceSequence[i],
                                     'Mutation nucleotide': nucleotide, 'location': i + 1,
                                     'nuc name': str(i + 1) + " " + referenceSequence[i] + " -> " + nucleotide,
                                     'protein': str(regionTitle), 'AAMutation': AAMutToCSv,
                                     'varname': str(regionTitle) + ":" + AAMutToCSv})


if __name__ == "__main__":
    main()
