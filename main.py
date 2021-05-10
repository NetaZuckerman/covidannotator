import csv

import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import re
from math import ceil
import pandas as pd
import datetime

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
        if filteredregionTitles:
            region_id = str([*filteredregionTitles][0])
        else:
            filteredregionTitles["none"] = [0, 0]
            region_id = str([*filteredregionTitles][0])
        return region_id, filteredregionTitles[region_id][0], filteredregionTitles[region_id][1]


def getTranslate(i, regionsList, referenceSequence, record, flag):
    regionTitle, regionStart, regionEnd = getRegion(i, regionsList)
    RefAA = translate(''.join(map(str, referenceSequence)), regionStart - 1, regionEnd,
                      codon_map)
    OtherSeqAA = translate(str(record.seq), regionStart - 1, regionEnd, codon_map)
    aaMutIndexr = (i + 1 - regionStart) / 3
    if aaMutIndexr % 1 == 0:
        aaMutIndex = int(aaMutIndexr + 1)
    else:
        aaMutIndex = ceil((i + 1 - regionStart) / 3)
    refaaseq = RefAA[aaMutIndex - 2:aaMutIndex + 2]
    otheraaseq = OtherSeqAA[aaMutIndex - 2:aaMutIndex + 2]
    if refaaseq:
        try:
            AAMutToCSv = str(RefAA[aaMutIndex - 1] + str(aaMutIndex) + OtherSeqAA[aaMutIndex - 1])
        except:
            print(str(aaMutIndex) + " " + str(len(RefAA)) + " " + str(len(OtherSeqAA)))
            AAMutToCSv = "None"

        if flag == 2:
            AAMutToCSv = AAMutToCSv[:-1] + "Del"
    else:
        AAMutToCSv = "None"
    return regionTitle, AAMutToCSv


def writeToCSV(writer, record, refnuc, mutnuc, i, regionTitle, AAMutToCSv, record_date, type):
    try:
        writer.writerow({'Sequence ID': record.id, 'Reference Nucleotide': refnuc,
                         'Mutation nucleotide': mutnuc, 'location': i + 1,
                         'nuc name': str(i + 1) + " " + refnuc + " -> " + mutnuc,
                         'protein': str(regionTitle), 'AAMutation': AAMutToCSv,
                         'varname': str(regionTitle) + ":" + AAMutToCSv, 'date': str(record_date), 'type': type})
    except:
        writer.writerow({'Sequence ID': record, 'Reference Nucleotide': refnuc,
                         'Mutation nucleotide': mutnuc, 'location': i + 1,
                         'nuc name': str(i + 1) + " " + refnuc + " -> " + mutnuc,
                         'protein': str(regionTitle), 'AAMutation': AAMutToCSv,
                         'varname': str(regionTitle) + ":" + AAMutToCSv, 'date': str(record_date), 'type': type})


def getDate(metadata, id):
    try:
        record_date = datetime.datetime.strptime(metadata.at[id, 'date'], '%Y-%m-%d').date()
    except:
        print(id)
        record_date = datetime.datetime(1111, 1, 1).date()
    return record_date


def filterByDate(record_date,argv):
    date1 = datetime.datetime(2021, int(argv[2]), int(argv[4])).date()
    date2 = datetime.datetime(2021, int(argv[2]), int(argv[3])).date()
    return date1 >= record_date >= date2


def addInsertions(argv, writer, regionList, ref, metadata):
    df = pd.read_csv(argv[5])
    for rowindex, row in df.iterrows():
        record_date = getDate(metadata, df.at[rowindex, 'strain'])
        if filterByDate(record_date,argv):
            for colindex, cell in enumerate(row):
                if cell is not np.nan:
                    location = df.columns.values[colindex]
                    try:
                        region = getRegion(int(location), regionList)
                        writeToCSV(writer, df.at[rowindex, 'strain'], ref[int(location) - 1], cell, int(location) - 1,
                                   region[0], "Ins", record_date, "Insertion")
                    except:
                        pass


def main(argv):
    month=argv[1]
    print("Starting...")
    # importing the Regions list
    metadata = pd.read_csv("metadata.tsv", sep='\t')
    metadata.set_index('strain', inplace=True)
    muttable = pd.read_csv("novelMutTable.csv")
    muttable = muttable.drop(muttable[muttable['type'] == 'Insertion'].index)
    all_mutations = set([x for x in muttable.AA])
    regiontable = "regions.csv"
    with open(regiontable, 'r') as f:
        reader = csv.reader(f)
        my_list = []
        for row in reader:
            my_list.append({'segment': row[0], 'id': row[1],
                            'region': row[2], 'start': row[3], 'end': row[4], 'function': row[5]})
        regionsList = my_list[1:]

    fieldnames = ['Sequence ID', 'Reference Nucleotide', 'Mutation nucleotide', 'location', 'nuc name', 'protein',
                  'AAMutation', 'varname', 'date', 'type']
    with open(f'unknownedMutations{month}.csv', 'w', newline='') as csvfile1:
        writer2 = csv.DictWriter(csvfile1, fieldnames=fieldnames)
        # Creating the output csv file
        with open(f'uk_{month}.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer2.writeheader()
            # parsing the pasta file
            records = list(SeqIO.parse(argv[0], "fasta"))
            # the first sequence is the reference
            referenceSequence = list(records[0])
            otherSequences = records[1:]
            # referenceSequence = list(SeqIO.parse("REF_NC_045512.2.fasta", "fasta"))[0].seq
            # otherSequences = list(SeqIO.parse("9747.fa", "fasta"))
            # Find mutations = iterating over each sequence and compare to the reference sequence
            for record in otherSequences:
                record_date = getDate(metadata, record.id)
                if filterByDate(record_date,argv):
                    for i, nucleotide in enumerate((record)):
                        if nucleotide != "-" and nucleotide != "N" and referenceSequence[i] != nucleotide:
                            regionTitle, AAMutToCSv = getTranslate(i, regionsList, referenceSequence, record, 1)
                            # Each mutation is written to the output csv file
                            writeToCSV(writer, record, referenceSequence[i], nucleotide, i, regionTitle, AAMutToCSv,
                                       record_date, "SNP")
                            NonSyn = AAMutToCSv[0] != AAMutToCSv[len(AAMutToCSv) - 1]
                            if AAMutToCSv not in all_mutations and NonSyn:
                                writeToCSV(writer2, record, referenceSequence[i], nucleotide, i, regionTitle,
                                           AAMutToCSv,
                                           record_date, "SNP")
                        elif nucleotide == "-" and i > 100:
                            regionTitle, AAMutToCSv = getTranslate(i, regionsList, referenceSequence, record, 2)
                            writeToCSV(writer, record, referenceSequence[i], nucleotide, i, regionTitle, AAMutToCSv,
                                       record_date, "Del")
            if argv[5]:
                addInsertions(argv, writer, regionsList, referenceSequence, metadata)


if __name__ == "__main__":
    main(sys.argv[1:])
