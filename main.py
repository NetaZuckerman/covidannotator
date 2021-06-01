#!/usr/bin/python
import pandas as pd
import covidAnnotator
import sys
import argparse


def main():
    all_mutations = pd.read_csv("all_mutations.csv")
    b117muts = pd.read_csv("b117muts.csv")
    uniqueIDs = pd.unique(all_mutations["Sequence ID"])
    uniqueMuts = pd.unique(all_mutations["nuc name"])
    numofSequences = len(uniqueIDs)
    mutNucCount = all_mutations['nuc name'].value_counts()
    freqs = (mutNucCount / numofSequences * 100)
    freqs = freqs.T.to_dict()

    freqTable = pd.DataFrame()
    freqTable['nuc name'] = uniqueMuts
    mapping = dict(all_mutations[['nuc name', 'type']].values)
    freqTable['type'] = freqTable['nuc name'].map(mapping)
    mapping = dict(all_mutations[['nuc name', 'protein']].values)
    freqTable['protein'] = freqTable['nuc name'].map(mapping)
    mapping = dict(all_mutations[['nuc name', 'AAMutation']].values)
    freqTable['AAMutation'] = freqTable['nuc name'].map(mapping)
    count_title = 'Count (' + str(numofSequences) + ')'
    freqTable[count_title] = freqTable['nuc name'].map(all_mutations['nuc name'].value_counts())
    freqTable['Freq'] = freqTable['nuc name'].map(freqs)
    freqTable.sort_values(by=['protein', 'Freq'], ascending=False, inplace=True)
    freqTable = freqTable.loc[freqTable['Freq'] >= 2]
    mapping = dict(b117muts[['nucleotide', 'lineage original']].values)
    freqTable['isUKLineage'] = freqTable['nuc name'].map(mapping)
    num_muts = len(freqTable.index) - 1
    if num_muts <= 0:
        num_muts = 0
    geneFreq = freqTable['protein'].value_counts()  # + " \ " + str(num_muts)
    geneFreq['total'] = num_muts
    geneFreq = pd.DataFrame(geneFreq)
    if not geneFreq.empty:
        geneFreq['Freq'] = geneFreq / num_muts * 100

    writer = pd.ExcelWriter("Freq_Table" + ".xlsx", engine='xlsxwriter')
    freqTable.to_excel(writer, sheet_name='Mutations Frequencies', index=False)
    geneFreq.to_excel(writer, sheet_name='gene count')
    writer.save()
    print("Freq_Table" + ".csv is ready")


if __name__ == '__main__':
    #  userChoose = loop.main()
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    parser.add_argument("-i", help="insertions files path", nargs=1, type=str,
                        dest='insPath')  # dest = name of variable
    parser.add_argument("-n", help="to add 'N' to excel", action='store_true')  # dest = name of variable
    args = parser.parse_args()
    covidAnnotator.main(args)
    main()
