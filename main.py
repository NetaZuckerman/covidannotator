#!/usr/bin/python
import pandas as pd
import covidAnnotator
import sys

def main():
    all_mutations = pd.read_csv("all_mutations.csv")
    uniqueIDs=pd.unique(all_mutations["Sequence ID"])
    uniqueMuts=pd.unique(all_mutations["nuc name"])
    numofSequences=len(uniqueIDs)
    mutNucCount=all_mutations['nuc name'].value_counts()
    freqs=(mutNucCount/ numofSequences*100)
    freqs=freqs.T.to_dict()

    freqTable=pd.DataFrame()
    freqTable['nuc name']=uniqueMuts
    mapping = dict(all_mutations[['nuc name', 'type']].values)
    freqTable['type'] = freqTable['nuc name'].map(mapping)
    mapping = dict(all_mutations[['nuc name', 'protein']].values)
    freqTable['protein'] = freqTable['nuc name'].map(mapping)
    mapping = dict(all_mutations[['nuc name', 'AAMutation']].values)
    freqTable['AAMutation'] = freqTable['nuc name'].map(mapping)
    count_title='Count ('+str(numofSequences)+')'
    freqTable[count_title]=freqTable['nuc name'].map(all_mutations['nuc name'].value_counts())
    freqTable['Freq'] = freqTable['nuc name'].map(freqs)
    freqTable.sort_values(by=['protein', 'Freq'],ascending=False,inplace=True)
    freqTable=freqTable.loc[freqTable['Freq'] >= 2]

    freqTable.to_csv("Freq_Table.csv", index=False)

    print("Freq_Table.csv is ready")



if __name__ == '__main__':
    covidAnnotator.main(sys.argv[1:])
    main()
