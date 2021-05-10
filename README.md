# covidannotator

This program find all the mutations in given sequences filtered by a given date

First activate conda - > conda activate nextstrain (or other env with pandas)

Input - 
in the wd folder have to be - "metadatat.tsv" "novelMutable.csv" "regions.csv"
as input in the CLI ->
1) the Fasta file, the first sequence must be the reference, other sequences need to be aligned.
2) the part of the month, i.e. jan2 (for the file name)
3) the required month in number (if march -> 3)
4) from day (i.e 1)
5) to day ( i.e. 14)
6) optional - Insertions file


Output -
in the wd file with all the mutations in the given dates - i.e. - uk_jan2.csv
             file with all the non-synonymous mutations that are not in the novelMuteTable i.e. - unkownMutationsJan2.txt

For example -
python main.py /data3/netanel_scripts/filterSequences/uk_aligned_ref.fasta jan2 1 1 14 CombinedCsv_20210502_3.csv
