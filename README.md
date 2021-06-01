# covidannotator

This program find all the mutations in given sequences.

First make sure you are on sudo mode, and the relevent conda enviorment is on - 

1) sudo -s

2) conda activate nextstrain (or other env with pandas)

3) make sure that "regions.csv" and "b117muts.csv" is in the WD (download from this github)


Input - 
the Multi-Fasta file, the first sequence must be the reference, other sequences need to be aligned.

Optional Input:
-i (PATH for insertions file) if you want to include insertions to the excel table
-n to include 'N' positions to the excel table

Output -
in the working directory "all_mutations" csv file with all the mutations in the sequences

Running example -

python main.py /data3/netanel_scripts/filterSequences/uk_aligned_ref.fasta
