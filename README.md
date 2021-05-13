# covidannotator

This program find all the mutations in given sequences.

First activate conda - > conda activate nextstrain (or other env with pandas)

Input - 
in the wd folder have to be -"regions.csv" - also attached to this git.
as input in the CLI -> the Fasta file, the first sequence must be the reference, other sequences need to be aligned.


Output -
in the working directory "all_mutations" csv file with all the mutations in the sequences

Running example -

conda activate nextstrain

python main.py /data3/netanel_scripts/filterSequences/uk_aligned_ref.fasta
