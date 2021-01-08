# bioinformatics-coding-CW
A collection of Python functions I wrote to automate common molecular biology analytical processes, as part of the Coding module in my Bioinformatics MSc. I put this up as an indication of my Python ability. Summary of code:
- `calculatePCRtemp`: reads a given FASTA file as input and returns the average melting temperature of the two primers of the sequence.
- `findLongestORF`: reads a DNA sequence as input and returns the candidate protein sequence corresponding to the longest open reading frame (ORF) in FASTA format.
- `classifyAA`: extracts amino acid sequence from multiple files as input and returns statistics on amino acid type.
- `clusterProteinsByAAUsage`: reads TSV file containing amino acid type statistics as input and returns Euclidean distance matrix.
