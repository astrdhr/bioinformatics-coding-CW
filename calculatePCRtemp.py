# Read DNA sequence from file
def readDNAsequence(fasta):
    """Read DNA sequence from FASTA file and return sequence
    
    Keyword arguments:
    fasta -- string containing name of an input fasta file"""
    
    fasta_content = "" # build FASTA content here
    
    # read the FASTA file
    with open(fasta, "r") as fasta_file:
        fasta_content = fasta_file.readlines()
        
    seq = "" # build sequence here

    for line in fasta_content:
        if line[0] != ">": # skip header for each sequence
            seq += line
    
    # replace all 'U' nucleotides with 'T' nucleotides, remove spaces and new lines in DNA sequence
    for pair in [["U","T"], [" ", ""], ["\n", ""]]:
        seq = seq.replace(pair[0], pair[1])
    
    # check if sequence only contains biological nucleotides
    for n in seq:
        if n not in ['A', 'C', 'T', 'G', 'U']:
            raise BadSequenceException("Not a nucleotide", n)

    return seq


# Compute sequence complement
def complement(dna_sequence):
    """Read DNA sequence and return complement
    
    Keyword arguments:
    dna_sequence -- string containing DNA sequence"""
    
    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # search through DNA sequence for complement
    new_dna = [dna_complement[x] for x in dna_sequence]
    string_new_dna = ''.join(new_dna) # build complement here
    
    # raise exception if sequence does not contain DNA nucleotides
    for n in string_new_dna:
        if n not in ['A', 'C', 'T', 'G']:
            raise BadSequenceException("Not a nucleotide", n)
            
    return string_new_dna


# Extract primers
def primer(sequence, length=20, forward=True):
    """Return the forward or reverse primer for a given DNA sequence 
    
    Keyword arguments:
    sequence -- string containing DNA sequence
    length -- integer length of primer, default is set to 20 nt
    forward -- boolean value whether to return forward or reverse primer, default set to True (i.e. forward primer)"""
    
    # assuming input sequence is sense strand
    if len(sequence) >= length:
        
        if forward == True:
            forward_primer = sequence[:length]
            return forward_primer     # forward primer is returned in 5' to 3' direction
        if forward == False:
            reverse_sequence = sequence[::-1]
            reverse_primer = complement(reverse_sequence[:length])
            return reverse_primer     # reverse primer is returned in 5' to 3' direction

    else:
        raise BadSequenceException("Sequence length is shorter than specified primer length")
        

# Calculate melting temperature
def meltingTemp(primer):
    """Compute melting temperature in degrees Celsius for given primer sequence
    
    Keyword arguments:
    primer -- string containing primer sequence"""
    
    for n in primer:
        if n not in ['A', 'C', 'T', 'G', 'U']:
            raise BadSequenceException("Not a nucleotide", n)
    
    # count each of the nucleotides in primer sequence
    count = {}
    count['A'] = primer.count('A')
    count['T'] = primer.count('T')
    count['C'] = primer.count('C')
    count['G'] = primer.count('G')
    
    # compute melting temperature
    melting_temp = 4*(count['G'] + count['C']) + 2*(count['A'] + count['T'])
    
    return melting_temp


# Calculate average melting temperature from FASTA file
def sequencePCRtemp(fasta):
    """Read FASTA file and return average melting temperature of two primers of the sequence
    
    Keyword arguments:
    fasta -- string containing name of an input fasta file"""
    
    # extract DNA sequence from FASTA file
    fasta_content = readDNAsequence(fasta)
    
    # extract primers from DNA sequence
    forward_primer = primer(fasta_content, length=20, forward=True)
    reverse_primer = primer(fasta_content, length=20, forward=False)
    
    # compute average melting temperatures of primers
    forward_tm = meltingTemp(forward_primer)
    reverse_tm = meltingTemp(reverse_primer)
    
    average_tm = (float(forward_tm) + float(reverse_tm))/2
    
    return average_tm