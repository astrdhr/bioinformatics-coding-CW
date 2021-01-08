# import libraries
import re

def translate(dna_seq):
    """Read DNA sequence and return a dictionary containing the translation
    of the sequence in all possible reading frames
    
    Keyword arguments:
    dna_seq -- string containing DNA sequence to be translated"""
    
    # dictionary of codons and their corresponding amino acids
    codon_dict = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
    
    dna_dict = {} # build all DNA reading frames here
    
    # get forward reading frames
    dna_dict["f1"] = dna_seq
    dna_dict["f2"] = dna_seq[1:]
    dna_dict["f3"] = dna_seq[2:]
    
    # get reverse reading frames
    rev_dna_seq = dna_seq[::-1]
    dna_dict["r1"] = rev_dna_seq
    dna_dict["r2"] = rev_dna_seq[1:]
    dna_dict["r3"] = rev_dna_seq[2:]
    
    # translate DNA reading frames to amino acid sequences   
    translated_dict = {}  
    for frame,seq in dna_dict.items():
        pos = 0
        translated_seq = ""
        while pos+3 <= len(seq):
            codon = seq[pos:pos+3]
            aa = codon_dict[codon]
            translated_seq += aa
            pos+=3
        translated_dict[frame] = translated_seq
            
    return translated_dict


def openReadingFrame(aa_seq):
    """Read amino acid sequence and return the sequence of the
    first open reading frame
    
    Keyword arguments:
    aa_seq -- string containing the amino acid sequence"""
    
    # read amino acid sequence
    for aa in aa_seq:
        if "M" not in aa_seq:
            aa_content = ''
        elif "*" not in aa_seq:
            aa_content = ''
        else:
            aa_content = re.search("M(.*?)(?=\*)", aa_seq).group(0) # search for ORF
            
    return aa_content


def candidateProtein(dna_seq):
    """Read DNA sequence and return translated amino acid sequence
    corresponding to the longest open reading frame
    
    Keyword arguments:
    dna_seq -- string containing DNA sequence"""
    
    # convert DNA sequence to amino acid sequence in all possible reading frames
    translated_frames = translate(dna_seq)
    
    # extract ORFs from each reading frame
    ORFs = [] # build ORFs here
    for frame,aa_seq in translated_frames.items():
        an_ORF = openReadingFrame(aa_seq)
        ORFs.append(an_ORF)
    
    # extract longest ORF
    longest_ORF = max(ORFs)
    
    return longest_ORF


def writeFASTA(sequence, description, filename):
    """
    Create a FASTA file given sequence, description and filename as input
    
    Keyword arguments:
    sequence    -- string containing amino acid sequence 
    description -- string containing description (e.g. name of protein, organism, etc)
    filename    -- string containing file name
    """
    
    # create file
    OUTF = open(filename, "w")    # assuming fasta extension is included in given filename
    
    # write description to file as FASTA header
    header = ">" + description + "\n"
    OUTF.write(header)
    
    # remove of tabs, newlines and spaces in sequence string
    for pair in [[" ", ""], ["\t", ""], ["\n", ""]]:
        sequence = sequence.replace(pair[0], pair[1])
    
    # output in 60 char lines for convenience
    line_width = 60
    pos = 0
    while sequence[pos:pos+line_width] != '':
        OUTF.write(sequence[pos:pos+line_width] + "\n")
        pos += line_width

    # close file
    OUTF.close()
    
    
def maximalORF(inputfile, outputfile, proteinname):
    """
    Read DNA sequence and write the candidate protein corresponding
    to the longest ORF in FASTA format
    
    Keyword arguments:
    inputfile -- string containing name of an input file
    outputfile -- string containing name of an output file
    proteinname -- string with a description of a candidate protein
    """
    
    # read and extract DNA sequence from specified input file
    file = open(inputfile, "r")
    dna_seq = "" # build DNA sequence here
    for line in file:
        if ">" not in line:
            dna_seq += line.rstrip()
    file.close()
    
    # convert DNA seq to amino acid seq and extract maximum ORF
    orf_seq = candidateProtein(dna_seq)
    
    # write to specified output file in FASTA format
    writeFASTA(orf_seq, proteinname, outputfile)