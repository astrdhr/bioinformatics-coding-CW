def readAAsequence(fasta):
    """Extract and return amino acid sequence from FASTA file
    
    Keyword arguments:
    fasta -- FASTA file containing amino acid sequence"""
    
    fasta_content = "" # build contents of FASTA file here
    
    # read FASTA file
    with open(fasta, "r") as fasta_file:
        fasta_content = fasta_file.readlines()
        
    aa_sequence = "" # build amino acid sequence here
    
    # extract amino acid sequence
    for line in fasta_content:
        if line[0] != ">": # skip header for each sequence
            aa_sequence += line.strip('\n')
    
    return aa_sequence


def AAtypes(aa_seq):
    """Return fraction of polar, small and hydrophobic amino
    acids present in amino acid sequence
    
    Keyword arguments:
    aa_seq -- string containing amino acid sequence"""
    
    # dictionary of amino acid types
    aa_cat = {'polar': ['C', 'S', 'T', 'N', 'Q', 'D', 'E', 'H', 'K', 'R', 'Y', 'W'],
              'small': ['P', 'A', 'G', 'C', 'S', 'T', 'N', 'D', 'V'],
              'hydrophobic': ['F', 'Y', 'W', 'H', 'K', 'M', 'T', 'C', 'A', 'I', 'L', 'V']
             }
    
    # initialise amino acid categories
    count_polar = 0
    count_small = 0
    count_hydrophobic = 0
    
    # count number of amino acids in sequence that are present for each category
    for aa in aa_seq:
        if aa in aa_cat['polar']:
            count_polar += 1
        if aa in aa_cat['small']:
            count_small += 1
        if aa in aa_cat['hydrophobic']:
            count_hydrophobic += 1
    
    # calculate fraction of amino acids from each category present in sequence
    ratio_polar = float(count_polar) / float(len(aa_seq))
    ratio_small = float(count_small) / float(len(aa_seq))
    ratio_hydrophobic = float(count_hydrophobic) / float(len(aa_seq))
    
    return [ratio_polar, ratio_small, ratio_hydrophobic]


def AAtypetable(filelist, outputfile):
    """Extract amino acid sequence from each file and return table
    with fraction of amino acid types along with corresponding file
    name
    
    Keyword arguments:
    filelist -- string containing list of FASTA files
    outputfile -- string specifying name of outputfile to save table in"""
    
    bad_files = [] # keep track of files that do not exist, or have bad data, etc
    aa_type_dict = {}
    
    for file in filelist:

        # extract amino acid sequences from file names and compute amino acid types
        try:
            aa_seq = readAAsequence(file)
            for n in aa_seq:
                if n not in ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']:
                    raise BadSequenceException
                    pass
            aa_type = AAtypes(aa_seq)
            aa_type_rounded = ["%.2f" % elem for elem in aa_type]
            aa_type_dict[file] = aa_type_rounded
        
        # continue building table if input file not found or file contains wrong data
        except (FileNotFoundError, BadSequenceException):
            bad_files.append(file)
            continue
                  
    # output to specified output file in TSV format
    with open(outputfile, "w") as OUTF:
        OUTF.write("# Filename" + "\t" + "Polar" + "\t" + "Small" + "\t" + "Hydro" + "\n")
        for k,v in aa_type_dict.items():
            row_tsv = k + "\t"
            row_tsv += "\t".join([str(aatype) for aatype in v]) + "\n"
            OUTF.write(row_tsv)
    
    # return list of only bad files
    if bad_files:
        return bad_files
    else:
        return []
