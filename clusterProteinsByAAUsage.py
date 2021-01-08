# import libraries
import csv

def distance(p, q):
    """Returns Euclidean distance between two given lists
    
    Keyword arguments:
    p -- list containing first set of floats
    q -- list containing second set of floats"""
    
    # raise exception if length of two lists not equal
    if len(p) != len(q):
        raise DimensionalityException("Input lists are not of equal length")
    # raise exception if list 'p' is empty
    elif not p:
        raise DimensionalityException("Empty list")
    # raise exception if list 'q' is empty
    elif not q:
        raise DimensionalityException("Empty list")
    # compute Euclidean distance between two lists
    else:
        sum_pq = 0
        for element1, element2 in zip(p,q):
            distance_pq = (element1 - element2)**2
            sum_pq += distance_pq
        
        return float((sum_pq)**(1/2))
        
        
def readTable(filename):
    """Read TSV file and return filename and corresponding fraction
    of amino acid types as a dictionary
    
    Keyword arguments:
    filename -- string containing name of TSV file to be read"""
    
    aa_type_dict = {} # build dictionary of amino acid types along with corresponding filenames here
    
    # open and read filename
    with open(filename, "r") as tsv_file:    
        next(tsv_file) # skip header
        for row in csv.reader(tsv_file, delimiter = "\t"):
            aa_type_dict[row[0]] = []
            for element in row[1:]:
                aa_type_dict[row[0]].append(float(element))
            
    return aa_type_dict


def distanceMatrix(inputfile, outputfile):
    """Read input TSV file containing fraction of amino acid
    types along with their corresponding file names, and return
    an output file containing a Euclidean distance matrix
    
    Keyword arguments:
    inputfile -- string containing name of input TSV file that
    contains filenames along with their amino acid type fractions
    outputfile -- string specifying name of outputfile to save
    distance matrix in"""
    
    input_dict = readTable(inputfile) # extract information from input file as a dictionary
    input_filenames = input_dict.keys() # extract file names from dictionary
    
    # build output file
    with open(outputfile, "w") as OUTF:
        OUTF.write("# Filename" + "\t" + "\t".join([file for file in input_filenames]) + "\n") # header
        matrix_dict = {} # build distance matrix here
        # compute Euclidean distance for each file pair  
        for k,v in input_dict.items():
            matrix_dict[k] = []
            for k1,v1 in input_dict.items():
                matrix_dict[k].append(distance(v,v1))         
        
        # build as matrix
        for k,v in matrix_dict.items():
            row_tsv = k + "\t"
            row_tsv += "\t".join(["%.2f" % elem for elem in v]) + "\n" # floats should be 2 decimal places
            OUTF.write(row_tsv)
