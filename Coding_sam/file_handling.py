### Author: Laura Gutierrez Funderburk
### Supervisor: Dr. Cedric Chauve
### Created on: Nov 18 2017
### Last edited: January 10 2018

# for handling all related to gene Sequence 
from Bio import SeqIO

# Numpy arrays, and performing operations on arrays
import numpy as np

# Defining panda dictionary for file access
import pandas as pd

# For getting all family names in synteny files without repetition while also keeping order of appearance
import ordered_set

import sys

# Input and Output locations

# sys.argv[0]
# sys.argv[1]

input_location = '/home/lgutierrezfunderburk/Documents/Test/FASTA/'
output_location = '/home/lgutierrezfunderburk/Documents/Test/Save_results/2018/synteny_in_x_and_n-x/synteny_in_3_and_n-3/'
#output_location = '/home/lgutierrezfunderburk/Documents/Test/Save_results/synteny_in_2_and_n-2_species/F2/'

########################################################################################################################################################################################################################################################################
"""Tools"""

#Define dictionary
file_dictionary = {'maculatus':'Anopheles-maculatus-Maculatus3_SCAFFOLDS_AmacM1.fa',
                  'epiroticus':'Anopheles-epiroticus-Epiroticus2_SCAFFOLDS.AepiE1.fa',
                  'atroparvus':'Anopheles-atroparvus-EBRO_SCAFFOLDS_AatrE1.fa',
                  'sinensis':'Anopheles-sinensis-SINENSIS_SCAFFOLDS_AsinS1.fa',
                  'melas':'Anopheles-melas-CM1001059_SCAFFOLDS_AmelC1.fa',
                  'merus':'Anopheles-merus-MAF1_SCAFFOLDS_AmerM1.fa',
                  'stephensi':('Anopheles-stephensi-SDA-500_SCAFFOLDS_AsteS1.fa',\
                               'Anopheles-stephensiI-Indian_SCAFFOLDS_AsteI2.fa'),
                  'darlingi':('Anopheles-darlingi-Coari_SCAFFOLDS_AdarC3.fa',\
                              'Anopheles-darlingi-Coari_SCAFFOLDS_AdarC2.fa'),
                  'gambiae':('Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa',\
                            'Anopheles-gambiae-PEST_SCAFFOLDS_AgamP3.fa'),
                  'minimus':'Anopheles-minimus-MINIMUS1_SCAFFOLDS_AminM1.fa',
                  'arabiensis':'Anopheles-arabiensis-Dongola_SCAFFOLDS_AaraD1.fa',
                  'farauti':'Anopheles-farauti-FAR1_SCAFFOLDS_AfarF1.fa',
                  'quadriannulatus':'Anopheles-quadriannulatus-SANGWE_SCAFFOLDS_AquaS1.fa',
                  'funestus':'Anopheles-funestus-FUMOZ_SCAFFOLDS_AfunF1.fa',
                  'dirus':'Anopheles-dirus-WRAIR2_SCAFFOLDS_AdirW1.fa',
                  'christyi':'Anopheles-christyi-ACHKN1017_SCAFFOLDS_AchrA1.fa',
                  'culicifacies':'Anopheles-culicifacies-A37_SCAFFOLDS_AculA1.fa',
                  'albimanus':'Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa'}
pd_dictionary = pd.DataFrame(file_dictionary)

# Key words
key_words= {'maculatus', 'epiroticus', 'atroparvus', 'sinensis', 'melas', 'merus', 'stephensi', 
'darlingi', 'gambiae', 'minimus', 'arabiensis', 'farauti', 'quadriannulatus', 'funestus', 'dirus', 
'christyi', 'culicifacies', 'albimanus'}

########################################################################################################################################################################################################################################################################

########################################################################################################################################################################################################################################################################

""" Function Definition Area, Part I: Automation and File Handling""" 

ALL_GENE_FILE_DIRECTORY = "./DATA/ALL_GENE_file"
SYNTENY_RESULTS_DIRECTORY= "./synteny_results/"

### STEP 0. Give the script a way to access adequate FASTA files
            
# Given a string of information, such as 
#  Anopheles_albimanus	KB672457	MZ22528743	AALB006655	-	4678482	4679472	3	4678482-4678753:4678832-4679255:4679290-4679472
# We want to be able to use this information once we start working with Species name, 
#       Fragment, Gene name, Gene orientation, and Exons Coordinates
# We define a function whose input is a string of the form:
# "Anopheles_albimanus	KB672457	MZ22528743	AALB006655	-	4678482	4679472	3	4678482-4678753:4678832-4679255:4679290-4679472"
# separated by tabs, and whose output is an array. In our example, the function will return:
# ['Anopheles_albimanus','KB672457,'MZ22528743','AALB006655,'-','4678482','4679472','3','4678482-4678753:4678832-4679255:4679290-4679472']
def get_information(string):
    
    """This function takes as input a string of characters such that each word is separated by tabs/space 
    and returns an array such that each entry corresponds to a word"""
    information = string.split("\t")
    return information    

# This function takes as input
#     a string of the form
#     "Species_name  Fragment    Family_name Gene_name   Orientation Start_exon  End_exon    Number_exons    Exon_coordinates"
# For example
#     "Anopheles_albimanus	KB672457	MZ22528743	AALB006655	-	4678482	4679472	3	4678482-4678753:4678832-4679255:4679290-4679472"
# And outputs a string with the file name of the appropriate FASTA file 
# In our example, the FASTA file that will be opened is called 
#     "Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa"

# This function uses our tool file_dictionary whose keys are species names, and values correspond to files. 
def get_FASTA_file(string):
    species_info = get_information(string)
    key_to_file = species_info[0].split("_")[1]
    FASTA_files = pd_dictionary[key_to_file]
    fragment = species_info[1]
    array = []
    k = 0
    while not array: 
        array = [input_location + FASTA_files[k] for key in get_sequence_from_file(input_location + FASTA_files[k]) if key==fragment]
        k +=1
    return (array[0])

### STEP ONE: WHERE IN THE ALL_GENE_FILE FILE IS THIS FAMILY NAME?

def get_families_from_synteny_file(synteny_file_name):
    num_lines = sum(1 for line in open(synteny_file_name))
    with open(synteny_file_name) as f:
        head = [next(f) for x in range(num_lines)]
    f.close()
    return head    

def row_families_from_synteny(synteny_file,row_num):
    families = [get_families_from_synteny_file(synteny_file)[row_num].split("\t")[i].split(":")[0] for i in range(4)]
    return families

def row_genes_from_synteny(synteny_file,row_num):
    genes = [get_families_from_synteny_file(synteny_file)[row_num].split("\t")[i].split(":")[1] for i in range(4)]
    return genes

    
def find_family_in_ALL_GENE_FILE(family_name):
    with open(ALL_GENE_FILE_DIRECTORY,'r') as inF:
        found_fam_in = [line for line in inF if family_name in line]
    return found_fam_in

###This function's purpose is to locate where in the ALL_GENE_FILE file 
#   a given family name is found, and upon obtaining that information, 
#   storing the entries where the family name was found.

def step_one(synteny_file_name,row_num,col_number):
    synteny_file_in_table = get_families_from_synteny_file(synteny_file_name)
    families = row_families_from_synteny(synteny_file_name,row_num)
    target_gene_family = families[col_number]
    return target_gene_family
    
def output_entries_in_ALL_GENE_FILE(target_family):
    get_family = find_family_in_ALL_GENE_FILE(target_family)
    return get_family



# STEP TWO
# GREAT! CAN I PLEASE GET THE FASTA HANDLE FOR ALL THOSE ENTRIES?

def get_header_FASTA(string):
    info = get_information(string)
    FASTA_file = get_FASTA_file(string)
    fragment = info[1]
    with open(FASTA_file,'r') as f:
        save_line = [line for line in f if fragment in line]
    f.close()
    contig = [save_line[i] for i in range(len(save_line)) if save_line[i].split(">")[1].split(" ")[0] == fragment]
    header = contig[0]
    return header
    #with open(output_location + "HEADERS/"+ info[2] + str("_") +info[1] + "_header",'w') as outfile:
    #    for letter in header:
    #        outfile.write(letter)
    #outfile.close()
########################################################################################################################################################################################################################################################################


""" Function Definition Area, Part II: Sequences"""


# Define complement function. Input is a string containing gene sequence, and output is a string containing
# the complement of the gene sequence

def complement(seq):
    
    """This function takes a sequence, and computes its complement"""
    
    # Define a dictionary with sequence symbols as keys, and their corresponding complement symbol as dictionary value
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n'} 
    
    # We define bases as a list where each element of the list corresponds to one letter in the seq string
    bases = list(seq) 
    
    # For each element in the list bases, extract the corresponding complement value from dictionary, and construct
    # a new list with the complement elements in orderly manner
    complement_bases = [complement[base] for base in bases] 
    
    # use join() method to join elements in complement_bases array, and return a string with the complement sequence
    return complement_bases


# Function from_array_to_string() takes as input an array whose elements are strings containing one letter
# and whose output is all elements of the array contatenated into a string
def from_array_to_string(sequence_array):
    """ This function will take all elements in the sequence array, and concatenate them into a single string"""
    return ''.join(sequence_array)

# Define reverse function. Input is a string containing a sequence, and output is a string with the sequence
# in reverse order

def reverse(comp_seq):
    """This function takes a sequence, and computes the sequence in reverse order"""
    
    # take input, and return the sequence in reverse order
    return comp_seq[::-1]

# Define get_sequence_from_file() function. 
# This function takes as input a FASTA file, and outputs an array whose elements are the letters in the sequence

def get_sequence_from_file(file_name):
    """This functions will read FASTA file once, and store the content of the sequence in an array"""
    # file_name is the file containing FASTA file
    
    # we parse the FASTA file, and obtain the sequence in it, storing it in array 
    sequences = {seq_record.name: seq_record.seq for seq_record in SeqIO.parse(file_name,'fasta')}
    
    # once we obtained the sequence, we want to store it in an array such that each element in the dictionary
    # is one letter from the sequence
    return(sequences)

# Define exon_coordinate_partition() a function whose input is an array containing elements of a sequence, and two
#    coordinates: start and end exons
#    The function will return an array containing elements in the sequence corresponding to the outlined coordinates
#
#  For example, if given a sequence=['a','a','g','g','c','t','a'] and the coordinates 3,6, this function
#   will return [g','g','c','t']
def exon_coordinate_partition(sequence,fragment, coordinate_array):
    
    """This function will partition our sequence as per given coordinates"""
    
    if (len(coordinate_array) == 1):
        parted = np.array(sequence[fragment][coordinate_array[0][0]-1:coordinate_array[0][1]])
    else:
        par = np.array([sequence[fragment][item[0]-1:item[1]] for item in coordinate_array])
        parted = np.sum(par)
    return parted



# We next define a function that allows us to use the exon coordinates. From our example, we are given coordinates
# 4678482-4678753:4678832-4679255:4679290-4679472
# And we want to get all three different exon coordinates. We define a function that takes as input a string containing
#    the exon coordinates, and which outputs an array with pairs of coordinates (X,Y) as elements. In our example
# we get [(4678482,4678753), (4678832,4679255),(4679290,4679472)]

# We define extract_coordinates() as a function which takes as input
#     string of the form exon_start_1 - exon_end_1 : exon_start_2 - exon_end_2 : ... :exon_start_n - exon_end_n
# And whose output is an array of the form
#     [(exon_start_1 , exon_end_1), (exon_start_2, exon_end_2) , ... ,(exon_start_n , exon_end_n)]
def extract_coordinates(string_array):
    
    """ This function takes a string with exon coordinates, and outputs an array whose entries are pairs of the form
    (X,Y) where X,Y are integers, such that X corresponds to exon start and Y to exon end"""
    
    
    # apply split() method on input string
    split_stage_one = string_array.split(":")
    
    # apply split() method on split_stage_one and store results in array
    split_stage_two = [item.split("-") for item in split_stage_one]
    
    # For each entry in split_stage_two, get pairs (X,Y) and transform values of X,Y into integer values
    int_split_final_stage = [ (int(split_stage_two[i][0]),int(split_stage_two[i][1])) \
                             for i in range(len(split_stage_two))]
    
    # return array [(X_1,Y_1),(X_2,Y_2),...,(X_n,Y_n)]
    return int_split_final_stage

# We next want to decide whether to compute the reverse complement of a given sequence, depending on the orientation
# Orientation is a string '+' or a string '-'
# Orientation can be found on the fifth position in the array we obtained previously, in our example
# Anopheles_albimanus	KB672457	MZ22528743	AALB006655	-	4678482	4679472	3	4678482-4678753:4678832-4679255:4679290-4679472
# We see the orientation is '-'. 
# If '-' we compute reverse complement, if '+' we leave the sequence as it is

# We define a function decide_for_reverse_complement() whose input is a 
#      string '-' or '+', and a sequence
# And whose output is either, respectively, the reverse complement of the sequence or the sequence itself
def decide_for_reverse_complement(string_sign,sequence):
    """This function will decide whether or not to compute the reverse complement of a given sequence"""
    
    # If orientation is negative
    if string_sign == '-':
        # Get complement of the sequence
        comp = complement(sequence)
        # Get reverse of the complement
        reverse_complement = reverse(comp)
        # Return reverse complement
        return from_array_to_string(reverse_complement)
    # If orientation is positive
    elif string_sign == '+':
        # Return original sequence
        return sequence
    # If orientation is neither '-' nor '+', return an error message
    #else:
    #    error = "Error, invalid value"
    #    return error
    
########################################################################################################################################################################################################################################################################


""" Function Definition Area, Part III: Storage"""

    
# This function takes as input a string of the form 
#  "Species_name  Fragment    Family_name Gene_name   Orientation Start_exon  End_exon    Number_exons    Exon_coordinates"
# For example
#   "Anopheles_albimanus    KB672457	MZ22528743	AALB006655	-	4678482	4679472	3	4678482-4678753:4678832-4679255:4679290-4679472"
# Gets appropriate FASTA file, extracts subsequence, decides whether to compute the reverse complement, and saves on file
#  for later analysis using MUSCLE or MAFFT
def store_sequences(string):
    
    # Get adequate FASTA file, as per our file_dictionary
    fasta_file = get_FASTA_file(string)
    
    # Prep input info for use
    info_array = get_information(string)
    
    # Get coordinates from adequate position
    coordinates = extract_coordinates(info_array[8])
    
    # Extract sequence from file
    Dict_seq = get_sequence_from_file(fasta_file)
    
    # Use exon coordinates and parition sequence
    concatenate_exons = exon_coordinate_partition(Dict_seq,info_array[1],coordinates)
    print(concatenate_exons)
    # Compute reverse complement if applicable
    final_sequence_arr = decide_for_reverse_complement(info_array[4],concatenate_exons)
    
    final_sequence  = from_array_to_string(final_sequence_arr)
    print(" ")
    print(final_sequence)
    
    # STORE IN ARRAY
    # Get Header
    head = get_header_FASTA(string)
    
    # Save all in array
    output_array = [head, final_sequence,"\n"]
    
    # Write on file
    #with open(output_location + info_array[2] + str("_") +info_array[1] + str("_")  + \
    #          info_array[3]+ str("_") + info_array[0] ,'w') as outfile:
    #    for item in write_array:
    #        outfile.write(item)
    #outfile.close()
    return output_array   

### This function takes as input an array like the one generated using the function above, and as output it yields a string
###  containing the concatenated elements found in the array

### Sample Suppose we feed the function store_sequences(string) the following three strings associated 
#     to the family MZ22514725 

#    Anopheles_gambiae	2L	MZ22514725	AGAP007354	+	46008993	46011382	2	46008993-46009709:46011119-46011382
#    Anopheles_melas	KI427838	MZ22514725	AMEC005293	-	3804	6159	2	3804-4067:5443-6159
#    Anopheles_merus	KI438994	MZ22514725	AMEM017160	+	192755	195221	2	192755-193522:194958-195221

# We can then use store_sequences(string) and store all genes into a single array of the form

# output_array = [['>2L chromosome:AgamP3:2L:1:49364325:1\n', 'ATGACGGTCACGTACGGTGTCGAGAAGTGCCCGGAAAAGTTCGACGAGTATCGGTGCAGCCTGCCCAGCTGGCGGCGGACCAGCTGCGCGAGGACGATCACATCCGCAAGCAGGCCATCGA', '\n'], ['>KI427838 dna:supercontig supercontig:GCA_000473525.1:KI427838:1:52148:1\n', 'ATGACGGTCACGTACGGTGTCGAGAAGTGTCCGGAAAAGTTCGACGAGTATCGGTGCAGCCTGCCCGAGCAGTACCGGAAGCTGGCCTGCGCGAGGACGATCACATCCGCAAGCAGGCCATCGA', '\n'], ['>KI438994 dna:supercontig supercontig:GCA_000473845.1:KI438994:1:1088242:1\n', 'TTGGAGGTGATCGTGAACACAACAGGAACAACACAACACAGCAGCAGCATCATGACGGTCACGTACGGTGTCGAGAAGTGTCCGGAAAAGTTCGACGAGTATCGGTGCAACCTGCCCGAAAGCT', '\n']]

# We can then turn this array into a single string that will later be written on a file for later sequence alignment

# We can then  use the function below to turn our array of arrays into a single string
def turn_to_string(output_array):
    """This functions takes as input an array containing FASTA headers and gene sequences and 
    outputs a string formed by the array elements"""
    
    # Empty string
    text = ""
    # Length of Array
    leng = len(output_array)
    
    # Perform string concatenation
    for i in range(leng):
        for item in output_array[i]:
            text += item
    return text

# This function will take as input a dictionary and a string indicating whether
#   we are working with F1, F2 or F1_F2. It will then store the gene sequences associated to MZ22514725 and store them into a file

# A word on this input dictionary:  keys are family names (i.e. MZ22514725) and values are strings associated to
#  that family. Strings associated are those that appear in the ALL_GENE_FILE file upon searching for that particular family name.
def store_seq_on_file(Fi_family_gene_dic,fam_name_str):
    
    """This function will store gene sequences for each family found in synteny files"""
    for key in Fi_family_gene_dic:
        # Part I GENE STRING (strings as they appear upon searching for family in ALL_GENE_FILE file)
        gene_string = Fi_family_gene_dic[key]
        
        # Some of the input strings have additional information which should not be there. This additional information is normally
        #   found in the coordinates. We will perform a check to see if such extra info is present, and if so, we will skip that family
        test = []
        for item in gene_string:
            info = get_information(item)
            value = " " in info[8]
            test.append(value)
            #print(test)
        if True in test:
            continue
            
        # Get family name for adequate file labeling 
        fam_name = info[2]
    
        # Part II GENE SEQUENCE (genes associated to Family, obtained from FASTA files)
        gene_seq = [store_sequences(item) for item in gene_string]
    
        # We now want to turn into a single string all those gene sequences we found for each family
        sequence = turn_to_string(gene_seq)
    
        # Time to store sequence in FASTA file
        with open(output_location +  fam_name_str + fam_name + ".fa",'w') as outfile:
            for item in sequence:
                outfile.write(item)
        outfile.close()
    
########################################################################################################################################################################################################################################################################    

########################################################################################################################################################################################################################################################################    



