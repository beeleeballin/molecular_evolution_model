#--------------------------------------------------------------------
# BIOL/CMPU-353
# Spring 2020
# Modeling Molecular Evolution
# Mutation Statistics (Starter code)
#
# Program: mutation_stats.py 
#
# Written by: Marc Smith and Jodi Schwarz
# Modified by: <your name here>
#
# Description:  This program
#
#   1. reads in six aligned DNA sequences of HSP70 from
#      Celegans, Rat, Drosophila, Mouse, Human, and Nematostella.
#      (we've given you a head start by completing this for you!)
#
#   2. calculates diversity of nucleotides in each codon position
#      (i.e., how many different nucleotides occur in each position)
#
#   3. calculates average diversity for codon positions 1, 2, and 3
#      (i.e., the average number of different nucleotides occurring
#       in position 1 of each codon, position 2 of each codon, and
#       position 3 of each codon)
# 
# Note: dealing with gaps ("-") in aligned sequences is not trivial.
#       you may reasonably ignore gaps when calculating the diversity
#       of nucleotides in each position of the aligned sequences. 
#
#--------------------------------------------------------------------

import re

# file containing six aligned nucleotide seqs of HSP70 
aligned_nuc_seqs = "multifasta_nuc_aligned.txt"

#############################################################################
# main function: called at bottom of file
#############################################################################
def main():

    # Open and ignore first line from aligned multi-fasta file
    in_file = open(aligned_nuc_seqs, 'r')
    in_file.readline()

    # Initialze the six HSP70 sequences from multi-fasta file
    seq1 = read_seq_from_file( in_file )
    seq2 = read_seq_from_file( in_file )
    seq3 = read_seq_from_file( in_file )
    seq4 = read_seq_from_file( in_file )
    seq5 = read_seq_from_file( in_file )
    seq6 = read_seq_from_file( in_file )
	
    # Close fasta file
    in_file.close()

    # 2. Calculate the diversity of nucs in each column of  
    #    a 6-sequence alignment 
    #    (uncomment the next line after you implement calc_col_div() 
    #div_counts = calc_col_div( seq1, seq2, seq3, seq4, seq5, seq6 )

    # 3. Calculate and print diversity statistics
    #    (uncomment the next line after you implement calc_codon_stats()
    #calc_codon_stats( div_counts )


#############################################################################
# calc_pos_div 
# 
# Calculates the nucleotide diversity in each column of a 6-sequence 
# alignment. The diversity count of a column is the number of different
# nucs that appear in that column.
# 
# INPUTS: seq1 - seq6
# OUTPUTS: list of diversity counts
#############################################################################
def calc_col_div( seq1, seq2, seq3, seq4, seq5, seq6 ):

    # initialize an empty list of diversity counts
    div_counts = []

    # for each nuc from the 6 aligned sequences
    # (n1 is from seq1, n2 is from seq2, etc.)
    for n1, n2, n3, n4, n5, n6 in zip(seq1, seq2, seq3, seq4, seq5, seq6):
	
        # make a string out of each nuc in current column of alignment
        # --col contains all the nucs (and gaps) in current column 
        # col = ...

        # count number of A's, C's, T's, and G's in this col
        #num_a = ...
        #num_c = ...
        #num_t = ...
        #num_g = ...

        # count total number of different nucleotides in this column. 
        # e.g.,
        # - if col contains all A's, unique_nucs is 1, 
        # - if col contains half A's and half C's, unique_nucs is 2
        # - etc.
        #unique_nucs = ...
        
        # append count of unique nucs in this column 
        # to the list of diversity counts
        #div_counts.append(unique_nucs)

    # return the list of diversity counts (# of unique nucs in each col)
    return div_counts


#############################################################################
# calc_codon_stats
# 
# Calculates and prints the average diversity by nucleotide position 
# within each codon
#
# INPUTS: list of diversity counts for each position in alignment
#############################################################################
def calc_codon_stats( div_counts ):

    # initialize sums for each position in codons
    sum_pos_0 = 0.0   # sum of diversity in codons' first position
    sum_pos_1 = 0.0   # sum of diversity in codons' second position
    sum_pos_2 = 0.0   # sum of diversity in codons' third position


    # loop through the list of diversity counts, and add each of the 
    # counts to the appropriate sum (sum_pos_0, sum_pos_1, sum_pos_2)


    # now that you have the sums for each codon position, 
    # calculate average diversity for each codon position
    # note: the average for each position is the respective sum 
    #       divided by the number of codons 
    #num_codons = ...
    #avg_div_0 = ...
    #avg_div_1 = ...
    #avg_div_2 = ...
    

    # print stats
    # (uncomment once you've computed the averages)
    #print("\nAverage diversity by nucleotide position within codons:")
    #print("\tCodon position 1: " + str(avg_div_0))
    #print("\tCodon position 2: " + str(avg_div_1))
    #print("\tCodon position 3: " + str(avg_div_2))
    #print("\n")


#############################################################################
# read_seq_from_file 
# Returns next DNA sequence from File
#############################################################################
def read_seq_from_file( file ):
    seq = ""
    line = file.readline().rstrip('\n')
    while (line != '' and not re.match('^>', line)):
        seq = seq + line
        line = file.readline().rstrip('\n')
    return seq

main()

