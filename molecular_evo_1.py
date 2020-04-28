# -----------------------------------------------------------------
# BIOL/CMPU-353
# Spring 2020
# Modeling Molecular Evolution
# Mutation Experiments
# Written by: Marc Smith and Jodi Schwarz
# Modified by: Brian Lee
#
# Description:
#   This program runs computational experiments to measure how
#   mutating a DNA sequence in different ways results in different
#   similarity measures between original and mutated sequences.
#
#   Both mutation strategies require mutating 15% of the original
#   DNA sequence, with the following differences:
#   - one strategy mutates nucleotides occurring anywhere in the
#     sequence
#   - the other strategy targets only nucleotides in the third
#     position of a codon
#
#   Which strategy will more greatly impact the similarity between
#   original and mutated protein sequences?
# -----------------------------------------------------------------


import re
import random

# Global variables

path = "/Users/beelee/PycharmProjects/molecular_evo/molecular_evo/"
aa_codon_file = path + "aa-codon-table.txt" # aa to codon conversion table
nuc_file = path + "hsp70_nuc.fasta.txt"  # DNA sequences
aa_file = path + "hsp70_aa.fasta.txt"  # protein sequences
align_file = path + "hsp70_alignments.fasta.txt"  # alignment files

codon_map = {}  # dictionary: codon -> AA
DNA_map = {}  # dictionary: animal -> DNA sequence
var_map = {}  # dictionary: pos on aligned sequence -> similarity rate

TRIALS = 100  # no. of experimental trials
MUTATION_RATE = .15  # mutation rate


#############################################################################
# main function: called at bottom of file
#############################################################################
def main():
    # Initialize DNA and protein sequences from data files
    orig_DNA_seq = read_fasta(nuc_file)
    orig_prot_seq = read_fasta(aa_file)

    create_codon_map()  # instantiate a reference for translating codons into aa
    animals = create_DNA_map()  # instantiate and save a reference to access animals' DNA
    # nuc_aligned = create_var_map(animals)  # instantiate and save a reference for the identity of each nuc pos
    create_var_map(animals)  # instantiate and save a reference for the variance of each nuc pos

    # calculate average percent identity at set mutation rate
    tot_id_1 = tot_id_2 = tot_id_3 = tot_id_4 = 0

    for trial in range(TRIALS):
        percent_id_1 = strategy1(orig_DNA_seq, orig_prot_seq)
        percent_id_2 = strategy2(1, orig_DNA_seq, orig_prot_seq)
        percent_id_3 = strategy2(2, orig_DNA_seq, orig_prot_seq)
        percent_id_4 = strategy2(3, orig_DNA_seq, orig_prot_seq)
        tot_id_1 += percent_id_1
        tot_id_2 += percent_id_2
        tot_id_3 += percent_id_3
        tot_id_4 += percent_id_4

    avg_id_1 = tot_id_1 / TRIALS
    avg_id_2 = tot_id_2 / TRIALS
    avg_id_3 = tot_id_3 / TRIALS
    avg_id_4 = tot_id_4 / TRIALS

    print("Using the human heatchuck protein as a sample gene, our model evaluated mutations at particular nucleotide "
          "positions to simulate non-synonymous mutation at the translation phase. It suggests the rate of amino acid "
          "mutation depends on the position of the mutated nucleotides in a codon:")
    print("\n---- after " + str(TRIALS) + " simulations at " + str(MUTATION_RATE * 100) + "% nucleotide mutation rate ----")
    print(" - " + str(round(avg_id_1*100)) + "% protein similarity when all nucleotides are subjected to mutation")
    print(" - " + str(round(avg_id_2*100)) + "% protein similarity when only the nucleotides in position 1 in a codon are subjected to mutation")
    print(" - " + str(round(avg_id_3*100)) + "% protein similarity when only the nucleotides in position 2 in a codon are subjected to mutation")
    print(" - " + str(round(avg_id_4*100)) + "% protein similarity when only the nucleotides in position 3 in a codon are subjected to mutation")

    # calculate the rate of having a similar nuc among the 6 animal DNA
    # sequences at each of the 3 positions on a possible codon
    percent_id = id_at_pos(id_indices(animals), only_aligned_indices(animals))
    print("\nThe 6 animal heatchuck proteins show different nucleotide matching percentages at each nucleotide position in a codon:")
    print(" - " + str(round(percent_id[0] * 100)) + "% of the time they will match at position 1\n - "
          + str(round(percent_id[1] * 100)) + "% of the time they will match at position 2\n - "
          + str(round(percent_id[2] * 100)) + "% of the time they will match at position 3")


    # calculate the variance of the nucleotide positions among the 6 animal DNA
    percent_var = var_at_pos(animals)
    print("\nIn addition, each nucleotide position in a codon has a different level of diversity, and position 3 is the most random.")
    print("If we hypothesize the chances of having any nucleotide at any position is as high as any other "
          + "(1/4 of the 6 aligned nucleotides in the 6 animal DNAs would be 'A', another forth with be 'T', and so on), "
          + "the Total Sum of Squares value represents the degree that our observed result for that position differs from the assumption. "
          + "The higher this value, the more similar those aligned nucleotides are, and 0 represents complete randomness")
    print(" - " + str(round(percent_var[0])) + " at position 1\n - "+ str(round(percent_var[1]))
          + " at position 2\n - " + str(round(percent_var[2])) + " at position 3\n - ")


#############################################################################
# create_codon_map
# initializes dictionary for aa lookups by codon
# it is wonderful python can just create a map like this!!
#############################################################################
def create_codon_map():
    in_file = open(aa_codon_file, 'r')

    for line in in_file:
        aa_list = line.split()
        aa = aa_list[0]
        for codon in aa_list[1:]:  # map all the codon variations to the aa
            codon_map[codon] = aa

    in_file.close()


#############################################################################
# create_DNA_map
# initializes dictionary to lookup the aligned DNA sequence given an animal
# returns a list of animal references
#############################################################################
def create_DNA_map():
    in_file = open(align_file, 'r')

    # hold the name and nuc sequence before save it to the dictionary
    title = ""
    content = ""
    title_ref = []

    # go line by line
    for line in in_file:
        match = re.search(r'^>(\w+)', line)
        # found a tag ine
        if match:
            content = ""  # reset sequence
            title = match.group(1)  # hold the new name
            title_ref.append(title)  # add ref to list for later
        # found a sequence line
        else:
            content += line.upper().rstrip('\n')  # reformat the content
            DNA_map[title] = content  # map ref to updated content

    in_file.close()

    return title_ref


#############################################################################
# create_var_map
# initializes a dictionary to reference the variance at the position where there is a nuc
#############################################################################
def create_var_map(ref_list):

    # get a list of indices for all DNA positions that have a nuc
    all_pos = only_aligned_indices(ref_list)

    nucs = {'A', 'T', 'C', 'G'}
    nuc_count = {}

    # the similarity at recorded pos
    for pos in all_pos:
        # 0 all counts for all nucs
        for nuc in nucs:
            nuc_count[nuc] = 0
        for ref in ref_list:
            nuc_count[DNA_map[ref][pos]] += 1
        TSS = 0
        for nuc in nucs:
            TSS += (nuc_count[nuc] - len(ref_list)/len(nucs))**2

        # # pick an animal from the list
        # for ref_1 in ref_list:
        #     # pick another animal downstream of the first animal
        #     for ref_2 in ref_list[ref_list.index(ref_1)+1:]:
        #         if DNA_map[ref_1][pos] != DNA_map[ref_2][pos]:
        #             var += 1
        # variance = var / (len(ref_list)*(len(ref_list)-1)/2)  # difference/total
        var_map[pos] = TSS


#############################################################################
# only_aligned_indices
# finds all the DNA positions that have a nuc and not a gap
# returns a list of position indices
#############################################################################
def only_aligned_indices(ref_list):

    # record the indices of which those DNA associated with the references
    # on the list have all nucleotides
    aligned_pos = []

    # length of any sequence (should all be the same)
    seq_length = len(DNA_map[ref_list[0]])  # 2282

    for pos in range(seq_length):
        # first suppose there is no gap on the position examined
        gap = False
        for ref in ref_list:
            # ignore the position if any animal has a gap in their DNA at that position
            if DNA_map[ref][pos] == '-':
                gap = True
                break
        if not gap:
            aligned_pos.append(pos)

    # print("this is the positions in total: "+str(seq_length))
    # print("this is the non-blank positions: " + str(len(aligned_pos)))

    return aligned_pos


#############################################################################
# id_ref
# finds all the DNA positions that have a matching nuc
# returns a list of position indices
#############################################################################
def id_indices(ref_list):

    # get a list of 1887 indices for all DNA positions that have a nuc
    potential_pos = only_aligned_indices(ref_list)

    id_pos = []

    # determine similarity at the recorded pos
    for pos in potential_pos:
        mismatch = False
        for ref in ref_list[1:]:
            if DNA_map[ref_list[0]][pos] != DNA_map[ref][pos]:
                mismatch = True
                break
        if not mismatch:
            id_pos.append(pos)

    # print("this is the number of matching positions: " + str(len(id_pos)))

    return id_pos


#############################################################################
# id_at_pos
# measures the frequency of the identical DNA positions among all the positions that have only nuc
# returns a list of percentages of identity
#############################################################################
def id_at_pos(id, sample):

    count_1 = total_1 = count_2 = total_2 = count_3 = total_3 = 0

    for pos in sample:
        if pos % 3 == 0:
            total_1 += 1
            if pos in id:
                count_1 += 1
        elif (pos - 1) % 3 == 0:
            total_2 += 1
            if pos in id:
                count_2 += 1
        else:
            total_3 += 1
            if pos in id:
                count_3 += 1

    # print("these are the numbers for each position: " + str(count_1) + "/" + str(total_1) + ", "
    #      + str(count_2) + "/" + str(total_2) + ", " + str(count_3) + "/" + str(total_3))
    # print("and these are the numbers for all positions: " + str(len(id)) + "/" + str(len(sample)))
    # occurance of identical nuc in the each positions of a codon
    id_pos = [count_1/total_1, count_2/total_2, count_3/total_3]

    return id_pos


#############################################################################
# id_at_pos
# calculates the rate of getting a similar nuc at a nuc site post DNA
# alignment with respect to its position on a codon
#############################################################################
def var_at_pos(ref_list):

    # get a list of indices for all DNA positions that have a nuc
    all_pos = only_aligned_indices(ref_list)

    count_1 = var_rate_1 = count_2 = var_rate_2 = count_3 = var_rate_3 = 0

    for pos in all_pos:
        if pos % 3 == 0:
            count_1 += 1
            var_rate_1 += var_map[pos]
        if (pos - 1) % 3 == 0:
            count_2 += 1
            var_rate_2 += var_map[pos]
        if (pos - 2) % 3 == 0:
            count_3 += 1
            var_rate_3 += var_map[pos]

    # average variance at each position in a codon
    avg_var = [var_rate_1/count_1, var_rate_2/count_2, var_rate_3/count_3]

    return avg_var


#############################################################################
# rand_nuc
# returns a random nucleotide from A, C, G, or T.
#############################################################################
def rand_nuc(choice):
    # offer a potentially new variant
    new_choice = random.choice(['A', 'C', 'T', 'G'])
    # prevent the new variant happen to be the same as the input
    while choice == new_choice:
        new_choice = random.choice(['A', 'C', 'T', 'G'])

    return new_choice


#############################################################################
# translate
# translates DNA strings to proteins
#############################################################################
def translate(DNA):
    # return protein sequence
    pro = ""
    # convert all nuc triplets in DNA
    for codon in re.findall('[A-Z]{3}', DNA):
        pro += codon_map[codon]
    return pro


#############################################################################
# calc_identity
# calculates the percent identity of the two given sequences.
#############################################################################
def calc_identity(seq1, seq2):
    # tally up the similar aa between 2 sequences
    # find percentage similarity
    same_codon = sum(codon1 == codon2 for codon1, codon2 in zip(seq1, seq2))
    same_percentage = same_codon/len(seq1)

    return same_percentage


#############################################################################
# read_fasta
# read fasta file into single line, return sequence
#############################################################################
def read_fasta(filename):
    in_file = open(filename, 'r')
    in_file.readline()  # ignore tag line

    sequence = ""
    # read and concatenate sequence lines without \n
    for line in in_file:
        sequence += line.rstrip('\n')

    in_file.close()
    return sequence


#############################################################################
# strategy1
# Perform a single experiment using strategy 1:
#   mutates 15% of nucleotides
#############################################################################
def strategy1(DNA_seq, prot_seq):
    # compare original DNA to a mutated new sequence
    new_DNA_seq = ""

    for nuc in DNA_seq:
        if random.random() < MUTATION_RATE:
            # randomly mutate nuc
            new_DNA_seq += rand_nuc(nuc)
        else:
            # append original nucleotide to new sequence
            new_DNA_seq += nuc

    # translate new DNA sequence to protein sequence
    new_prot_seq = translate(new_DNA_seq)

    # calculate the percent identity
    identity = calc_identity(new_prot_seq, prot_seq)

    return identity


#############################################################################
# strategy2
# Perform a single experiment using strategy 2:
#   Mutates 15% of nucleotides in third positions of codons
#############################################################################
def strategy2(pos, DNA_seq, prot_seq):
    # compare original DNA to a mutated new sequence
    new_DNA_seq = ""

    if pos == 3: pos = 0

    for count, nuc in enumerate(DNA_seq, 1):
        # mutation rate apply to all 3rd position nuc
        if count % 3 == pos:
            # the overall mutation should still be MUTATION_RATE
            # it is higher in this position to maintain the overall rate
            if random.random() < MUTATION_RATE*3:
                # randomly mutate nuc
                new_DNA_seq += rand_nuc(nuc)
            else:
                # append original nucleotide to new sequence
                new_DNA_seq += nuc
        else:
            new_DNA_seq += nuc

    # translate new DNA sequence to protein sequence
    new_prot_seq = translate(new_DNA_seq)

    # calculate the percent identity
    identity = calc_identity(new_prot_seq, prot_seq)

    return identity


# run main()
main()