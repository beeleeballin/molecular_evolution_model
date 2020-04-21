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
alignment_map = {}  # dictionary: animal -> DNA sequence
sim_map = {}  # dictionary: pos on aligned sequence -> similarity rate

TRIALS = 2  # no. of experimental trials
MUTATION_RATE = .15  # mutation rate


#############################################################################
# main function: called at bottom of file
#############################################################################
def main():
    # Initialize DNA and protein sequences from data files
    orig_DNA_seq = read_fasta(nuc_file)
    orig_prot_seq = read_fasta(aa_file)

    create_codon_map()  # instantiate a reference for translating codons into aa
    animals = create_alignment_map()  # instantiate a reference for animals' DNA sequences
    sim_nuc = create_sim_map(animals)  # instantiate a reference for the mutation rate of each nuc pos


    # # Conduct TRIALS no. of experiments using strategy 1 & 2
    # print("Results for random mutation strategy\n")
    # print("Trial, % identity\n")
    # for trial in range(TRIALS):
    #     percent_id = strategy1(orig_DNA_seq, orig_prot_seq)
    #     print("Trial "+str(trial + 1) + ", " + str(percent_id))
    # print("\n")
    # print("Results for 3rd codon position mutation strategy\n")
    # print("Trial, % identity\n")
    # for trial in range(TRIALS):
    #     percent_id = strategy2(orig_DNA_seq, orig_prot_seq)
    #     print(str(trial + 1) + ", " + str(percent_id))
    # print("\n")

    # calculate average percent identity at set mutation rate
    tot_id_1 = 0
    tot_id_2 = 0

    for trial in range(TRIALS):
        percent_id_1 = strategy1(orig_DNA_seq, orig_prot_seq)
        percent_id_2 = strategy2(orig_DNA_seq, orig_prot_seq)
        tot_id_1 += percent_id_1
        tot_id_2 += percent_id_2

    avg_id_1 = tot_id_1/TRIALS
    avg_id_2 = tot_id_2/TRIALS

    print("\nOur model suggests that the prevalence of non-synonymous mutation depends on the nucleotide position in a codon:")
    print("---- after " + str(TRIALS) + " trials at a " + str(MUTATION_RATE*100) + "% rate of mutation ----\n")
    print(" - a " + str(round(avg_id_1*100)) + "% identity similarity when all positions are subjected to mutation")
    print(" - a " + str(round(avg_id_2*100)) + "% identity similarity when only the 3rd positions are subjected to mutation")

    # calculate the rate of having a similar nuc among the 6 animal DNA
    # sequences at each of the 3 positions on a possible codon
    print("\nBased on the animal HSP40 proteins, our model suggests that the rate of mutation at each position in a codon varies:\n")
    percent_id = sim_at_pos(sim_nuc)
    print(" - the identity similarity in position 1 is " + str(round(percent_id[0]*100))
          + "%\n - the identity similarity in position 2 is " + str(round(percent_id[1]*100))
          + "%\n - the identity similarity in position 3 is " + str(round(percent_id[2]*100))
          + "%")


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
# create_alignment_map
# initializes dictionary to lookup the aligned DNA sequence given an animal
#############################################################################
def create_alignment_map():
    in_file = open(align_file, 'r')

    # hold the name and nuc sequence before save it to the dictionary
    name = ""
    sequence = ""
    ref = []

    for line in in_file:
        match = re.search(r'^>(\w+)', line)
        # found tag ine
        if match:
            sequence = ""
            name = match.group(1)  # save the new animal name
            ref.append(name)  # add it to animals[] for reference
        # found sequence line
        else:
            sequence += line.upper().rstrip('\n')  # reformat DNA
            alignment_map[name] = sequence  # map name to updated DNA

    in_file.close()
    # for animal in animals:
    #     print(animal)
    #     print(alignment_map[animal])

    return ref


#############################################################################
# create_mutation_rate_map
# initializes a dictionary to reference the rate of getting a similar nuc
# given the nuc position
#############################################################################
def create_sim_map(align_ref):

    sim_ref = []
    # if there exists a nuc on a certain pos across DNAs of all animals
    aligned_seq_length = len(alignment_map[align_ref[0]])  # 2282

    for pos in range(aligned_seq_length):
        gap = False
        for ref in align_ref:
            if alignment_map[ref][pos] == '-':
                # print(animal+" "+str(pos))
                gap = True;
                break;
        if not gap:
            sim_ref.append(pos)

    # print(not_gap)  # the indices to be measured
    # print(str(len(not_gap))+" pos out of "+str(aligned_seq_length)+" will be measured")
    # print(animals)  # all animals

    # the similarity at recorded pos
    for pos in sim_ref:
        sim = 0
        for ref_1 in align_ref:
            for ref_2 in align_ref[align_ref.index(ref_1)+1:]:
                #print(ani_1+" compared to "+ani_2+" at pos "+str(pos))
                if alignment_map[ref_1][pos] == alignment_map[ref_2][pos]:
                    sim += 1
        # print("there are "+str(sim)+" out of "+str(len(animals)*5/2)+" nucleotides that are same at "+str(pos))
        sim_rate = sim / (len(align_ref)*5/2)  # similar/total
        sim_map[pos] = sim_rate

    return sim_ref


#############################################################################
# rand_nuc
# returns a random nucleotide from A, C, G, or T.
#############################################################################
def rand_nuc(choice):
    new_choice = random.choice(['A', 'C', 'T', 'G'])
    while choice == new_choice:
        new_choice = random.choice(['A', 'C', 'T', 'G'])
    return new_choice


#############################################################################
# translate
# translates DNA strings to proteins
#############################################################################
def translate(DNA):
    pro = ""
    # all nuc triplets in DNA
    for codon in re.findall('[A-Z]{3}', DNA):
        pro += codon_map[codon]
    return pro


#############################################################################
# calc_identity
# calculates the percent identity of the two given sequences.
#############################################################################
def calc_identity(seq1, seq2):
    # tally up the similar aa between 2 sequences, find percentage similarity
    same_codon = sum(codon1 == codon2 for codon1, codon2 in zip(seq1, seq2))
    same_percentage = same_codon/len(seq1)

    return same_percentage


#############################################################################
# id_at_pos
# calculates the rate of getting a similar nuc at a nuc site post DNA
# alignment with respect to its position on a codon
#############################################################################
def sim_at_pos(ref):
    count_1 = sim_rate_1 = count_2 = sim_rate_2 = count_3 = sim_rate_3 = 0

    for pos in ref:
        if pos % 3 == 0:
            count_1 += 1
            sim_rate_1 += sim_map[pos]
        if (pos - 1) % 3 == 0:
            count_2 += 1
            sim_rate_2 += sim_map[pos]
        if (pos - 2) % 3 == 0:
            count_3 += 1
            sim_rate_3 += sim_map[pos]
    per_id = [sim_rate_1/count_1, sim_rate_2/count_2, sim_rate_3/count_3]

    return per_id


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
            # randomly mutate nuc to a new nucleotide
            new_nuc = rand_nuc(nuc)
            print("s1: old " + nuc + " new " + new_nuc)
            new_DNA_seq += new_nuc
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
def strategy2(DNA_seq, prot_seq):
    # compare original DNA to a mutated new sequence
    new_DNA_seq = ""

    for count, nuc in enumerate(DNA_seq, 1):
        # mutation rate apply to all 3rd position nuc
        if count % 3 == 0:
            if random.random() < MUTATION_RATE:
                # randomly mutate nuc to a new nucleotide (A, C, G, T)
                new_nuc = rand_nuc(nuc)
                print("s2: old "+nuc+" new "+new_nuc)
                new_DNA_seq += new_nuc
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