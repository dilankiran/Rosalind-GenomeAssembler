# Genome Assembly as Shortest Superstring problem solution

# Problem:
# For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.

# By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

# Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).

# The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

# Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).

# Sample Data Set:
# >Rosalind_56
# ATTAGACCTG
# >Rosalind_57
# CCTGCCGGAA
# >Rosalind_58
# AGACCTGCCG
# >Rosalind_59
# GCCGGAATAC

# Sample Output: 
# ATTAGACCTGCCGGAATAC

#########################

# Solution

conda install -c conda-forge biopython #installing biopython

from Bio import SeqIO # importing standard Sequence Input/Output interface
import io # python's facility to deal with input/outputs

def construct_superstring(reads_list, superstring=''): 
    if len(reads_list) == 0:
        return superstring

    elif len(superstring) == 0:
        superstring = reads_list.pop(0)
        return construct_superstring(reads_list, superstring)

    else:
        for current_read_index in range(len(reads_list)):
            current_read = reads_list[current_read_index]
            current_read_length = len(current_read)

            for trial in range(current_read_length // 2):
                overlap_length = current_read_length - trial

                if superstring.startswith(current_read[trial:]):
                    reads_list.pop(current_read_index)
                    return construct_superstring(reads_list, current_read[:trial] + superstring)

                if superstring.endswith(current_read[:overlap_length]):
                    reads_list.pop(current_read_index)
                    return construct_superstring(reads_list, superstring + current_read[overlap_length:])


if __name__ == "__main__": 

    given_dataset = ">Rosalind_56 \nATTAGACCTG\n>Rosalind_57\nCCTGCCGGAA\n>Rosalind_58\nAGACCTGCCG\n>Rosalind_59\nGCCGGAATAC\n"
    handle = io.StringIO(given_dataset) # defined as string io 

    reads = []
    for record in SeqIO.parse(handle, 'fasta'):
        reads.append(str(record.seq)) # reading as str
    handle.close()
    result = construct_superstring(reads) 
    print(result)
