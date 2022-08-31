Genome Assembly as Shortest Superstring problem solution

Problem:
For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.

By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).

The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).



Sample Data Set:
>Rosalind_56
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC

Sample Output: 
ATTAGACCTGCCGGAATAC

----------------------------------------------------------------

Solution

conda install -c conda-forge biopython (ilk olarak biopyhton yuklenmeli)

from Bio import SeqIO # importing standard Sequence Input/Output interface
import io # pyhton's facility to deal with input/outputs

def superstring_olusturma(reads_list, superstring=''): #fonksiyon yazılır
    if len(reads_list) == 0:
        return superstring

    elif len(superstring) == 0:
        superstring = reads_list.pop(0)
        return superstring_olusturma(reads_list, superstring)

    else:
        for mevcut_read_index in range(len(reads_list)):
            mevcut_read = reads_list[mevcut_read_index]
            mevcut_read_uzunlugu = len(mevcut_read)

            for trial in range(mevcut_read_uzunlugu // 2):
                overlap_uzunlugu = mevcut_read_uzunlugu - trial

                if superstring.startswith(mevcut_read[trial:]):
                    reads_list.pop(mevcut_read_index)
                    return superstring_olusturma(reads_list, mevcut_read[:trial] + superstring)

                if superstring.endswith(mevcut_read[:overlap_uzunlugu]):
                    reads_list.pop(mevcut_read_index)
                    return superstring_olusturma(reads_list, superstring + mevcut_read[overlap_uzunlugu:])


if __name__ == "__main__": #datasetini problemde verildigi sekilde tanimliyoruz

    verilen_dataset = """>Rosalind_56 
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC
"""
    handle = io.StringIO(verilen_dataset) #string io olarak tanimlanir 

    reads = []
    for record in SeqIO.parse(handle, 'fasta'):
        reads.append(str(record.seq)) #str olarak okutuyoruz
    handle.close()
    sonuc = superstring_olusturma(reads) #sonucumuzu tanimliyoruz
    print(sonuc) #sonucu yazdirmak icin


