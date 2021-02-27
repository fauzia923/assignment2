from foziafun import missMatch, kmers
from Bio.Seq import Seq




primer = "TTTC"
seq = "TTTCTATAAATTTCTATAATGGGGGCCCCCCTT"


my_dna = Seq("TTTCTATAAATTTCTATAATGGGGGCCCCCCTTT")

print(my_dna)
print(my_dna.complement())

seqComp = my_dna.complement()

def primerAttNo(primer, seq):
     kmersList = kmers(seq, primer)
     count = 0
     for kmer in kmersList:
          diff = missMatch(kmer,primer)
          if diff < 2:
               count = count + 1
     return count


kmersList = kmers(seqComp, primer)          
totalBindings = primerAttNo(primer, seq)
totalKmers = len(kmersList)

accu = totalBindings/totalKmers
print(accu)







     
     #calculate mismatches for each kmer with primer and store it in a variable
     #add an if statement whcih checks if the distanc3 is less tahn 2.
     # +1 inn he counter
##     print(kmer) 
