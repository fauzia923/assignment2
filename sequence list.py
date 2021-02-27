from foziafun import missMatch, kmers, getNamesSeqs, primerAttNo
from Bio.Seq import Seq


geneNames, genes = getNamesSeqs("test_sequences.fasta")
primerNames, primers = getNamesSeqs("primer.fasta")


oPrimer = primers[0]

primers.append(str(Seq(oPrimer).reverse_complement()))
primerNames.append(">reverse_complement_primer")

primers.append(oPrimer[::-1])
primerNames.append(">reverse_primer")


gene = genes[2]
geneName = geneNames[2]

atachNumber = []


for i in range(len(primers)):
    primer = primers[i]
##    print(primer)
##    print(gene)
    
    att = primerAttNo(primer, gene)
    print(att)
    

##    atachNumber.append(att)
##

    
##print(atachNumber)



##print(primerNames, primers)
##print(names, genes)
