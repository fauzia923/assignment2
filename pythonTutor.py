def reverseComp(dna):
    dna = dna.upper()
    comp = ''
    for i in dna:
        if i == "A":
            comp += "T"
        elif i == "T":
            comp += "A"
        elif i == "G":
            comp += "C"
        elif i == "C":
            comp += "G"

    return ''.join(reversed(comp))

def comp(dna):
    dna = dna.upper()
    comp = ''
    for i in dna:
        if i == "A":
            comp += "T"
        elif i == "T":
            comp += "A"
        elif i == "G":
            comp += "C"
        elif i == "C":
            comp += "G"

    return ''.join(comp)

def getNamesSeqs(file):
    nameList = []
    genList = []

    with open(file,"r") as file_in:
        for lines in file_in:
            if not lines.startswith(">"):
                sequence = lines.strip()
                genList.append(sequence)
            elif lines.startswith(">"):
                name = lines.strip()
                nameList.append(name)
    return nameList, genList

def kmers(seq, primer):
   seq = seq.upper()
   primer = primer.upper()
   lenSeq = len(seq)
   lenPrimer = len(primer) 
   kmers = []
   for i in range(0, lenSeq-lenPrimer+1):
       kmer = seq[i:i+lenPrimer]
       kmers.append(kmer)
   return kmers

def missMatch(seq_a, seq_b):
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    count = 0
    length = len(seq_a)
    for base in range(0,length):
        if seq_a[base] != seq_b[base]:
            count = count + 1
    return count


def ifComp(b1, b2):
    basesDict = {"A": "T",
                 "T": "A",
                 "G": "C",
                 "C": "G"}
    return b1 == basesDict[b2]




def compMissMatch(seq_a, seq_b):
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    count = 0
    length = len(seq_a)
    for i in range(0,length):
        seq_aBase = seq_a[i]
        seq_bBase = seq_b[i]
        if not ifComp(seq_aBase, seq_bBase):
            count += 1
    return count



def primerCompAttNo(primer, seq):
    seq = seq.upper()
    primer = primer.upper()
##     print(primer)
    kmersList = kmers(seq, primer)
    kmerNumber = 0
    count = 0
    for kmer in kmersList:
        diff = compMissMatch(kmer,primer)
        if diff <= 2:
            count = count + 1
    return count

  
    
def primerAttNo(primer, seq):
    seq = seq.upper()
    primer = primer.upper()
##     print(primer)
    kmersList = kmers(seq, primer)
    kmerNumber = 0
    count = 0
    for kmer in kmersList:
        diff = missMatch(kmer,primer)
        if diff <= 2:
            count = count + 1
    return count


def primerAttKmers(primer, seq):
    seq = seq.upper()
    primer = primer.upper()
##     print(primer)
    kmersList = kmers(seq, primer)
    attKmers = []
    for kmer in kmersList:
        diff = missMatch(kmer,primer)
        if diff <= 2:
            attKmers.append(kmer)
    return attKmers



def indsKmers(seq, primer):
    seq = seq.upper()
    primer = primer.upper()
    lenSeq = len(seq)
    lenPrimer = len(primer) 
    kmers = []
    inds = []
    for i in range(0, lenSeq-lenPrimer+1):
        kmer = seq[i:i+lenPrimer]
        diff = compMissMatch(kmer,primer)
        if diff <= 2:
            kmers.append(kmer)
            inds.append(i)
    return kmers, inds

######## Trying Primer attachment for 3 genes
##############################################


## getting genes and primers 




geneNames, genes = getNamesSeqs("test_sequences.fasta")
primerNames, primers = getNamesSeqs("primer.fasta")

oPrimer = primers[0] #Orignal Primer

primers.append(oPrimer[::-1])
primerNames.append(">reverse_primer")

primers.append(reverseComp(oPrimer))
primerNames.append(">reverse_complement_primer")

### Checking attachemnts 

for i in range(len(primers)):
    primer = primers[i]
    for index, gene in enumerate(genes):
        att = primerCompAttNo(primer, gene)
        print("primer", [i+1], ": gives, ", att, " attachments for gene", [index+1])






gene = genes[2]
geneName = geneNames[2]
primer = primers[0]
print(geneName)
print(gene)



kmersList, inds = indsKmers(gene, primer)
print(kmersList, inds)

print(primer.upper())

print("Mismatches:",compMissMatch(primer, kmersList[0]))

##
##
####for i in range(len(primers)):
####    primer = primers[i]
####    for index, gene in enumerate(genes):
####
######        print(len(gene))
######        print(index)
######        print(gene)
####        att = primerAttNo(primer, gene)
####        print("primer", [i+1], ": gives, ", att, " attachments for gene", [index+1])
####    
##    
##
##
####
####kmersList = kmers(gene, primer)
####print(len(kmersList))
####
######
####att = primerAttNo(primer, gene)  
####print(att)
##
####print(primers)
####
####    
######print(primers[0].upper())
########[print(len(genes[i])) for i in genes]
######for gene in genes:
######    print(len(gene))
######    
######print(primerAttKmers(primers[0].upper(), genes[2]))
######
######
##
