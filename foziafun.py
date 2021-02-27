def gc_contentCounter(seq):
     g = seq.count("G")
     c = seq.count("C")
     gc_content = ((g+c)/len(seq))*100
     return gc_content

##print(gc_contentCounter("ACGC"))

####    Second function for counting lines in a file 
    
def lineCounter(file):
    count = 0
    with open (file, "r") as file_in:
      for line in file_in:
           count = count+1
    return count

##f = "workshop_data.fasta"
##lineCounter(f)


### Write a function that validates a DNA seq

def dna_seqvaliadate(seq):
    c = seq.count("C")
    t = seq.count("T")
    g = seq.count("G")
    a = seq.count("A")
    lenght = len(seq)
    Dna = c+t+g+a
    sequence_validation = Dna == lenght
    return(sequence_validation)

#print(dna_seqvaliadate("seq"))

seq_a = "TATAAATTTCTA"
seq_b = "TTTCTATAAATT"
def missMatch(seq_a, seq_b):
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()

     
    count = 0
    length = len(seq_a)
    for base in range(0,length):
        if seq_a[base] != seq_b[base]:
            count = count + 1
    return count
    
##print(missMatch(seq_a, seq_b))
##    
##primer = "TTTC"
##seq = "TTTCTATAAATTTCTATAATGGGGGCCCCCCTT"
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



def primerAttNo(primer, seq):
     seq = seq.upper()
     primer = primer.upper()
     print(primer)
     kmersList = kmers(seq, primer)
     count = 0
     for kmer in kmersList:
          diff = missMatch(kmer,primer)
          if diff < 2:
               count = count + 1
     return count


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




       
       
   
