

      
#Reading FASTA sequences with one compact function
def read_FASTA(filename):
      with open(filename) as file:
            return [(part[0].split('|'),
                     part[2].replace('\r', ''),
                     part[2].replace('\n', ''))
                    for part in
                    [entry.partition('\n')
                     for entry in file.read().split('>')[1:]]]



def Merge(infilename):
    infilepath = Getfilepath(infilename)
    os.chdir(infilepath)

    infilename_fmt = fileformat(infilename)
    infilename_nf = filenamenoformat(Getfilename(infilename))
    outfilename = "%s_mg.%s" % (infilename_nf, infilename_fmt)
    with open(outfilename, "w") as fout:
        with open(infilename, "r") as fin:
            rfin = csv.reader(fin, delimiter="\t")
            for row in rfin:
                if "#" in "".join(row):
                    fout.writelines(row)
                    fout.writelines("\n")
                elif ">" in "".join(row):
                    fout.writelines(row)
                    fout.writelines(";")
                else:
                    fout.writelines(row)
                    fout.writelines("\n")


def gccontent(input):
    #input=raw_input(' Enter Sequence: ')
    count=float()
    for i in input:
        if i=='G' or i=='C':
            count+=1
    content=(count/len(input)*100)
    return content

def GC_from_FASTA(fasta):
      fasta=read_FASTA(fasta)
      GCcontent=[]
      name=[]
      for i in range(0, len(fasta)):
            GCcontent+= [gccontent(fasta[i][1])]
            name+= [fasta[i][0]]
      print name[GCcontent.index(max(GCcontent))]
      print str(max(GCcontent)) + str('%')

def overlap(fasta):
      '''
      This code will take multiple fasta sequences from a single file, and return the longest
      consensus sequence between them.
      '''
      inputs=read_FASTA(fasta)
      total=[]
      for i in range(0, len(inputs)):
            total+=[inputs[i][1]]
      start=len(min(total))
      trial=min(total)[0:start]
      collection=[]
      for i in range(0,start+1):
            for j in range(0,start+1):
                  collection+=[trial[j:i]]
      new_collect=[]
      for i in collection:
            if all(i in j for j in total) is True:
                  new_collect+=[i]
      print (max(new_collect))
      

def revcomp(rev):
    # type: (object) -> object
    output=''
    error=0
    place=[]
    x2=rev[::-1]
    x3=x2
    for i in x2:
        if i=='A':
            output+='T'
        if i=='T':
            output+='A'
        if i=='G':
            output+='C'
        if i=='C':
            output+='G'
        elif i!='C' and i!='G' and i!='A' and i!='T':
            output+='N'
            error+=1
            place+=[x3.find(i)]
            x3=x3[:x3.find(i)] + x3[(x3.find(i)+1):]

    #print 'Number of errors: ' + str(error) + ', at position(s): ' + str(place)
    #print 'Reverse complement: ' + str(output)
    return output




def consensusseq(x,y):
    count=0
    count2=0
    consensus=''
    for i in x:
        if x[count:count+15] not in y[0:15]:
            consensus+=x[count]
            count+=1
    for i in y:    
        consensus+=y[count2]
        count2+=1
      
    print ('Length of assembled sequence: ' + str(len(consensus)))
    return str(consensus)




def primerwalk(filename):
    seqx=fasta_reader(filename) #Create an empty list
    walk=[]
    for i in range(0,len(seqx)):
        walk+=[seqx[i][1]]
    newseq=''
    for i in range(0,len(walk)):
        newseq=consensusseq(newseq,walk[i])
    count=0
    null_list=[]
    for i in range(0,len(newseq)):
        if newseq[i]=='N':
            count+=1
            null_list+=[i]
    print ('Number of Ns: ' + str(count) + ', at positions: ' + str(null_list))
    return newseq
                     




bases = ['T', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))


#RNA to Protein
def translate(seq):
      peptide=''
      for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break
                
      return peptide

#RNA splicing
def RNA_splice(filename):
      '''This code will take a fasta file with the first sequence being the unspliced RNA,
      and the remaining sequences being the introns, it will remove those introns from the RNA,
      produce a final mRNA, and translate it'''
      RNA=read_FASTA(filename)
      mRNA=RNA[0][1]
      intron=[]
      for i in range(1,len(RNA)):
            intron+=[RNA[i][1]]
      for i in intron:
            if i in mRNA:
                  mRNA=mRNA.replace(i,'')
      return translate(mRNA)
      
      
      
##Get ORFS from FASTA with their length, strand, and coordinates

from Bio import SeqIO
#record = SeqIO.read("NC_005816.gb","genbank")
table = 11
min_pro_len = 100

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer


#orf_list1 = find_orfs_with_trans(seq1, table, min_pro_len)
#for start, end, strand, pro in orf_list:
#    print("%s...%s - length %i, strand %i, %i:%i" \
#          % (pro[:30], pro[-3:], len(pro), strand, start, end))


def Uracil(x):
    for i in range(0, len(x)):
        if x[i] == 'T':
            x = x.replace(x[i], 'U')


    return Uracil('CCTCTTTTCAAGCAGCAGTTCATCTGCCTACACCAATATGCC')

'''
Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A CONSENSUS string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
'''

from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio import motifs
instances = []
for entry in fasta_reader('rosalind_cons.txt'):  #example fasta file with multiple entries
    instances+=[str(entry.seq)]
instances2=[]
for i in range(0,len(instances)):
    instances2+=[Seq(instances[i])]
    
m = motifs.create(instances)
print (m.consensus)
print (str('A:') + str(m.counts['A']).replace(',',' '))
print (str('C:') + str(m.counts['C']).replace(',',' '))
print (str('G:') + str(m.counts['G']).replace(',',' '))
print (str('T:') + str(m.counts['T']).replace(',',' '))



def fasta_reader(filename):
    from Bio.SeqIO.FastaIO import FastaIterator
    input=[]
    with open(filename) as handle:
        for record in FastaIterator(handle):
            input+=[[str(record.id),str(record.seq)]]
    return input



def six_translate(seq):
    seqlist=[]
    for i in range(0,len(seq), 3):
        if seq[i:i+3]=='ATG':
            seqlist+=[[translate(seq[i:]),1]]
    for i in range(1,len(seq), 3):
        if seq[i:i+3]=='ATG':
            seqlist+=[[translate(seq[i:]),2]]
    for i in range(2,len(seq), 3):
        if seq[i:i+3]=='ATG':
            seqlist+=[[translate(seq[i:]),3]]
    revy=revcomp(seq)
    for i in range(0,len(revy), 3):
        if revy[i:i+3]=='ATG':
            seqlist+=[[translate(revy[i:]),-1]]
    for i in range(1,len(revy), 3):
        if revy[i:i+3]=='ATG':
            seqlist+=[[translate(revy[i:]),-2]]
    for i in range(2,len(revy), 3):
        if revy[i:i+3]=='ATG':
            seqlist+=[[translate(revy[i:]),-3]]
            
    unique_setlist=[]
    for x in range(0,len(seqlist)):
        if len(seqlist[x][0])>50:
            unique_setlist+=[seqlist[x]]
    return unique_setlist


# The following snippet will take from the generator and make a list that contains the ID and its associated sequence as an entry within the list
input=[] #Create an empty list
for entry in fasta_reader('rosalind_cons.txt'):  #example fasta file with multiple entries
    input+=[[str(entry.id),str(entry.seq)]]

for i in range(0,len(input)):
    for j in range(0,len(input)):
        if input[i][1][:3] == input[j][1][-3:] and input[i][1] != input[j][1]:
            print (str(input[j][0] + str(' ')+ str(input[i][0])))
    
    
def primerwalk2(walk):
    #This is a primer walk that already contains the list (called 'walk') with the sequences in order
    newseq=''
    for i in range(0,len(walk)):
        newseq=consensusseq(newseq,walk[i])
    count=0
    null_list=[]
    for i in range(0,len(newseq)):
        if newseq[i]=='N':
            count+=1
            null_list+=[i]
    #print ('Number of Ns: ' + str(count) + ', at positions: ' + str(null_list))
    return newseq


def de_novo_assembly(filename):
    input=fasta_reader('rosalind_cons.txt') #import the fasta file using the fasta_reader function
    walk=[] #create an empty list
    
    for i in range(0,len(input)):
        walk+=[input[i][1]] #add only the sequences to be assembled to the empty list
    assemble=[] 
    count=0
    for i in range(0,len(walk)):
        walk2=walk.copy()
        query=str(walk[i])
        if str(query) not in str(assemble):
            walk2.remove(str(query))
            #print(str('walk: ') + str(len(walk)))
            #print(str('walk2: ') + str(len(walk2)))
        mid_length = int(len(query)/2) #get the mid length of the first string in the list 'walk'
        if query[:mid_length] not in str(walk2):
            assemble+=[query]
            if query[mid_length:] not in str(walk2): #this will find the first sequence in the assembly because the first half of the string is not in any other string
                assemble.remove(query)
            
            
            
            

    print (str('walk:') + str(len(walk)))
    print (str('walk2:') + str(len(walk2)))            
    print (str(assemble))
    while len(assemble)<len(walk):
        for j in range(0,len(walk)): #This is to assemble to the remaining sequences after the first sequence is found. This is done by matching the last half of the string to the rest of the remaining list
            walk2=walk
            mid_length = int(len(walk[j])/2) #get the mid length of the first string in the list 'walk'
            if str(walk[j][mid_length:]) in str(assemble[-1]): #finds the string within j whose last half is in the current assembly
                query=str(walk[j])
                walk2.remove(walk2[j])
                assemble+=[query]
            
    print (str(assemble))
    return primerwalk2(assemble)
                    
        
   
