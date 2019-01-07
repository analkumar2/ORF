# This file writes without the peptide sequence and with only the first gene instance numbered
# Make changes to line 11, 114, and 129 according to where your files are stored
# For making a fasta file, uncomment lines 130-137
# For numbering every instance of the TTG found, comment out line 121 and replace the '1' in line 122 to '0'.


import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

list_of_genes = list(SeqIO.parse("C:\Google Drive\study\Computers\Anup bhaiya ORF program\orf_genomic_1000_all.fasta", "fasta"))

    

class gene:                                         #defines a class 'gene' which has the attributes 'position', 'context' etc.
    def __init__(self, gene_idd, seqq):
        self.gene_idd = gene_idd
        self.seqq = seqq
    def gene_id(self):
        return self.gene_idd
    def ORF_seq(self):
        return self.seqq[1001:]
    def ORF_length(self):
        return len(self.ORF_seq())
    def UTR_seq(self):
        return self.seqq[:1000]
    def UTR_length(self):
        return len(self.UTR_seq())
    def position(self, string = 'TTG', start = -50):
        list_of_positions = []
        for gggg in re.finditer(string, self.UTR_seq()[start:]):
            if gggg.start()!=-1:
                list_of_positions.append(gggg.start()+start)
        return list_of_positions[::-1]
    def isitkozak(self):
        list_of_isitkozak =[]
        for i in self.position():
            if i == -3 and self.UTR_seq()[-6:-3]=='AAA':
                list_of_isitkozak.append('Y')
            elif i ==-4 and self.UTR_seq()[-7]=='A' and self.UTR_seq()[-1]=='T':
                list_of_isitkozak.append('Y')
            elif self.UTR_seq()[i-3]=='A' and self.UTR_seq()[i+3]=='T':
                list_of_isitkozak.append('Y')
            elif self.UTR_seq()[i-3:i]=='AAA':
                list_of_isitkozak.append('Y')
            else:
                list_of_isitkozak.append('N')
        return list_of_isitkozak[::-1]
    def context(self):
        list_of_contexts = []
        for i in self.position():
            if i == -3 or i == -4:
                list_of_contexts.append(self.UTR_seq()[i-3:])
            else:
                list_of_contexts.append(self.UTR_seq()[i-3:i+4])
        return list_of_contexts[::-1]
    def IO(self):
        list_of_IO = []
        for i in self.position():
            if i%3 ==0:
                list_of_IO.append('I')
            else:
                list_of_IO.append('O')
        return list_of_IO[::-1]            
    def peptide_length(self, tostop = True):
        list_of_peptide_length=[]
        for i in self.translated_peptide_sequence(tostop):
            list_of_peptide_length.append(len(i))
        return list_of_peptide_length[::-1]
    def translated_peptide_sequence(self, tostop = True):
        list_of_peptide = []
        for i in range(0,len(self.position())):
            peptide = Seq(self.seqq[1000+self.position()[i]:],IUPAC.unambiguous_dna)
            list_of_peptide.append(str(peptide.translate(to_stop=tostop)))
        return list_of_peptide[::-1]
    def isuORF(self):
        list_of_isuORF =[]
        for i in range(0,len(self.position())):
            if self.peptide_length()[i]<self.peptide_length(False)[i]:
                list_of_isuORF.append('Y')
            else:
                list_of_isuORF.append('N')
        return list_of_isuORF[::-1]
    def oORF(self):
        list_of_oORF=[]
        for i in range(0,len(self.isuORF())):
            if self.IO()[i]=='O' and self.isuORF()[i]=='N':
                list_of_oORF.append('Y')
            else:
                list_of_oORF.append('N')
        return list_of_oORF[::-1]
    def pre_termination(self):
        pass
    def parallel_to_ORF(self):
        pass
        
list_of_object_genes=[]                                                            #stores the information

for i in range(0,len(list_of_genes)):                                    #looping to separate the gene_id and the sequence from fasta file and making them class 'gene' objects
    oo = gene(list_of_genes[i].description.split()[1] , str(list_of_genes[i].seq))
    list_of_object_genes.append(oo)
    

#for i in range(0,len(list_of_object_genes)):
#    if len(list_of_object_genes[i].position()) > 1:
#        print('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<}    {:<}    {:<}\n'.format(i+1, list_of_object_genes[i].gene_id(), list_of_object_genes[i].position()[0], list_of_object_genes[i].context()[0], list_of_object_genes[i].IO()[0], list_of_object_genes[i].translated_peptide_sequence()[0], list_of_object_genes[i].isitkozak()[0], list_of_object_genes[i].kozakcontext()[0]))
#        for k in range(1,len(list_of_object_genes[i].position())):
#            print('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<}    {:<}    {:<}\n'.format('' ,'' , list_of_object_genes[i].position()[k], list_of_object_genes[i].context()[k], list_of_object_genes[i].IO()[k], list_of_object_genes[i].translated_peptide_sequence()[k], list_of_object_genes[i].isitkozak()[k], list_of_object_genes[i].kozakcontext()[k]))
#    elif len(list_of_object_genes[i].position()) == 1:
#        print('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<}    {:<}    {:<}\n'.format(i+1, list_of_object_genes[i].gene_id(), list_of_object_genes[i].position()[0], list_of_object_genes[i].context()[0], list_of_object_genes[i].IO()[0], list_of_object_genes[i].translated_peptide_sequence()[0], list_of_object_genes[i].isitkozak()[0], list_of_object_genes[i].kozakcontext()[0]))            #finds position of 'TTG' -50 nucleotides upstream
#    else:
#        print('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<}    {:<}    {:<}\n'.format(i+1, list_of_object_genes[i].gene_id(), '--', '--', '--', '--', '--', '--'))

output_text_file = open("C:\Google Drive\study\Computers\Anup bhaiya ORF program\ORF outputwithout peptide seq.txt",'w')

output_text_file.write('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<14}    {:<23}    {:<10}    {:<5}\n'.format('S.No.', 'Gene ID', 'Position of TTG', 'Context of TTG', 'I\O', 'Peptide length', 'Is it a Kozak Sequence?', 'Is it uORF', 'oORF'))
for i in range(0,len(list_of_object_genes)):
    if len(list_of_object_genes[i].position()) > 1:
        output_text_file.write('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<14}    {:<23}    {:<10}\n'.format(i+1, list_of_object_genes[i].gene_id(), list_of_object_genes[i].position()[0], list_of_object_genes[i].context()[0], list_of_object_genes[i].IO()[0], list_of_object_genes[i].peptide_length()[0], list_of_object_genes[i].isitkozak()[0], list_of_object_genes[i].isuORF()[0]))
        for k in range(1,len(list_of_object_genes[i].position())):
            output_text_file.write('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<14}    {:<23}    {:<10}    {:<5}\n'.format('' ,list_of_object_genes[i].gene_id() , list_of_object_genes[i].position()[k], list_of_object_genes[i].context()[k], list_of_object_genes[i].IO()[k], list_of_object_genes[i].peptide_length()[k], list_of_object_genes[i].isitkozak()[k], list_of_object_genes[i].isuORF()[k], list_of_object_genes[i].oORF()[k]))
    elif len(list_of_object_genes[i].position()) == 1:
        output_text_file.write('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<14}    {:<23}    {:<10}    {:<5}\n'.format(i+1, list_of_object_genes[i].gene_id(), list_of_object_genes[i].position()[0], list_of_object_genes[i].context()[0], list_of_object_genes[i].IO()[0], list_of_object_genes[i].peptide_length()[0], list_of_object_genes[i].isitkozak()[0], list_of_object_genes[i].isuORF()[0], list_of_object_genes[i].oORF()[0]))            #finds position of 'TTG' -50 nucleotides upstream
    else:
        output_text_file.write('{:<5}    {:<9}    {:<15}    {:<14}    {:<5}    {:<14}    {:<23}    {:<10}    {:<5}\n'.format(i+1, list_of_object_genes[i].gene_id(), '--', '--', '--', '--', '--', '--', '--'))    

output_text_file.close()

#fastafile = open("C:\Google Drive\study\Computers\Anup bhaiya ORF program\infast.fasta",'w')
#recordobject=[]
#for i in range(0, len(list_of_object_genes)):
#    for k in range(0,len(list_of_object_genes[i].position())):
#        recordobject.append(SeqRecord(Seq(list_of_object_genes[i].translated_peptide_sequence()[k], IUPAC.protein), list_of_object_genes[i].gene_id(), description = list_of_object_genes[i].IO()[k] ))
#        
#SeqIO.write(recordobject, fastafile, "fasta")
#fastafile.close()        
        
        
        
        
        
        
        
        