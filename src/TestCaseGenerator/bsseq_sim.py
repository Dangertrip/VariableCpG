'''
usage:  python bsseq_sim.py readsname
Generate certain length bisulfite reads from certain genome.
2 files will be generated after every run. bed file is the original position of every reads.

Generated TestData_50000_hg19.fastq for testing. 50000 reads from hg19 genome.
TestData_50000_hg19.fastq.bam is geneerated by bsmap.
bamToBed + bedtools intersect with all hg19 CpG can find all CpGs on mapped reads, which
can be used to test the haplotype's correctness.
'''


import numpy as np
import math
import random
import sys

#==============================Variable setting===================================
readsname = sys.argv[1]
conversion_ratio = 0.96
methylation_ratio = 0.7
readsnum = 50000
readlength = 100
#================================================================================


#==============================Get Genome fasta dictory=========================
def fasta_loader(path='/data/dsun/ref/humanigenome/hg19.fa'):
    '''
    Given certain path of fasta file, load all genome into memory. No check of fasta format!
    param:
        path: string, path to fasta file.
    return:
        chr_content_dict: dict, {chr_name:chr_content}. For example:{'chr1':'aatcc....'}
        chr_names: string[], chromsome names in loading order.
        chr_length: int[], chromsome names in chr_names' order.
    '''
    chr_name = ''
    chr_content_dict = {}
    mark = False
    chr_content = ''
    genome_length = 0
    chr_length = []
    chr_names = []
    with open(path) as f:
        for line in f:
            line_content = line.strip()
            # In the begining of the chromsome, fasta has ">chr1"
            if line_content[0]=='>':
                if chr_name!='':
                    chr_content_dict[chr_name] = chr_content
                    chr_length.append(len(chr_content))
                    genome_length += chr_length[-1]
                    chr_content = ''
                chr_name = line_content[1:]
                chr_names.append(chr_name)
            else:
                chr_content += line_content
    chr_content_dict[chr_name] = chr_content
    chr_length.append(len(chr_content))
    genome_length += chr_length[-1]
    return chr_content_dict, chr_names, chr_length


def reverse(read):
    dic={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    r=''
    for rr in read:
        r=dic[rr.upper()]+r
    return r


def bisulfite(conv_r,meth_r,read):
    r=''
    l = len(read)
    for i in range(1,l-1):
        base=read[i].upper()
        if base=='C':
            c = random.random()
            if c<conv_r:
                if read[i+1].upper()=='G':
                    m = random.random()
                    if m>meth_r:
                        base='T'
                else:
                    base='T'
        r=r+base
    return r

if __name__=="__main__":
    chr_content_dict, chr_names, chr_length = fasta_loader()
    genome_length = sum(chr_length)
    finalbed=[]
    for i in range(readsnum):
        read_start_pos = random.randint(0,genome_length)
        chr=0
        f=0
        t=0
        while read_start_pos>chr_length[chr]:
            read_start_pos -= chr_length[chr]
            chr+=1
        start = read_start_pos-1
        end = read_start_pos + readlength + 1
        if start<1: continue
        if end>chr_length[chr]: continue
        read = chr_content_dict[chr_names[chr]][start:end]
        r = read
        a1 = random.random()
        a2 = random.random()
        if a1>0.5:
            r=reverse(read)# Get reads from +/- strand
        r = bisulfite(conversion_ratio,methylation_ratio,r)
        if a2>0.5:
            r=reverse(r) # PCR +/-
        quality = 'E' * readlength
        print('@'+str(i)+'_'+readsname)
        print(r)
        print('+')
        print(quality)
        finalbed.append(chr_names[chr]+'\t'+str(start+1)+'\t'+str(end-1)+'\t'+str(i)+'_'+readsname+'\t'+str(f)+'\t'+str(t)+'\n')
    with open(readsname+'_simulation.bed','w') as f:
        f.writelines(finalbed)

