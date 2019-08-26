import sys
import subprocess
import pysam
import os
from utils import binary_search

def load_cpg_index(cpg_file_path):
    '''
    Generate cpg index.
    Param:
        cpg_file_path: string, file name of cpg bed file.
    Return:
        cpg_index: dicionary, {'chromsome_name(chr1)':[1,2,3(cpg_position)]}
    '''
    index = {}
    with open(cpg_file_path) as f:
        for line in f:
            chrom, cpg_start, cpg_end = line.strip().split()[:3]
            if chrom in index:
                index[chrom].append(int(cpg_start))
            else:
                index[chrom]=[int(cpg_start)]
    
    for chrom in index:
        index[chrom].sort()
    return index
    

def get_certain_range_cpgs(cpg_starts_pos, cpg_list_chrom, start, end):
    '''
    Get cpgs between start and end in certain chrom
    Param:
        cpg_starts_pos: int, Location(index) of starting cpg in cpg_list_chrom.
        cpg_list_chrom: int[], cpg list of certain chrom. Numbers are the starting points of cpgs at that chromsome.
        start: int, read starting point.
        end: int, read ending point.
    '''
    cpg_arr = []
    for i in range(cpg_starts_pos, cpg_starts_pos+100):
        # need to test whether marginal cpg will be included(the first base and the last base)
        if i==len(cpg_list_chrom) or  cpg_list_chrom[i]>end:
            break
        num=cpg_list_chrom[i]
        if num>=start and num<end:
            cpg_arr.append(num)
    return cpg_arr


def extract_cpgs_from_read(reads, read_start, cpg_list, strand_specific_table):
    '''
    Extract all cpg status from one read.
    Param:
        reads: string, reads information: ATCGATGTCGTA......
        read_start: int, read starting position.
        cpg_list: int[], list of cpg site position on genome, only cpg on reads will be included.
        strand_specific_table: string[], marker for deciding whether cpg_status should be C or T.

    '''
    cpg_status=''
    cpg_pos=''
    seq_err=False
    for pos in cpg_list: 
        relative_pos = pos - read_start
        if relative_pos<0 or relative_pos>len(reads): continue
        # This is about the marginal CpG processing. I should double check the marginal process in testing.
        cpg_in_reads=reads[relative_pos:relative_pos+2]
        if cpg_in_reads==strand_specific_table[0]:
            cpg_status=cpg_status+'C'
            cpg_pos = cpg_pos+str(pos)+','
        else:
            if cpg_in_reads==strand_specific_table[1]:
                cpg_status=cpg_status+'T'
                cpg_pos = cpg_pos+str(pos)+','
                #Position should be different in +/- strand for cpg. I will deal with this later.
            else:
                #print(cpg_in_reads,strand_specific_table,num,s)
                seq_err=True
                break
    return cpg_status, cpg_pos[:-1], seq_err # remove ',' from cpg_pos
 

def get_haplotype(cpg_index, bam_file_path):
    '''
    Generate bam_file_path+'.status' as haplotype file.
    Param:
        cpg_index: cpg_index from load_cpg_index
        bam_file_path: string, bam file path.
    '''
    
    os.system('rm '+bam_file_path+'.status')
    strand_map={
            0:{'+':['CG','TG'],'-':['CG','CA']},
            1:{'-':['CG','TG'],'+':['CG','CA']}
    }
    #'+'(++,+-) or '-'(-+,--) means strand, +CG,-CG means methylated +C/-C; +TG/-CA means unmethylated +C/-C
    #position of C in CG in different strand do not change! Because bsmap transfer all reads into ++ and __ which all 
    #in the + strand. 


    bam_file = pysam.AlignmentFile(bam_file_path, "r")
    result=[]
    lastreads=[]
    dup_dic=set()

    for line in bam_file:
        bam_line = line.to_string().strip().split()
        name = bam_line[0]
        reverse_flag = int((int(bam_line[1])&(0x10))>0) # 0x10 SEQ being reverse complemented. 
        chrom = bam_line[2]
        read_start = int(bam_line[3])
        reads = bam_line[9]
        read_end = read_start+len(reads)
        quality = bam_line[10]
        mismatch = bam_line[11][bam_line[11].rfind(':'):]
        strand = bam_line[12][-2:]

        # We don't use chrM reads or reads in undifined chromsome in our haplotype analysis pipeline.
        if chrom == 'chrM' or chrom not in cpg_index: continue
        if strand[0]!='+' and strand[0]!='-' or strand[1]!='+' and strand[1]!='-':
            strand = bam_line[13][-2:]

        cpg_starting_pos = binary_search(read_start,cpg_index[chrom])

        cpg_list = get_certain_range_cpgs(cpg_starting_pos, cpg_index[chrom], read_start, read_end)

        strand_specific_table = strand_map[reverse_flag][strand[1]]
        cpg_status, cpg_pos, seq_err = extract_cpgs_from_read(reads, read_start, cpg_list, strand_specific_table)
        
        if seq_err: continue
        
        if cpg_status:

            duplicate_param = [chrom, str(read_start), str(read_end), cpg_status, strand, cpg_pos]
            duplicate_marker = ''.join(duplicate_param)
            
            if duplicate_marker in dup_dic:
                continue
            else:
                dup_dic.add(duplicate_marker)
                lastreads.append(duplicate_marker)
                if len(lastreads)>30:
                    remove_s = lastreads[0]
                    lastreads=lastreads[1:]
                    dup_dic.remove(remove_s)
            output_param = [chrom, str(read_start), str(read_end), cpg_status, strand, cpg_pos, str(reverse_flag), mismatch]
            output_str = '\t'.join(output_param) + '\n'
            result.append(output_str)
            if len(result)>100000:
                with open(bam_file_path+'.status','a') as ff:
                    ff.writelines(result)
                result=[]
    bam_file.close()
    with open(bam_file_path+'.status','a') as ff:
        ff.writelines(result)


if __name__=="__main__":
    #print(search(10542,[10469, 10471, 1048 10489, 10493, 10497, 10525, 10542, 10563, 10571]))
    
    cpg_file_path = "/data/yyin/data/ref/cpg/hg19_cpg.bed"
    bam_file_path = "TestCaseGenerator/TestData_50000_hg19.fastq.bam"

    cpg_index = load_cpg_index(cpg_file_path)

    #finish build the binary look up table

    get_haplotype(cpg_index, bam_file_path)


