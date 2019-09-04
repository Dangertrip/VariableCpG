import numpy as np
import sys
from Experiments.utils import *

def bed_reader(filename, depth_limit=10):
    ratio = {} #chr_pos:ratio
    with open(filename) as f:
        for line in f:
            if line[0]=='#': continue
            line_content = line.strip().split()
            depth = int(line_content[4])
            if depth<depth_limit: continue
            r = float(line_content[3])
            pos = line_content[0]+'_'+str(int(line_content[1])+1)
            ratio[pos] = r
    return ratio

def status_reader(filename, methy_dict, depth_limit=10):
    result = []
    with open(filename) as f:
        for line in f:
            for line in f:
                if 'pos' in line: continue
                line_content = line.strip().split()
                methy_ratio = []
                chrom = line_content[0]
                cpg_pos = line_content[1].split(',')
                for cpg in cpg_pos:
                    r = methy_dict.get(chrom+'_'+cpg,-1)
                    if r<0:
                        break
                    methy_ratio.append(r)
                if len(methy_ratio)<3: continue
                vector = np.array(list(map(lambda x:int(x),line_content[2:10])))
                content = [chrom, line_content[1], np.array(methy_ratio), vector]
                result.append(content)
    return result

def main(bed_filename, status_filename):
    methy_ratio_dict = bed_reader(bed_filename)
    status_list = status_reader(status_filename, methy_ratio_dict)
    result = []
    for content in status_list:


            

