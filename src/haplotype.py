import sys
import subprocess
import pysam
import os

def search(s,arr):
    l=0
    r=len(arr)-1
    #print(s,arr)
    while l<r:
        mid=(l+r)//2
        if arr[mid]<s:
            l = mid+1
            continue
        if arr[mid]>=s:
            r = mid-1
        #print(l,r)
    return (l+r)//2

def main():
    '''
    target_flanking_len=80
    rcTable={}
    rcTable{'A'}='T'
    rcTable{'T'}='A'
    rcTable{'G'}='C'
    rcTable{'C'}='G'
    rcTable{'N'}='N'
    rcTable{'R'}='Y'
    rcTable{'Y'}='R'
    rcTable{'M'}='K'
    rcTable{'K'}='M'
    rcTable{'S'}='S'
    rcTable{'W'}='W'
    '''
    strand_map={
            0:{'+':['CG','TG'],'-':['CG','CA']},
            1:{'-':['CG','TG'],'+':['CG','CA']}
            }
    #'+'(++,+-) or '-'(-+,--) means strand, +CG,-CG means methylated +C/-C; +TG/-CA means unmethylated +C/-C
    #position of C in CG in different strand do not change! Because bsmap transfer all reads into ++ and __ which all 
    #in the + strand. 


    cpgfile = sys.argv[1]
    bamfile = sys.argv[2]
    fn=bamfile
    with open(cpgfile) as f:
        lines = f.readlines()
    dic={}
    os.system('rm '+fn+'.status')
    for line in lines:
        chr,s,e = line.strip().split()[:3]
        if chr in dic:
            dic[chr].append(int(s))
        else:
            dic[chr]=[int(s)]
    
    for chr in dic:
        dic[chr].sort()
    

    #finish build the binary look up table

    f=pysam.AlignmentFile(fn, "r")
    result=[]
    lastreads=[]
    dup_dic=set()
    for line in f:
        temp = line.to_string().strip().split()
        name = temp[0]
        flag = int((int(temp[1])&(0x10))>0)
        chr = temp[2]
        if chr=='chrM':
            continue
        s = int(temp[3])
        reads = temp[9]
        qua = temp[10]
        mismatch = temp[11][temp[11].rfind(':'):]
        strand = temp[12][-2:]
        if strand[0]!='+' and strand[0]!='-' or strand[1]!='+' and strand[1]!='-':
            strand = temp[13][-2:]
        if not chr in dic: continue
        pos=search(int(s),dic[chr])
        e = s+len(reads)
        cpg_arr=[]
        for i in range(pos,pos+100):
            if i==len(dic[chr]) or  dic[chr][i]>e:
                break
            num=dic[chr][i]
            if num>=s and num<e:
                cpg_arr.append(num)
        #print(cpg_arr)
        #print(s)
        #print(reads)
        cpg_status=''
        containc=0
        containt=0
        cpg_pos=''
        #Here's a new method about strand specific haplotype
        strand_specific_table = strand_map[flag][strand[1]]
        seq_err=False
        for num in cpg_arr: 
            if num-s+1>=len(reads) or num-s<=1:
                continue
            cpg_in_reads=reads[num-s:num-s+2]
            if cpg_in_reads==strand_specific_table[0]:
                cpg_status=cpg_status+'C'
                cpg_pos = cpg_pos+str(num)+','
                containc=1
            else:
                if cpg_in_reads==strand_specific_table[1]:
                    cpg_status=cpg_status+'T'
                    cpg_pos = cpg_pos+str(num)+','
                    containt=1
                    #Position should be different in +/- strand for cpg. I will deal with this later.
                else:
                    #print(cpg_in_reads,strand_specific_table,num,s)
                    seq_err=True
                    break
        #print(cpg_status,seq_err)
        #print(temp,flag)
        if seq_err:
            continue
        if cpg_status:
            dup_s=chr+str(s)+str(e-1)+cpg_status+strand+cpg_pos[:-1]
            if dup_s in dup_dic:
                continue
            else:
                dup_dic.add(chr+str(s)+str(e-1)+cpg_status+strand+cpg_pos[:-1])
                lastreads.append(dup_s)
                if len(lastreads)>30:
                    remove_s = lastreads[0]
                    lastreads=lastreads[1:]
                    dup_dic.remove(remove_s)
            result.append(chr+'\t'+str(s)+'\t'+str(e-1)+'\t'+name+'\t'+cpg_status+'\t'+strand+'\t'+cpg_pos[:-1]+'\t'+str(flag)+'\t'+mismatch+'\n')
            if len(result)>100000:
                with open(fn+'.status','a') as ff:
                    ff.writelines(result)
                result=[]
    f.close()
    with open(fn+'.status','a') as ff:
        ff.writelines(result)
        #print(fn,sum,mix,mix/sum)


if __name__=="__main__":
    #print(search(10542,[10469, 10471, 10484, 10489, 10493, 10497, 10525, 10542, 10563, 10571]))
    main()
