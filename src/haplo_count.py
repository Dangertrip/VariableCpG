import sys
import subprocess
import pysam


def main():
    haplofile = sys.argv[1]
    lookup_num = int(sys.argv[2])
    if lookup_num==2:
        lookup_table = {'CC':0,'CT':1,'TC':2,'TT':3}
    if lookup_num==3:
        lookup_table = {'CCC':0,'CCT':1,'CTC':2,'CTT':3,'TCC':4,'TCT':5,'TTC':6,'TTT':7}
    if lookup_num==4:
        lookup_table = {'CCCC':0,'CCCT':1,'CCTC':2,'CCTT':3,'CTCC':4,'CTCT':5,'CTTC':6,'CTTT':7,'TCCC':8,'TCCT':9,'TCTC':10,'TCTT':11,'TTCC':12,'TTCT':13,'TTTC':14,'TTTT':15}    
    

    #finish build the binary look up table
    dic={}
    with open(haplofile) as f:
        for line in f:
            temp = line.strip().split()
            cpg = temp[3]
            pos = temp[5].split(',')
            chr = temp[0]
            if len(cpg)<lookup_num:
                continue
            start = cpg[:lookup_num]
            for i in range(lookup_num,len(cpg)):
                p=chr
                for j in range(i-lookup_num,i):
                    #print(j)
                    p = p+'_'+pos[j]
                if p in dic:
                    dic[p][lookup_table[start]]+=1
                else:
                    dic[p] = [0]*(2**lookup_num)
                    dic[p][lookup_table[start]]+=1
                start=cpg[i-lookup_num+1:i+1]

    t=['']*(2**lookup_num)
    for key in lookup_table:
        s=''
        for k in key:
            if k=='C': s+='1'
            else: s+='0'
        t[lookup_table[key]]=s

    title='chr\tpos\t'+'\t'.join(t)+'\n'
    result=[title]
    for p in dic:
        temp = p.split('_')
        arr=list(map(lambda x:str(x),dic[p]))
        sarr = '\t'.join(arr)
        s=temp[0]+'\t'+','.join(temp[1:])+'\t'+sarr+'\n'
        result.append(s)
    with open(haplofile+'.haplo_comb_'+str(lookup_num),'w') as f:
        f.writelines(result)


                

if __name__=="__main__":
    #print(search(10542,[10469, 10471, 10484, 10489, 10493, 10497, 10525, 10542, 10563, 10571]))
    main()
