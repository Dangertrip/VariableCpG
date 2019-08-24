import sys
import os

filename = sys.argv[1]
tempfile = filename + '.pairend.temp'
os.system('cat '+filename+' | sort -T ./ -k4,4 -k2n,2 >'+tempfile)
if os.path.exists(filename+'.pair.temp'):
    os.system('rm '+filename+'.pair.temp')
with open(tempfile) as f:
    lastreadsname=''
    lastline=''
    result=[]
    for line in f:
        if len(result)==0:
            result.append(line)
            continue
        new = line.strip().split()
        old = result[-1].strip().split()
        if new[3]==old[3]:#readsname
            if new[0]!=old[0] or abs(int(new[1])-int(old[1]))>5000:#reads chr and position
                result.append(line)
            else:
                s=old[0]+'\t'+old[1]+'\t'+new[2]+'\t'+old[3]
                dic={}
                arr=[]
                t_cpg = old[4]
                t_pos = old[6].split(',')
                for i in range(len(t_cpg)):
                    arr.append((t_cpg[i],t_pos[i]))
                    dic[t_pos[i]]=t_cpg[i]
                t_cpg = new[4]
                t_pos = new[6].split(',')
                mark=False
                for i in range(len(t_cpg)):
                    if t_pos[i] in dic:
                        if t_cpg[i]!=dic[t_pos[i]]:
                            mark=True
                            break
                    else:
                        arr.append((t_cpg[i],t_pos[i]))
                if mark:
                    continue
                cpg=''
                pos=''
                for c,p in arr:
                    cpg = cpg+c
                    pos = pos+p+','
                newline = old
                newline[2]=new[2]
                newline[4]=cpg
                newline[6]=pos[:-1]  #get rid of the last ,
                newline[5]=old[5][0]
                s='\t'.join(newline)+'\n'
                result[-1]=s
                if len(result)>1000000:
                    with open(filename+'.pair.temp','a') as ff:
                        ff.writelines(result)
                    result=[]
        else:
            result.append(line)

    if len(result)>0:
        with open(filename+'.pair.temp','a') as ff:
            ff.writelines(result)
os.system('sort -T /mt1/yyin/Heart_haplotype/haplo_all/ -k1,1 -k2n,2 '+filename+'.pair.temp >'+filename+'.pair.result')

