import sys

fn = sys.argv[1]
num = int(sys.argv[2])

num = 2**num
if num==8:
    lookup_table = {0:'111',1:'110',2:'101',3:'100',4:'011',5:'010',6:'001',7:'000'}
if num==16:
    lookup_table = {'CCCC':0,'CCCT':1,'CCTC':2,'CCTT':3,'CTCC':4,'CTCT':5,'CTTC':6,'CTTT':7,'TCCC':8,'TCCT':9,'TCTC':10,'TCTT':11,'TTCC':12,'TTCT':13,'TTTC':14,'TTTT':15}
def count1(filename):
    f=open(filename)
    result=[]
    for line in f:
        temp = line.strip().split()
        if temp[1]=='pos': continue
        chr = temp[0]
        s = temp[1]
        _sum = 0
        _count = list(map(lambda x:int(x),temp[2:num+2])) 
        if sum(_count)<3: continue
        count_num=0
        r=[chr,s]
        order=0
        for c in _count:
            if c>=2: 
                count_num+=1
                r.append(order)
                r.append(lookup_table[order])
            order+=1
        if count_num>1: continue
        result.append(r)
    return result

def specific_haplo(results):
    dic={}
    result_str=[]
    for r in results:
        temp = map(lambda x:'_'.join(x))
        result_str.append(temp)
        for t in temp:
            if t in dic:
                dic[t]=dic[t]+1
            else:
                dic[t]=1
    specific=[]
    for i in range(len(results)):
        r = results[i]
        rr = result_str[i]
        for j in range(len(r)):
            if dic[rr[j]]==

    

if __name__=="__main__":
    print(count1(fn))

        
    
