import sys

filename = sys.argv[1]
num = int(sys.argv[2])

num = 2**num


#CC,CT  CC,TC   CC,TT   CT,TC   CT,TT   TC,TT
split_count=[0]*6

count_arr=[0]*num
homocount=0
homocount_meth=0
count2=[]
f=open(filename)
for line in f:
    temp = line.strip().split()
    if temp[1]=='pos': continue
    chr = temp[0]
    s = temp[1]
    _sum=0
    count=0
    for i in range(2,num+2):
        t = int(temp[i])
        temp[i]=t
        if t==0:
            continue
        _sum+=t
        if t>1: count+=1
    #print(sum,count)
    if _sum>=3 and count>0:
        #print(line,count_arr,count)
        if count>1:
            tri_start,_,tri_end=temp[1].split(',')
            count2.append(temp[0]+'\t'+tri_start+'\t'+tri_end+'\n')
        count_arr[count-1]+=1
        if count==1:
            if temp[2]!=0 or temp[-1]!=0: homocount+=1
            if temp[2]!=0: homocount_meth+=1
        if num==4 and count==2:
            if (temp[2]*temp[3]!=0):
                split_count[0]+=1
            if (temp[2]*temp[4]!=0):
                split_count[1]+=1
            if (temp[2]*temp[5]!=0):
                split_count[2]+=1
            if (temp[3]*temp[4]!=0):
                split_count[3]+=1
            if (temp[3]*temp[5]!=0):
                split_count[4]+=1
            if (temp[4]*temp[5]!=0):
                split_count[5]+=1
if sum(count_arr)!=0:
	result=[filename[:filename.find('.')]]
	result.extend(count_arr)
	result.append(count_arr[0]/float(sum(count_arr)))
	result.append((homocount)/float(count_arr[0]))
	result.append(homocount_meth/float(homocount))
	if num==4: result.extend(split_count)
	rr = '\t'.join(list(map(lambda x:str(x),result))) 
	print(rr)
#print(count_arr)
#print(count_arr[0]/float(sum(count_arr)))
#print((homocount)/float(count_arr[0]))
#print(homocount_meth/float(homocount))
f.close()
with open(filename+'.count2','w') as f:
    f.writelines(count2)
