import sys
import numpy as np

def read_file(name):
    with open(name) as f:
        lines = f.readlines()
    dict = {}
    for line in lines[1:]:
        line_content = line.strip().split()
        vector = list(map(lambda x:int(x),line_content[2:]))
        vector_marker = line_content[0]+'_'+line_content[1]
        dict[vector_marker] = vector
    return dict

def common_status(names):
    lookup_table = {0:'CCC',1:'CCT',2:'CTC',3:'CTT',4:'TCC',5:'TCT',6:'TTC',7:'TTT'}
    wt = read_file(names[0])
    cancer = read_file(names[1])
    label = names[2]
    result = ['name\twt_status\tcancer_status\tdepth_wt\tdepth_cancer\tstatus_number_wt\tstatus_number_cancer\twt_most_common\tcancer_most_common\tcommon_table\tcommon_sum\n']
    for marker in wt:
        if marker in cancer:
            v_wt = np.array(wt[marker])
            v_cancer = np.array(cancer[marker])
            status_number_wt = np.sum(v_wt>0)
            status_number_cancer = np.sum(v_cancer>0)
            s_wt = np.sum(v_wt)
            s_cancer = np.sum(v_cancer)
            v_wt = v_wt/s_wt
            v_cancer = v_cancer/s_cancer
            common = np.abs(v_wt - v_cancer)
            s_common = np.sum(common)
            # s_wt = sum(v_wt)
            # s_cancer = sum(v_cancer)
            # s_common = sum(common)
            v_wt_str = ','.join(list(map(lambda x:str(x),v_wt)))
            v_cancer_str = ','.join(list(map(lambda x:str(x),v_cancer)))
            common_str = ','.join(list(map(lambda x:str(x),common)))
            r = [marker,v_wt_str,v_cancer_str,str(s_wt),str(s_cancer),str(status_number_wt),str(status_number_cancer),np.argmax(v_wt),np.argmax(v_cancer),common_str,s_common]
            r = '\t'.join(list(map(lambda x:str(x),r)))+'\n'
            result.append(r)
    with open(label+'common.3CpG.status.txt','w') as f:
        f.writelines(result)


def common_counts(names):
    wt = read_file(names[0])
    cancer = read_file(names[1])
    label = names[2]
    result = ['name\twt_total_count\tcancer_total_count\tcommon_count\twt_ratio\tcancer_ratio\n']
    for marker in wt:
        if marker in cancer:
            v_wt = wt[marker]
            v_cancer = cancer[marker]
            common = []
            for i in range(8):
                common.append(min(v_wt[i], v_cancer[i]))
            s_wt = sum(v_wt)
            s_cancer = sum(v_cancer)
            s_common = sum(common)
            r = [marker,s_wt,s_cancer,s_common,float(s_common)/s_wt,float(s_common)/s_cancer]
            r = '\t'.join(list(map(lambda x:str(x),r)))+'\n'
            result.append(r)
    with open(label+'common.3CpG.txt','w') as f:
        f.writelines(result)

if __name__=="__main__":
    common_status(sys.argv[1:])

