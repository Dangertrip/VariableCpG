import numpy as np

def vector_generator():
    sec = np.random.random_sample((7,))
    sec = np.sort(sec)
    vector = []
    lastnum = 0
    for s in sec:
        vector.append(s-lastnum)
        lastnum = s
    vector.append(1-lastnum)
    return np.array(vector)

def linear_distance(v1,v2):
    '''
    Compute the distacne between two vectors.
    Param:
        v1,v2: numpy.array, 1*N vector.
    return:
        linear distance between v1 and v2.
    '''
    return np.sqrt(np.sum(np.power(v1-v2,2)))

def target_CpG_methy_ratio(v):
    '''
    Compute the methylation of tri-bases.
    Since:
        {'CCC':0,'CCT':1,'CTC':2,'CTT':3,'TCC':4,'TCT':5,'TTC':6,'TTT':7}
    methylation of target CpG is v[0]+v[1]+v[4]+v[5]
    Param:
        v: numpy.array, 1*8 vector. 
    
    '''
    if (len(v)!=8):
        raise Exception('target_CpG_methy_ratio only support tri-bases!')
    return v[0]+v[1]+v[4]+v[5]

def information_entropy(v):
    return -1*np.sum(v*np.log(v, where=v>0))


def methy_entropy(v):
    ratio = target_CpG_methy_ratio(v)
    return -1*(ratio*np.log(ratio)+(1-ratio)*np.log(1-ratio))


def mas(v):
    return (information_entropy(v)+1)**methy_entropy(v)


def conditional_information_entropy(v):
    pass
    # methy_ratio  = target_CpG_methy_ratio(v) + 1e-10
    # entropy = v[0]*np.log(v[0]/methy_ratio)


if __name__=="__main__":
    testcases = [[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125],
                 [0,0,0,0,0,0,0,1],
                 [0.5,0,0,0.5,0,0,0,0]]
    for testcase in testcases:
        print(information_entropy(np.array(testcase)))
