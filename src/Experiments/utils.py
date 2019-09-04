import numpy as np

def fixed_vector_generator():
    pass


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


def prob_status(methylation_ratio, status):
    # 111     110     101     100     011     010     001     000
    status_map = {0:'111', 1:'110', 2:'101', 3:'100', 4:'011', 5:'010', 6:'001', 7:'000'}
    s = status_map[status]
    status_prob = 1
    for i,c in enumerate(s):
        methy_ratio = methylation_ratio[i]
        if c=='0': methy_ratio = 1 - methy_ratio
        status_prob *= methy_ratio
    return status_prob

def prob_matrix(methylation_ratio):
    transfer_path = [[0,1,2,4],
                     [0,1,3,5],
                     [0,2,3,6],
                     [1,2,3,7],
                     [0,4,5,6],
                     [1,4,5,7],
                     [2,4,6,7],
                     [3,5,6,7]]
    matrix = np.zeros((8,8))
    for i,nodes in enumerate(transfer_path):
        probs = []
        for node in nodes:
            probs.append(prob_status(methylation_ratio,i)*prob_status(methylation_ratio,node))
        probs = np.array(probs)
        probs = probs/np.sum(probs)
        for j,node in enumerate(nodes):
            matrix[i,node] = probs[j]
    return matrix
        
        
def status_transfer_step_list(start, methylation_ratio=[0.3,0.3,0.2]):
    matrix = prob_matrix(methylation_ratio)
    steps = []
    for end_node in range(8):
        prob = np.zeros(8)
        prob[start] = 1
        step = 0
        count = 1
        while np.sum(prob)>0.05:
            prob = prob.dot(matrix)
            if prob[end_node]!=0:
                step += count * prob[end_node]
                prob[end_node] = 0
            count += 1
        steps.append(step)
    steps = np.array(steps)/np.min(steps)
    return steps



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
    if ratio<0.0001 or ratio>0.9999: return 0
    return -1*(ratio*np.log(ratio)+(1-ratio)*np.log(1-ratio))


def mas(v):
    return (information_entropy(v)+1)**methy_entropy(v)


def conditional_information_entropy(v):
    pass
    # methy_ratio  = target_CpG_methy_ratio(v) + 1e-10
    # entropy = v[0]*np.log(v[0]/methy_ratio)


if __name__=="__main__":
    # testcases = [[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125],
    #              [0,0,0,0,0,0,0,1],
    #              [0.5,0,0,0.5,0,0,0,0]]
    # for testcase in testcases:
    #     print(information_entropy(np.array(testcase)))
    print(status_transfer_prob_list())
