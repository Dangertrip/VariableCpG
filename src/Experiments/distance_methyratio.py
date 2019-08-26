import numpy as np
from utils import target_CpG_methy_ratio, linear_distance

def distance_methyratio_random_vs_fixed(sample_size=10000, fixed_point=[0,0,0,0,0,0,0,1]):
    result = []
    for i in range(sample_size):
        v = vector_generator()
        result.append([target_CpG_methy_ratio(v),linear_distance(v,fixed_point)])
    result = np.array(result)
    return result[:,0], result[:,1]

def distance_methyratio_random_vs_random():
    pass
