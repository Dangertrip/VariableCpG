import numpy as np
from utils import *

def distance_entropy(sample_num=10000, certain_point=[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]):
    result = []
    for i in range(sample_num):
        v1 = vector_generator()
        result.append([information_entropy(v1),linear_distance(v1, certain_point)])
    result = np.array(result)
    return result[:,0],result[:,1]

if __name__=="__main__":
    print(information_entropy(np.array([0,0,0,0,0,0,0,1])))
    print(distance_entropy().shape)
