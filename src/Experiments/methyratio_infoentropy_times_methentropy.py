import numpy as np
from utils import vector_generator, target_CpG_methy_ratio, linear_distance, information_entropy


def methyratio_infoentropy_times_methentropy(sample_size=10000):
    result = []
    for i in range(sample_size):
        v1 = vector_generator()
        # v2 = vector_generator()
        # methyratio_diff = target_CpG_methy_ratio(v1)-target_CpG_methy_ratio(v2)
        methy_v1 = target_CpG_methy_ratio(v1)
        # methy_v2 = target_CpG_methy_ratio(v2)
        entropy_v1 = -1*(methy_v1*np.log(methy_v1) +
                         (1-methy_v1)*np.log(1-methy_v1))
        mass = np.array([0.125, 0.125, 0.125, 0.125,
                         0.125, 0.125, 0.125, 0.125])
        entropy = (methy_v1/np.abs(methy_v1))*np.abs(methy_v1-0.5)**((information_entropy(v1)+1))
        # ((information_entropy(v1)+1)**(entropy_v1+1))**(1+methy_v1)
        # entropy_v2 = -1*(methy_v2*np.log(methy_v2)+(1-methy_v2)*np.log(1-methy_v2))
        result.append([methy_v1, entropy])
    result = np.array(result)
    return result[:, 0], result[:, 1]
