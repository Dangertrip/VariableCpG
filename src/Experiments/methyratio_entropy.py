import numpy as np
from utils import vector_generator, target_CpG_methy_ratio, linear_distance, information_entropy


def methyratio_entropy(sample_size=10000):
    result = []
    for i in range(sample_size):
        v1 = vector_generator()
        v2 = vector_generator()
        methy_v1 = target_CpG_methy_ratio(v1)
        entropy = -1 * (methy_v1*np.log(methy_v1)+(1-methy_v1)*np.log(1-methy_v1))
        result.append([methy_v1, entropy])
    result = np.array(result)
    return result[:, 0], result[:, 1]
