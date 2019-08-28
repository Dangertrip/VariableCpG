import numpy as np
from utils import vector_generator, target_CpG_methy_ratio, linear_distance, information_entropy, mas


def methyratio_dist_times_mas(sample_size=10000):
    result = []
    for i in range(sample_size):
        v1 = vector_generator()
        v2 = vector_generator()
        methy_v1 = target_CpG_methy_ratio(v1)
        methy_v2 = target_CpG_methy_ratio(v2)
        dist = linear_distance(v1,v2)
        # entropy_v1 = -1*(methy_v1*np.log(methy_v1)+(1-methy_v1)*np.log(1-methy_v1))
        # entropy_v2 = -1*(methy_v2*np.log(methy_v2)+(1-methy_v2)*np.log(1-methy_v2))
        result.append([methy_v1 - methy_v2, dist/mas(v1)-mas(v2)])
    result = np.array(result)
    return result[:, 0], result[:, 1], result[:, 2]
