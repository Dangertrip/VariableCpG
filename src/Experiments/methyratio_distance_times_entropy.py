import numpy as np
from utils import vector_generator, target_CpG_methy_ratio, linear_distance, information_entropy


def methyratio_distance_times_entropy(sample_size=10000):
    result = []
    for i in range(sample_size):
        v1 = vector_generator()
        v2 = vector_generator()
        methyratio_diff = target_CpG_methy_ratio(v1)-target_CpG_methy_ratio(v2)
        methy_v1 = target_CpG_methy_ratio(v1)
        methy_v2 = target_CpG_methy_ratio(v2)
        entropy_v1 = -1*(methy_v1*np.log(methy_v1)+(1-methy_v1)*np.log(1-methy_v1))
        entropy_v2 = -1*(methy_v2*np.log(methy_v2)+(1-methy_v2)*np.log(1-methy_v2))
        neg_methyratio_diff = 1 - methyratio_diff
        entropy = -1*(methyratio_diff*np.log(methyratio_diff, where=methyratio_diff>0)+neg_methyratio_diff*np.log(neg_methyratio_diff,where=neg_methyratio_diff>0))
        result.append([methyratio_diff,
                       linear_distance(v1*(information_entropy(v1)*entropy_v1+1), v2*(information_entropy(v2)*entropy_v2+1))])
    result = np.array(result)
    return result[:, 0], result[:, 1]
