from distance_entropy import distance_entropy 
from distance_methyratio import distance_methyratio_random_vs_fixed
from methyratio_distance_times_entropy import methyratio_distance_times_entropy
from plot_utils import scatter_plot, scatter_plot_3d
from methyratio_entropy import methyratio_entropy
from methyratio_infoentropy_times_methentropy import methyratio_infoentropy_times_methentropy
from methyratio_dist_times_mas import methyratio_dist_times_mas

if __name__=='__main__':
    # distance(certain point, random points) vs information entropy
    # x,y = distance_entropy(100000)
    # scatter_plot(x, y, {'figure_name':'distance_vs_entropy.png'})

    # methylation(centain point's methy ratio=0) vs distance(certain point, random point)
    # x,y = distance_methyratio_random_vs_fixed(100000)
    # scatter_plot(x, y, {'figure_name':'distance_vs_methyratio_fixed_point.png'})

    # x,y = methyratio_distance_times_entropy(200000)
    # scatter_plot(x, y, {'figure_name':'methyratiodiff_vs_distance_times_entropy.png'})

    # x,y = methyratio_entropy(200000)
    # scatter_plot(x, y, {'figure_name':'methyratio_vs_methyentropy.png'})

    x,y = methyratio_infoentropy_times_methentropy(500000)
    scatter_plot(x, y, {'figure_name':'methyratio_infoentropy_times_methentropy.png'})
    
    # x,y,z = methyratio_dist_times_mas(50000)
    # scatter_plot_3d(x, y, z, {'figure_name':'methyratio_dist.png'})

