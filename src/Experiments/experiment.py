from distance_entropy import distance_entropy 
from distance_methyratio import distance_methyratio_random_vs_fixed
from plot_utils import scatter_plot 

if __name__=='__main__':
    # distance(certain point, random points) vs information entropy
    # x,y = distance_entropy(100000)
    # scatter_plot(x, y, {'figure_name':'distance_vs_entropy.png'})

    # methylation(centain point's methy ratio=0) vs distance(certain point, random point)
    x,y = distance_methyratio_random_vs_fixed(100000)
    scatter_plot(x, y, {'figure_name':'distance_vs_methyratio_fixed_point.png'})


