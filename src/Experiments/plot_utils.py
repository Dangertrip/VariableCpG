import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import time

def scatter_plot(x,y,args={}):
    plt.figure()
    color = args.get('color', 'red')
    alpha = args.get('alpha', 0.3)
    plt.scatter(x, y, color = color, alpha = alpha)
    filename = args.get('figure_name','default.'+str(time.time())+'.png') 
    plt.savefig(filename)
