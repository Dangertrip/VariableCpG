import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def scatter_plot(x,y,args={}):
    plt.figure()
    color = args.get('color', 'red')
    alpha = args.get('alpha', 0.3)
    plt.scatter(x, y, color = color, alpha = alpha)
    filename = args.get('figure_name','default.'+str(time.time())+'.png') 
    plt.savefig(filename)

def scatter_plot_3d(x,y,z,args={}):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    color = args.get('color', 'red')
    alpha = args.get('alpha', 0.3)
    for xx,yy,zz in zip(x,y,z):
        ax.scatter(xx,yy,zz,alpha = alpha, color = color)
    ax.set_xlabel('Methy_ratio1')
    ax.set_ylabel('Methy_ratio2')
    ax.set_zlabel('dist')
    filename = args.get('figure_name','default.'+str(time.time())+'.png') 
    plt.savefig(filename)
