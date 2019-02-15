import matplotlib.pyplot as plt
import numpy as np
import matplotlib
data = np.random.randint(0,6,9216).reshape((96,96))


cmap = matplotlib.colors.ListedColormap(['white', 'forestgreen', 'blue', 'purple','coral', 'red','black'])
norm = matplotlib.colors.BoundaryNorm([0,1,2,3,4,5,6,7],cmap.N)

data_R = np.load('./PM-R.npy')
data_F = np.load('./PM-F.npy')

 
plt.figure(dpi=600)
def plot_heatmap(data):
    im = plt.imshow(data, cmap=cmap, norm=norm)
    ax = plt.gca()

    cbar = ax.figure.colorbar(mappable=im, ax=ax)
    cbar.ax.set_ylabel('No. of BESs', rotation=-90, va='bottom')
    cbar.ax.get_yaxis().set_ticks([])
    for j, lab in enumerate(['$0$','$1$','$2$','$3$','$4$','$5$','$6$']):
        cbar.ax.text(.5, (2 * j + 1) / 14.0, lab, ha='center', va='center')


    #cbar.ax.set_xticks[np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])]
    #cbar.ax.set_xticklabels(np.array(['0','1','2','3','4','5']))

    ax.set_xticks(np.array([0,8,16,24,32,40,48,56,64,72,80,88]))
    ax.set_xticklabels(np.array(['Y01','Y09','Y17','Y25','Y33','Y41','Y49','Y57','Y65','Y73','Y81','Y89']),fontdict={'fontsize':6})

    ax.set_yticks(np.array([0,8,16,24,32,40,48,56,64,72,80,88]))
    ax.set_yticklabels(np.array(['X01','X09','X17','X25','X33','X41','X49','X57','X65','X73','X81','X89']),fontdict={'fontsize':6})

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)


plot_heatmap(data_F)
plt.savefig('BES_F.png',dpi=600)
plt.show()
