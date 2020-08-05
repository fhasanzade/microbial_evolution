
import numpy, gzip, pickle, dill, sys
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.pylab import *
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    fname = sys.argv[1]
    
    fh = gzip.open(fname, "rb")
    params = dill.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T, DIP, H, V, I = [], [], [], [], []
    max_I=0
    for epoch in range(params["epochs"]):
        dip = dill.load(fh)
        hosts = dill.load(fh)
        viruses = dill.load(fh)
        interactions = dill.load(fh)
        T.append(epoch)
        DIP.append(dip)
        V.append(viruses)
        H.append(hosts)
        I.append(interactions)

        #Calculate the Max number of Interaction for fixing the Color Bar in Animation results
        if len(interactions) > max_I:
            max_I= len(interactions)

    fh.close()

    # bin widths for interaction (based on the bins for host)
    bin_width_I_H= float(input('input bin size! (default: 0.001): '))
    bin_width_I_V= float(input('input bin size! (default: 0.001): '))

    # create figure object
    fig, (ax0) = plt.subplots(1,1, figsize=(8, 6))
    ax0=plt.gca()
    # for updating the colorbar in each iteration
    ax0_divider = make_axes_locatable(ax0)
    cax0 = ax0_divider.append_axes("right", size="2%", pad="2%")

    #def updateHist(frame):
    def updateHist(frame):

        t=frame        
        print(" Epoch : %d" %t)
        
        I_virus, I_host = [], []

        ep_I_V, ep_I_H= 0, 0
        
        for i in range(len(I[t])):        
            #Interaction.append(I[t][i])
            I_virus.append(I[t][i][0].genotype)
            I_host.append(I[t][i][1].genotype)
        
        # Plotting, clear axis  
        ax0.clear()    
        cax0.cla()
                
        if len(I_virus)==0:    # in case of no interaction 
            ax0.plot([])
        else: 
            if min(I_virus)== max(I_virus):
                ep_I_V= min(I_virus)
                
            if min(I_host)== max(I_host):
                ep_I_H= min(I_host)
                
            # hist Plotting
            bins_I = [np.arange(min(I_host), max(I_host)+ ep_I_H + bin_width_I_H, bin_width_I_H), np.arange(min(I_virus), max(I_virus)+ep_I_V + bin_width_I_V, bin_width_I_V)]
            
            norm = matplotlib.colors.Normalize(vmin=0, vmax=max_I)   #fixed color bar for all iterations
            f0=ax0.hist2d(I_host, I_virus, bins_I, norm=norm, cmap='jet', edgecolor='')
            plt.colorbar(f0[3], cax=cax0)                    

        ax0.set_ylabel("virus, genotype (beta)")    
        ax0.set_xlabel("host, genotype (mu_max)")    
        ax0.set_title('Epoch= {}, bin width- I_V= {}, I_H = {}'.format(t, bin_width_I_V, bin_width_I_H) )
        ax0.set_xlim([0,0.2])
        ax0.set_ylim([0,0.2]) 
        
        return fig, 

    simulation = animation.FuncAnimation(fig, updateHist, frames=300, blit=True)
    simulation.save('interaction_Plot_hist2D.gif', writer='PillowWriter', fps=1)


if __name__ == "__main__":
    main()