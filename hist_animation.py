# Fixed bin width
import gzip, pickle, sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  
from matplotlib.pylab import *
import matplotlib.animation as animation

def main():
    fname = sys.argv[1] 
    fh = gzip.open(fname, "rb")
    params = pickle.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T, V, H, I = [], [], [], []
    for epoch in range(params["epochs"]):
        viruses = pickle.load(fh)
        hosts = pickle.load(fh)
        interactions = pickle.load(fh)
        T.append(epoch)
        V.append(viruses)
        H.append(hosts)
        I.append(interactions)
 
    fh.close()

    bin_width_V= float(input('input bin size! : '))   
    bin_width_mumax= float(input('input bin size! : '))
    bin_width_mass= float(input('input bin size! : ')) 

    # create figure object
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(8, 6))

    def updateHist(frame):           
        t= frame 

        V_beta, H_muMax, H_mass = [], [], []        
        for v in range(len(V[t])):        
            V_beta.append(V[t][v].beta)

        for h in range(len(H[t])):
            H_muMax.append(H[t][h].mu_max)
            H_mass.append(H[t][h].mass)
        
        ep_V, ep_Hm, ep_HM= 0, 0, 0

        for ax in (ax0, ax1, ax2):
            ax.clear()

        # subplot(3, 1, 1) 
        if len(V_beta)==0:
            ax0.plot([])
        else:
            if min(V_beta)== max(V_beta):
                ep_V= bin_width_V
            bins_V = np.arange(min(V_beta), max(V_beta) + ep_V + bin_width_V, bin_width_V)         
            ax0.hist(V_beta, bins_V, color='red', alpha=0.5)
        ax0.set_ylabel("virus, beta")
        ax0.legend([f'Epoch {t}']) 
        ax0.set_title('bin width- beta= {}, muMax= {}, mass= {}'.format(bin_width_V, bin_width_mumax, bin_width_mass) )
        
        # subplot(3, 1, 2)
        if len(H_muMax)==0: 
            ax1.plot([])
        else:
            if min(H_muMax)== max(H_muMax):
                ep_HM = bin_width_mumax
            bins_H_mumax = np.arange(min(H_muMax), max(H_muMax) + ep_HM + bin_width_mumax, bin_width_mumax)                  
            ax1.hist(H_muMax, bins_H_mumax, color='blue', alpha=0.5)
        ax1.set_ylabel("Host mu_max")
        ax1.legend([f'Epoch {t}'])
        
        # subplot(3, 1, 3)   
        if len(H_mass)==0:
            ax2.plot([]) 
        else:
            if min(H_mass)== max(H_mass):
                ep_Hm= bin_width_mass        
            bins_H_mass = np.arange(min(H_mass), max(H_mass) + ep_Hm + bin_width_mass, bin_width_mass)
            ax2.hist(H_mass, bins_H_mass, color='green', alpha=0.5)
        ax2.set_ylabel("Host Mass")
        ax2.legend([f'Epoch {t}'])

        print(" Epoch : %d" %t)
        return fig,

    simulation = animation.FuncAnimation(fig, updateHist, frames=240, blit=True)
    simulation.save("%s" % (fname.replace("pkl", "gif")), writer='PillowWriter', fps=1)


if __name__ == "__main__":
    main()
