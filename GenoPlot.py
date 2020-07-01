import numpy, gzip, pickle, dill
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import pandas as pd
import numpy as np  

def main():
    fname = sys.argv[1]

    fh = gzip.open(fname, "rb")
    params = dill.load(fh)
    for param in params.keys():
        print(f"{param} : {params[param]}")
    T, DIP, H, V, I = [], [], [], [], []
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
    fh.close()

    # using bins for each time step
    # bin widths
    bin_width_V= float(input('input bin size! (default: 0.01): '))   
    bin_width_genotype= float(input('input bin size! (default: 0.01): '))
    bin_width_mass= float(input('input bin size! (default: 0.1): ')) 

    #for fixed axes plots and calculating bin width
    V_beta_all, H_mass_all, H_genotype_all =[], [], []  #y axis for a compelete plot 
    X_V, X_H = [], []   #x axis for a compelete plot 
    Vn_beta, Hn_mass, Hn_genotype = [], [], []
    X_Vn, X_Hn, X_Hng = [], [], []  # x axis for bin plot  
    bins_V_all, bins_H_genotype_all, bins_H_mass_all =[], [], []   # y axis for bin plot
            
    #def updateHist(frame):
    for t in range(len(T)):
        
        V_beta, H_genotype, H_mass = [], [], []
        bins_V, bins_H_genotype, bins_H_mass =[], [], []
        ep_V, ep_Hm, ep_HM= 0, 0, 0
        
        for v in range(len(V[t])):        
            V_beta.append(V[t][v].genotype)

        if len(V_beta)> 0:    
            V_beta_all.extend(V_beta)
            X_V.extend(t * np.ones(len(V[t])))

            # actual bins valuse + same size x axis for time 
            # for time t
            bins_V= np.arange(min(V_beta), max(V_beta) + ep_V + bin_width_V, bin_width_V)
            nV, _ = np.histogram(V_beta, bins=bins_V)  #abundance in each bin
            # for all times
            bins_V_all.extend(bins_V[0:-1])
            X_Vn.extend(t * np.ones(len(nV)))
            #all abundance in all times
            Vn_beta.extend(nV)      

        for h in range(len(H[t])):
            H_genotype.append(H[t][h].genotype)
            H_mass.append(H[t][h].mass)

        if len(H_genotype)> 0:    
            H_mass_all.extend(H_mass)
            H_genotype_all.extend(H_genotype) 
            X_H.extend(t * np.ones(len(H[t])))

            # actual bins valuse + same size x axis for time 
            # for time t
            bins_H_genotype= np.arange(min(H_genotype), max(H_genotype) + ep_HM + bin_width_genotype, bin_width_genotype)
            bins_H_mass= np.arange(min(H_mass), max(H_mass) + ep_Hm + bin_width_mass, bin_width_mass)
            nHgen, _ = np.histogram(H_genotype, bins=bins_H_genotype)  #abundance in each bin
            nHmass, _ = np.histogram(H_mass, bins=bins_H_mass)  #abundance in each bin
            # for all times
            bins_H_genotype_all.extend(bins_H_genotype[0:-1])
            bins_H_mass_all.extend(bins_H_mass[0:-1])
            X_Hng.extend(t * np.ones(len(nHgen)))
            X_Hn.extend(t * np.ones(len(nHmass)))
        
            #all abundance in all times 
            Hn_mass.extend(nHmass)
            Hn_genotype.extend(nHgen)
        
        print(" Epoch : %d" %t)


    #now scatter plots using bins
    # X asix: time, Y axis: actual values for genotypes/mass, colorbar: abundance

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(8, 6))
    #colormap = plt.cm.bwr

    # subplot(3, 1, 1)
    f0=ax0.scatter(X_Vn, bins_V_all, c=Vn_beta , s=5, cmap='jet', edgecolor='', marker='.')
    ax0.set_ylabel("virus, genotype (beta)")    
    ax0.set_title('bin width- V.genotype= {}, H.genotype= {}, H.mass= {}'.format(bin_width_V, bin_width_genotype, bin_width_mass) )
    fig.colorbar(f0, ax=ax0)

    # subplot(3, 1, 2)  
    f1=ax1.scatter(X_Hng, bins_H_genotype_all,  c=Hn_genotype, s=5, cmap='jet', edgecolor='', marker='.')
    ax1.set_ylabel("Host genotype")
    fig.colorbar(f1, ax=ax1)

    # subplot(3, 1, 3)  
    f2= ax2.scatter(X_Hn, bins_H_mass_all,  c=Hn_mass, s=5, cmap='jet', edgecolor='', marker='.')
    fig.colorbar(f2, ax=ax2)
    ax2.set_ylabel("Host Mass")
    ax2.set_xlabel("time")

    # save as pdf
    fig.savefig("GenoPlot.pdf", bbox_inches='tight')
    #plt.savefig("%s" % (fname.replace("pkl", "pdf")), format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()