## ------
## Mathieu Bourdenx 2021
## -----


import numpy as np
import matplotlib.pyplot as plt

## Graphical setup
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

## Helper functions
def lambdaf(nu, D, k):
    return (np.sqrt(nu**2 + 4*D*k) - nu) / (2*D)


## RNA distribution
def r_ss_den(x, nu, D, k, betaR):
    ''' Steady state distribution of RNA in dendrites (Equation 1)'''
    lambdaR = lambdaf(nu, D, k)

    return (betaR * lambdaR / k) * np.exp(-lambdaR * x)


## Protein distribution
def p_ss_den(x, nuR, DR, kR, betaR, nuP, DP, kP, betaP):
    ''' Steady state distribution of dendritic protein dynamic (Equation 2)'''
    lambdaR = lambdaf(nuR, DR, kR)
    lambdaP = lambdaf(nuP, DP, kP)
    
    p_den = ((betaP * betaR * lambdaR) / (kR * (DP * lambdaR**2 + nuP * lambdaR - kP)))*(-np.exp(-lambdaR * x) + ((DP * lambdaP * lambdaR + nuP * lambdaP)/kP) * np.exp(-lambdaP * x))

    return p_den

def p_ss_som(x, nuP, DP, kP, betaP):
    ''' Steady state distribution of somatically 
    localized proteins (Equation 3)'''

    lambdaP = lambdaf(nuP, DP, kP)
    
    p_som = (betaP * betaR * lambdaP) / (kR * kP) * np.exp(-lambdaP * x)

    return p_som

def p_ss_total(x, nuR, DR, kR, betaR, nuP, DP, kP, betaP):
    ''' Mixture of somatic and dendritic protein distribution '''

    # Fraction of somatic mRNA
    s = 0.557
    # Fraction of somatically generated proteins in the dendrites
    z = 0.2231802348794625
    # Translation efficiency ratio
    Tr = 1.217066737341612
    
    return s * z * p_ss_som(x, nuP, DP, kP, betaP) + (1-s) * p_ss_den(x, nuR,DR, kR, betaR, nuP, DP, kP, betaP) * Tr

    
if __name__ == '__main__':

    ### Constant for RNA ###
    # Velocity
    nuR_high = 1.3
    nuR_low = 5.8e-3

    # Diffusion constant
    DR_high = 3.8e-3
    DR_low = 3e-3

    # Transcription rate
    betaR = 0.001

    # Degradation rate
    kR = np.log(2) / (16 * 3600)


    x = np.linspace(0, 100, num=100)

    # Computation for upper and lower bounds 
    rna_low = r_ss_den(x, nuR_low, DR_low, kR, betaR)
    rna_high = r_ss_den(x, nuR_high, DR_high, kR, betaR)

    # Normalization to somatic values (x=0)
    rna_low_normalized = rna_low / rna_low[0]
    rna_high_normalized = rna_high / rna_high[0]

    # Plotting
    plt.figure(figsize=(6, 4), constrained_layout=True)
    plt.plot(x, rna_low_normalized, c='lightblue')
    plt.plot(x, rna_high_normalized, c='darkblue')
    plt.xlim(0, 100)
    plt.xticks([0, 50, 100])
    plt.xlabel('Dendritic distance (µm)', fontdict={'size':14})
    plt.ylim(0, 1.1)
    plt.yticks([0, 0.5 , 1])
    plt.ylabel('mRNA density', fontdict={'size': 14})
    plt.text(x=101, y=rna_low_normalized[-1],
             s=np.round(rna_low_normalized[-1], decimals=3),
             fontdict={'size': 14})
    plt.text(x=101, y=rna_high_normalized[-1],
             s=np.round(rna_high_normalized[-1], decimals=3),
             fontdict={'size': 14})

    plt.savefig('../Plots/Fig1A.png', dpi=150)
    plt.show()


    ## -------------------------------------------------------------------------
    ## PROTEIN DISTRIBUTION
    ## -------------------------------------------------------------------------
    
    ## Constant for protein
    # Velocity
    nuP = 0

    #Diffusion
    DP_high = 4.5
    DP_low = 0.023

    #Translation rate
    #betaP = 0.021
    betaP = 0.02301
    
    # Degradation rate
    kP = np.log(2) / (6.64 * 24 * 3600)



    # Computation for upper (high) and lower (low) bounds
    protein_low = p_ss_total(x,
                             nuR_low, DR_low, kR, betaR,
                             nuP, DP_low, kP, betaP)
    
    protein_high = p_ss_total(x,
                              nuR_low, DR_low, kR, betaR,
                              nuP,DP_high, kP, betaP)
    
    #Normalization to somatic values (x=0)
    protein_low_normalized = protein_low / protein_low[0]
    protein_high_normalized = protein_high / protein_high[0]


    # Plotting
    plt.figure(figsize=(6, 4), constrained_layout=True)
    plt.plot(x, protein_high_normalized, c='red')
    plt.plot(x, protein_low_normalized, c='orange')
    plt.xlim(0, 100)
    plt.xticks([0, 50, 100])
    plt.xlabel('Dendritic distance (µm)', fontdict={'size': 14})
    plt.ylim(0, 1.2)
    plt.yticks([0, 0.5, 1])
    plt.ylabel('Protein density', fontdict={'size': 14})

    plt.text(x=101, y=protein_low_normalized[-1],
             s=np.round(protein_low_normalized[-1], decimals=3),
             fontdict={'size': 14})
    
    plt.text(x=101, y=protein_high_normalized[-1],
             s=np.round(protein_high_normalized[-1], decimals=3),
             fontdict={'size': 14})

    plt.savefig('../Plots/Fig1b.png', dpi=300)

    plt.show()  
