## my first function ##

## import stuff ##
import numpy as np

class GCPRE:
    """
    
    """

    def GCPRE_profile(idx, label, res_start, res_end):
        ## set some variables; these are largely unchanged by user ##
        Kcm = 1.23e-32    # cm^6 / s^2
        K = Kcm * 1e42    # nm^6 / s^2

        l = 0.38    # monomer length, in nm

        alpha = 0.8    # cos of bond angle supplements for freely rotating chain model

        tau_c = 4e-9   # effective correlation time, in seconds

        W_H = 600e6     # Larmor frequency of proton, in Hz

        R2D = float(14)    # transverse relaxation rate, in 1/sec

        t = 0.012   # INEPT delay time, in sec

        ## denote the residue numbers & generate array of inter-residue monomer numbers (nVal) ##
        #res_start = 268 
        #res_end = 414
        diff = res_end - res_start

        #label = 368

        #res_list = list(range(res_start, res_end))

        nVal = []

        for res in range(res_start, res_end):
            nVal.append(float(abs(res - label)))


        ## calculate <r^2> between residues ##
        r = []

        for n in nVal:
            num = (2 * alpha * (1 - (alpha ** n))) 
            denom = ((n * (1 - alpha) ** 2))
            temp = np.sqrt(n * l * l * (((1 + alpha) / (1 - alpha)) - np.divide(num, denom)))
            r.append(temp)

        #print(r)
        r6 = np.power(r, 6)

        ## calculate PRE contribution to transverse relaxation rate (gamma) ##
        paren = K * ((4 * tau_c) + (3 * tau_c / (1 + (W_H * W_H * tau_c * tau_c))))
        gamma_GC = np.divide(paren, r6)

        ## calculate intensity ratios ##

        top = R2D * np.exp(-gamma_GC * t)
        bottom = R2D + gamma_GC

        ratio = np.divide(top, bottom)

        ## NB - the labeled residue has a distance of n = 0 monomers from the label position
            ## although this means that the r distance is div/0, we bypass this in the previous cells and then set it to zero below (for plotting purposes)
        ratio[np.isnan(ratio)] = 0

        return ratio