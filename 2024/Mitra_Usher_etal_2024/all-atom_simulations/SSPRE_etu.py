## import stuff ##

import mdtraj as md
import numpy as np
import scipy


## my makeshift PRE simulator ##

# define variables #

original_K = 1.2300e-32     # in cm^6 / s^2
Kcon = original_K * 1e42    # in nm^6 / s^2

class SSPRE:
    """
    
    """
    
    def __init__(self, ProtObj, tau_c, t_delay, R_2D, W_H):
    
        """
        tau_c - effective correlation time, in ns
        t_delay - INEPT delay time, in ms
        R_2D - transverse relaxation rate of backbone amide protons in diamagnetic protein form, in Hz
        W_H - Larmor frequency for proton, in Hz
    
        """
    
        self.SSPO = ProtObj
    
        self.t_delay = float(t_delay)   # in ms
        self.R_2D = float(R_2D)    # in Hz
        self.tau_c = tau_c    # in ns
        self.W_H = W_H   # in Hz
    
        # perform some conversions #
        tau_c = float(tau_c) / 10e9    # now in seconds
        tau_c_sq = tau_c * tau_c
    
        # compute 'prefactor' term used later
        W_H_sq = W_H * W_H
        PREFACTOR = ((3 * tau_c) / (1 + W_H_sq * tau_c_sq)) + (4 * tau_c)
        self.PREFACTOR = PREFACTOR * Kcon
        
    def __repr__(self):
        """

        """
        return "["+hex(id(self))+ "]: SSPRE OBJ - (R_2D = %3.2f Hz, t_delay = %3.2f ms, tau_c = %3.2f ns, 1H Larmor = %3.3e Hz)" % (self.R_2D, self.t_delay, self.tau_c, self.W_H)

    ## makes PRE and gamma profiles ##

    def generate_PRE_profile(self, label_position, spin_label_atom = 'CB', target_relaxation_atom = 'N'):
        """

        label_position - int, position in sequence where spin label is located
        spin_label_atom - str (default is CB)
        target_relaxation_atom - str (default is N)

        >>returns a 2-place tuple (pos 0 is PRE intensity profile, pos 1 is PRE 1H relaxation profile

        """

        # first get index value of all residues
        residue_list = sorted(list(self.SSPO._SSProtein__CA_residue_atom.keys()))

        # calculate the mean rij distance for residue pairs

        gamma = []

        for idx in residue_list:
            r_6_nm = np.power(0.1 * self.SSPO.get_inter_residue_atomic_distance(label_position, idx, A1 = spin_label_atom, A2 = target_relaxation_atom), 6)
            gamma.append(np.mean(self.PREFACTOR / r_6_nm))


        # convert delay to seconds
        t_delay_sec = self.t_delay / 1000

        # for each gamma get intensity ratio
        profile = []
        for g in gamma:
            profile.append((self.R_2D * np.exp(-g * t_delay_sec)) / (self.R_2D + g))

        return (profile, gamma)