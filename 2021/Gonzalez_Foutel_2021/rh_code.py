import numpy as np

from afrc import AnalyticalFRC

def compute_rg_re_rh(seq, count=2000):
        """
        Returns the apparent hydrodynamic radius as calculated based on the approximation
        derived by Nygaard et al. [1]. Returns a hydrodynamic radius in Angstroms.

        Parameters (alpha1/2/3 should not be altered to recapitulate behaviour defined
        by Nygaard et al.

        NOTE that this approximation holds for fully flexible disordered proteins, but will likely
        become increasingly unreasonable in the context of larger and larger folded domains.

        References:
        -------------
        [1] Nygaard M, Kragelund BB, Papaleo E, Lindorff-Larsen K. An Efficient
        Method for Estimating the Hydrodynamic Radius of Disordered Protein
        Conformations. Biophys J. 2017;113: 550–557.

        Radius of gyration is returned in Angstroms.


        Parameters
        ---------------
        rg : float
            The radius of gyration, passed in Angstroms

        alpha1 : float {0.216}
           First parameter in equation (7) from Nygaard et al.

        alpha2 : float {4.06}
           Second parameter in equation (7) from Nygaard et al.

        alpha3 : float {0.821}
           Third parameter in equation (7) from Nygaard et al.

        Returns
        -----------
        np.ndarray
            Returns a numpy array with per-frame instantaneous hydrodynamic radii


        """

        alpha1 = 0.216
        alpha2 = 4.06
        alpha3 = 0.821

        A = AnalyticalFRC(seq)
        nres = len(seq)

        rg = A.sample_rg_distribution(count)
        re = A.sample_re_distribution(count)

        
        


        # precompute
        N_033 = np.power(nres, 0.33)
        N_060 = np.power(nres, 0.60)

        Rg_over_Rh = ((alpha1*(rg - alpha2*N_033)) / (N_060 - N_033)) + alpha3

        rh =  (1/Rg_over_Rh)*rg
        return [rg, re, rh]
