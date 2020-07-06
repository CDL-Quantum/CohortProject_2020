from math import factorial
import numpy as np
from numpy.polynomial.hermite import hermval
from numpy.polynomial.legendre import leggauss

class FCFSpec():
    # convention: hbar = 1.0  
    '''
    This code takes as input:
        - Upper bound of ground and excited vibrational states allowed for 
          transitions
        - Reduced mass of H2 molecule (same as H2+)
        - Frequency of the H2 molecule  (ground state, 0)
        - Frequency of the H2+ molecule (excited state, p)
        - Equilibrium bond lengths (the difference is the displacement)
        - Difference between potential minima of the electronic states 
          (ionization energy)

    This code outputs:
        - FCFs for each pair of ground and excited vibrational states up to a 
          threshold, using the n_H2 = 0, n_H2+ = 0 Franck-Condon Factor as the 
          reference.
        - Spectral positions, SP
    '''
 
    def __init__(self, n_0, n_p):
        self.n_0 = n_0
        self.n_p = n_p
        
        # unit analysis
        self.invcm_to_invEh = 100.0*5.29177210903*pow(10.,-11.)*137.036*2*np.pi
        self.amu_to_me = 1822.888486209
        self.ang_to_bohr = 1.88973     
        
        self.initialize_constants()
        self.initialize_integration_params()

    def initialize_constants(self):
        # hard set for H2 (0) and H2+ (p)
        self.reduced_mass = self.amu_to_me*1.00784/2.
        self.omega_0 = self.invcm_to_invEh*4401.0
        self.omega_p = self.invcm_to_invEh*2322.0
        self.x_eq_0 = self.ang_to_bohr*0.742
        self.x_eq_p = self.ang_to_bohr*1.057
        self.ionization_energy = self.invcm_to_invEh*124418.457

    def initialize_integration_params(self):
        # to calculate overlap integral numerically
        self.R_min = self.ang_to_bohr*(-5.0)
        self.R_max = self.ang_to_bohr*5.0
        self.quadrature_points = 5000

    def H2_energy(self, n):
        return self.omega_0*(n + 0.5) / self.invcm_to_invEh

    def H2p_energy(self, n):
        return self.omega_p*(n + 0.5) / self.invcm_to_invEh + self.ionization_energy

    def H2_psi(self, x, n):
        # n'th harmonic wavefunction for H2  
        tmp = self.reduced_mass*self.omega_0
        prefactor = pow(tmp / np.pi, 0.25) / np.sqrt(2.0**n * factorial(n))
        coeffs = np.zeros(n+1, float)
        coeffs[n] = 1.0
        hermpoly = hermval(np.sqrt(tmp)*(x - self.x_eq_0), coeffs)
        return prefactor*np.exp(-0.5*tmp*(x - self.x_eq_0)**2)*hermpoly
    
    def H2p_psi(self, x, n):
        # n'th harmonic wavefunction for H2+  
        tmp = self.reduced_mass*self.omega_p
        prefactor = pow(tmp / np.pi, 0.25) / np.sqrt(2.0**n * factorial(n))
        coeffs = np.zeros(n+1, float)
        coeffs[n] = 1.0
        hermpoly = hermval(np.sqrt(tmp)*(x - self.x_eq_p), coeffs)
        return prefactor*np.exp(-0.5*tmp*(x - self.x_eq_p)**2)*hermpoly

    def spectrum_analysis(self):
        # calculate FCFs
        GLquad = leggauss(self.quadrature_points)

        all_data = []

        for k in range(self.n_0+1):
            E0 = self.H2_energy(k)
            for l in range(self.n_p+1):
                Ep = self.H2p_energy(l)

                overlap = 0
                for p in range(self.quadrature_points):
                    
                    new_point = 0.5*(self.R_max - self.R_min)*GLquad[0][p] + 0.5*(self.R_max + self.R_min)
                    new_weight = 0.5*(self.R_max - self.R_min)*GLquad[1][p]
                    overlap += self.H2_psi(new_point, k) * self.H2p_psi(new_point, l) * new_weight

                FCF = overlap**2.
        
                # n_0   n_p   FCF
                data = np.zeros(3)

                if (k==0 and l==0):
                    reference = FCF
                
                FCF /= reference

                data[0] = k
                data[1] = l
                data[2] = FCF

                all_data.append(data)

        return np.array(all_data)
