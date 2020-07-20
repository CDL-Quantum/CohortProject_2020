import abc

import numpy as np

class AbstractIsing(abc.ABC):
    @abc.abstractmethod
    def energy(self):
        """Returns the energy of the current spin configuration"""

    @abc.abstractmethod
    def energy_diff(self, *coords):
        """Returns the energy difference resulting from flipping the site at the given coordinates"""
        
    @abc.abstractmethod
    def rand_site(self):
        """Selects a site in the lattice at random"""
    
    def mc_step(self, T):
        """Performs a full update of the given Ising model using the Metropolis-Hastings algorithm"""
        current_energy = self.energy()
        for _ in range(self.num_spins):
            site = self.rand_site()
            dE = self.energy_diff(*site)

            if (dE < 0) or (np.random.rand() < np.exp(-dE / T)):
                current_energy += dE
                self.spins[site] *= -1

        return current_energy