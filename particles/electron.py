import numpy as np
import math
from interfaces import IParticle
import random as rd
m_e = 511
class electron(IParticle):
    def __init__(self, px, py, pz, theta, phi, x, y, z, t):
        self.px = px
        self.py = py
        self.pz = pz
        self.theta = theta
        self.phi = phi
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.m = m_e
        self._dx, self._dy, self._dz = self._calculate_direction()
        
    def get_energy(self, fano):
        true_energy = (self.px**2 + self.py**2 + self.pz**2)/(2*self.m)
        sigma = math.sqrt(true_energy) * fano #The true standard deviation for a poisson (or gaussian) dist is the square root of the mean. However, creating charge carriers/exciting atoms in a scintillator crystal has less variation than would be predicted by a poisson distribution. Fortunately, the deviation from this exected standard deviation is defined by the fano factor-an inherent quality for each detector.
        return rd.gauss(true_energy, sigma)
    def _calculate_direction(self):
        if self.pz < 0:
            dz = -1
            dx = -np.cos(self.theta) * np.cos(self.phi) / np.sin(self.phi)
            dy = -np.sin(self.theta) * np.cos(self.phi) / np.sin(self.phi)
        else:
            dz = 1
            dx = np.cos(self.theta) * np.cos(self.phi) / np.sin(self.phi)
            dy = np.sin(self.theta) * np.cos(self.phi) / np.sin(self.phi)
        return dx, dy, dz