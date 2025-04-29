import numpy as np
import math
from interfaces import IParticle

class photon(IParticle):
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
        self._dx, self._dy, self._dz = self._calculate_direction()

    def batch_calculate_direction(pz, theta, phi):
        # im pretty certain we avoid division by zero if we avoid tan(phi) entirely, dz=0 should have a valid direction vector right?
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        dx = cos_theta * cos_phi 
        dy = sin_theta * cos_phi
        dz = sin_phi
        
        return np.column_stack([dx, dy, dz])

    def get_origin(self):
        return (self.x, self.y, self.z)

    def get_direction(self):
        return (self._dx, self._dy, self._dz)

    def energy(self):
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)