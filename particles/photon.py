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

    def get_origin(self):
        return (self.x, self.y, self.z)

    def get_direction(self):
        return (self._dx, self._dy, self._dz)

    def energy(self):
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)