import numpy as np
import math
from interfaces import IParticle
from particles.electron import electron
m_e = 511

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
        self._dx, self._dy, self._dz = self.get_direction()

    def batch_calculate_direction(theta, phi):
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

    def get_angle(self):
        return (self.theta, self.phi)

    def energy(self):
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)
    
    def get_direction(self):
        return (self._dx, self._dy, self._dz)
        
    
    def photoelectric(thetas, phis, energy, x, y, z, t, fano): #photons will have a chance to undergo the photoelectric effect when detected. They create an electron with the same energy as the incident photon (-the binding energy but that is negligible so will be ignored for now)
        cos_theta = np.cos(thetas)
        sin_theta = np.sin(thetas)
        cos_phi = np.cos(phis)
        sin_phi = np.sin(phis)
        P_e = np.sqrt(2*m_e*energy)
        px_e = P_e * cos_theta * cos_phi
        py_e = P_e * sin_theta * cos_phi
        pz_e = P_e * sin_phi
        photoelectron = electron(px_e, py_e, pz_e, thetas, phis, x, y, z, t) #This is assuming that the electron recoils in the same direction as the photon. In reality it will be more complex but I will come to that later
        return photoelectron.get_energy(fano)
        
        
        