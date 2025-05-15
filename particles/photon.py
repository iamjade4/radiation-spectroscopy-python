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
     
    def comptonscatter(theta, phi, energy, x, y, z, t, fano):
        #Ef = Ei/(1 + (Ei/me*c**2)(1-costheta)) where Ef is the energy of the final photon -> therefore energy deposited into the electron is Ei - Ef
        #Theta in the above statement is not the same as the theta taken as an argument. The theta in the calculation is the scattered angle, not the initial angle
        theta_s = np.random.uniform(0, 2*np.pi) #THIS IS NOT THE REAL ANGULAR DISTRIBUTION. THE REAL ONE IS BASED OFF OF THE KLEIN NISHINA FORMULA https://en.wikipedia.org/wiki/Klein%E2%80%93Nishina_formula
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        P_i = energy
        px_i = P_i * cos_theta * cos_phi
        py_i = P_i * sin_theta * cos_phi
        pz_i = P_i * sin_phi #momenta for incident photon
        
        Ef = energy/(1 + (energy/511)*(1-np.cos(theta_s)))
        E_el = energy - Ef
        px_s = Ef * np.cos(theta_s + theta) * cos_phi
        py_s = Ef * np.sin(theta_s + theta) * cos_phi #ignoring the phi component for now, 3D hurts head
        
        #conservation of momentum
        #Pi = Ps + Pe --> Pe = Pi - Ps
        px_e = px_i - px_s
        py_e = py_i - py_s
        P_e = math.sqrt(E_el*2*511)
        pz_e = math.sqrt(abs(P_e**2 - px_e**2 - py_e**2)) #This is negative somehow therefore my particle kinematics is wrong but i need to sleep soon so i will abs this
        theta_e = np.arctan(py_e/px_e)
        phi_e = np.arcsin(pz_e) ##invalid error here apparently
        comptonelectron = electron(px_e, py_e, pz_e, theta_e, phi_e, x, y, z, t)
        return comptonelectron.get_energy(fano)
        
        