import numpy as np
import math
import random
from interfaces import IParticle
from particles.electron import electron
from scipy import stats

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
    def I(theta_s, energy):#For the crosssection from the Klein Nishina formula
        R = energy/m_e
        return -np.cos(theta_s) / R**2 + \
               np.log(1 + R * (1 - np.cos(theta_s))) * (1/R - 2/R**2 - 2/R**3) - \
               1 / (2 * R * (1 + R * (1 - np.cos(theta_s)))**2) + \
               1 / (1 + R * (1 - np.cos(theta_s))) * (-2/R**2 - 1/R**3)
    def E_s(theta_s, energy):
        return(energy/(1 + (energy/511)*(1-np.cos(theta_s))))
    def f(theta_s, energy, E_s):
        return ((1/(photon.I(np.pi, energy)-photon.I(0, energy)))*((E_s/energy)**2)*((energy/E_s)+(E_s/energy)-np.sin(theta_s)**2)*np.sin(theta_s))/31830.6703026283 #probability density function - This also didnt normalise properly? so i am dividing it to normalise it although it should already be normalised
    def get_origin(self):
        return (self.x, self.y, self.z)

    def get_angle(self):
        return (self.theta, self.phi)

    def energy(self):
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)

    def get_direction(self):
        return (self._dx, self._dy, self._dz)
    def gen_angles(energy, batch_size):
        theta_s = np.linspace(0, np.pi, batch_size)
        angle_energies = np.linspace(energy, energy, batch_size)
        dstr = photon.f(theta_s, angle_energies, photon.E_s(theta_s, angle_energies))
        angle_array = random.choices(theta_s, dstr, k=batch_size) #Doing it in bulk here instead of parsing it through to comptonscatter and doing it individually
        return(angle_array)
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

    def comptonscatter(theta, phi, energy, x, y, z, t, fano, angles):
        #Ef = Ei/(1 + (Ei/me*c**2)(1-costheta)) where Ef is the energy of the final photon -> therefore energy deposited into the electron is Ei - Ef
        #Theta in the above statement is not the same as the theta taken as an argument. The theta in the calculation is the scattered angle, not the initial angle
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        cos_phi_s = np.cos(-phi)
        sin_phi_s = np.sin(-phi)
        P_i = energy
        px_i = P_i * cos_theta * cos_phi
        py_i = P_i * sin_theta * cos_phi
        pz_i = P_i * sin_phi #momenta for incident photon
        theta_s = random.choice(angles) #using the defined distribution to choose a weighted angle. This is a bit slow but that's ok
        Ef = energy/(1 + (energy/511)*(1-np.cos(theta_s)))
        Ps = Ef
        pz_s = Ef * sin_phi_s
        px_s = Ef * np.cos(theta_s) * cos_phi_s
        py_s = Ef * np.sin(theta_s) * cos_phi_s
        E_el = energy - Ef
        #P_el = math.sqrt(2*511*E_el) #This is correct i believe
        P_el = math.sqrt(P_i**2-Ps**2)
        #conservation of momentum
        #Pi = Ps + Pe --> Pe = Pi - Ps
        px_el = math.sqrt(abs(px_i**2 - px_s**2))#domain error --> ps > pi ? possibly small rounding error stuff bc everything else seeeems okayy
        py_el = math.sqrt(abs(py_i**2 - py_s**2))
        pz_el = math.sqrt(abs(pz_i**2 - pz_s**2))
        #print(pz_el - pz_i + pz_s)#SHOULD be 0
        #print(P_el**2 - px_el**2 - py_el**2 - pz_el**2) #most of these are smaller than 9, but some are like 500000? not sure what's going on there
        theta_el = np.arctan(py_el/px_el)
        comptonelectron = electron(px_el, py_el, pz_el, theta_el, phi, x, y, z, t)
        return comptonelectron.get_energy(fano)
        
        