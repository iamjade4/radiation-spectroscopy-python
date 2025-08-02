import numpy as np
import math
import random
from interfaces import IParticle
from particles.electron import electron
from scipy import stats

m_e = 511 #restmass energy

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
    def f(theta_s, energy, E_s):
        energy = energy*1000
        E_s = E_s*1000 #putting them in eV as I believe that's correct
        return ((1/(photon.I(np.pi, energy)-photon.I(0, energy)))*((E_s/energy)**2)*((energy/E_s)+(E_s/energy)-np.sin(theta_s)**2)*np.sin(theta_s)) #probability density function - This also didnt normalise properly? so i am dividing it to normalise it although it should already be normalised
    def I(theta_s, energy):#For the crosssection from the Klein Nishina formula
        R = energy/m_e
        return -np.cos(theta_s) / R**2 + \
               np.log(1 + R * (1 - np.cos(theta_s))) * (1/R - 2/R**2 - 2/R**3) - \
               1 / (2 * R * (1 + R * (1 - np.cos(theta_s)))**2) + \
               1 / (1 + R * (1 - np.cos(theta_s))) * (-2/R**2 - 1/R**3)
    def E_s(theta_s, energy):
        return(energy/(1 + (energy/511)*(1-np.cos(theta_s))))
    def get_origin(self):
        return (self.x, self.y, self.z)

    def get_angle(self):
        return (self.theta, self.phi)

    def energy(self):
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)

    def get_direction(self):
        return (self._dx, self._dy, self._dz)
    def comptonscatter_csc(Ei):
        epsilon = Ei/m_e
        r_es = 0.074 #B
        #FIXED THIS
        return 2*np.pi*(r_es)*((1+epsilon)/epsilon**3 * ((2*epsilon*(1 + epsilon))/(1+2*epsilon)-math.log(1+2*epsilon))+math.log(1+2*epsilon)/2*epsilon -(1+3*epsilon)/((1+2*epsilon)**2))
    def gen_angles(energy, batch_size):
        theta_s = np.linspace(0, np.pi, batch_size)
        angle_energies = np.linspace(energy, energy, num = batch_size)
        dstr = photon.f(theta_s, angle_energies, photon.E_s(theta_s, angle_energies))
        dstr = dstr/sum(dstr) #normalise
        angle_array = np.random.choice(theta_s, size=batch_size, p=dstr) #changed to np.random.choice from random.choices due to an internal value error in random
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
    def photoelectric_csc_mid(Z, Ei):#crosssection for photoelectric effect
        epsilon = Ei/m_e    
        gamma = (Ei + m_e)/m_e #ignoring E_b
        Ei = Ei*1000 #keV to eV
        return 3/2 * 0.66526 * (1/137**4) * ((Z/epsilon)**5)*((gamma**2 - 1)**3/2)*(4/3 + (gamma*(gamma-2))/(gamma+1) * (1 - 1/(2*gamma*math.sqrt(gamma**2 - 1)) *math.log((gamma + math.sqrt(gamma**2-1))/(gamma - math.sqrt(gamma**2-1)))))
    def photoelectric_csc_high(Z, Ei):
        epsilon = Ei/m_e 
        a = 1.6268*10**-9, 1.5274*10**-9, 1.1330*10**-9, -9.12*10**-11
        b = -2.683*10**-12, -5.110*10**-13, -2.177*10**-12, 0
        c = 4.173*10**-2, 1.027*10**-2, 2.013*10**-2, 0
        p = 1, 2, 3.5, 4
        csc = 0
        for i in range(4):
            csc += (Z**5)*((a[i]+b[i]*Z)/(1+c[i]*Z)*epsilon**-p[i])
        return csc
    def comptonscatter(theta, phi, energy, x, y, z, t, fano, angles):
        #Ef = Ei/(1 + (Ei/me*c**2)(1-costheta)) where Ef is the energy of the final photon -> therefore energy deposited into the electron is Ei - Ef
        #Theta in the above statement is not the same as the theta taken as an argument. The theta in the calculation is the scattered angle, not the initial angle
        if theta < 0:
            theta+=2*np.pi
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        P_i = energy
        px_i = P_i * cos_theta * cos_phi
        py_i = P_i * sin_theta * cos_phi
        pz_i = P_i * sin_phi #momenta for incident photon
        theta_s = random.choice(angles) #using the defined distribution to choose a weighted angle. This is a bit slow but that's ok
        Ef = energy/(1 + (energy/511)*(1-np.cos(theta_s)))
        Ps = Ef

        pX_s = Ef * np.cos(theta_s)#These are not on the same axis as px_i and py_i. These are on a plane that follows P_i. P_i is parallel to pX_s
        pY_s = Ef * np.sin(theta_s)

        E_el = (energy - Ef) + m_e #E**2 = P**2 + m**2 where P here is energy and Ef
        KE_el = energy-Ef
        P_el = math.sqrt((KE_el)**2+2*m_e*KE_el) #This is correct i believe

        px_s = pX_s*cos_theta*cos_phi - pY_s*sin_theta#converting from the cme frame to the origin frame
        py_s = pX_s*sin_theta*cos_phi + pY_s*cos_theta
        pz_s = pX_s*sin_phi
        mag = np.sqrt(px_s**2 + py_s**2 + pz_s**2)
        
        px_s = px_s/mag #normalising here
        py_s = py_s/mag
        pz_s = pz_s/mag
        
        px_el = px_i - px_s
        py_el = py_i - py_s
        pz_el = pz_i - pz_s #linear conservation of momentum
        theta_el = np.arctan(py_el/px_el)
        comptonelectron = electron(px_el, py_el, pz_el, theta_el, phi, x, y, z, t)
        theta_s_ = np.arctan(py_s/px_s)
        phi_s = np.arcsin(pz_s/pX_s)
        #print(py_s, px_s, theta_s, phi_s)
        #if abs(P_el**2 - px_el**2 - py_el**2 - pz_el**2) > 1: 
        #    return KE_el, px_s, py_s, pz_s, Ef, theta_s_, phi_s#, total, bad
        return KE_el, px_s, py_s, pz_s, Ef, theta_s_, phi_s#, total, bad #returning KE_el since that's the energy that gets deposited by the electron

        
        