import numpy as np
import random as rd
import math

# abstraction layers for particle and detector (i for interface) to follow dependency inversion
class IParticle:
    def get_origin(self):
        pass
    
    def get_direction(self):
        pass

class IDetector:
    def detects(self, particle: IParticle):
        pass

photons = []
detected = 0
#2darray containing the information about all of the photons. Unsure as to whether I should even be using classes or if I should directly append particles to the array

#Energies are in keV, natural units are assumed, distances are in mm for the time being (for the sake of dimensions of detectors)
E = 662 
#Defining the photon class
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


class NaITl(IDetector):
    def __init__(self, x, y, z, X, Y, Z):
        #Lowercase are the minumum bound position and uppercase are the dimensions
        #So far this only allows for cuboidal detectors, I will have to redefine it for more complex geometry
        self.x = x
        self.y = y
        self.z = z
        self.X = X
        self.Y = Y
        self.Z = Z
        self._bounds = self._calculate_bounds()

    def _calculate_bounds(self):
        return (
            self.x, self.y, self.z,
            self.x + self.X,
            self.y + self.Y,
            self.z + self.Z
        )

    def detects(self, particle: IParticle):
        x0, y0, z0 = particle.get_origin()
        dx, dy, dz = particle.get_direction()
        xc, yc, zc, xf, yf, zf = self._bounds

        distc = []
        distf = []
        for origin, bound, direction in zip((x0, y0, z0), (xc, yc, zc), (dx, dy, dz)):
            distc.append((bound - origin) / direction) if direction != 0 else None
        for origin, bound, direction in zip((x0, y0, z0), (xf, yf, zf), (dx, dy, dz)):
            distf.append((bound - origin) / direction) if direction != 0 else None

        if not distc or not distf:
            return False

        tclose = max(distc)
        tfar = min(distf)

        return tclose <= tfar


#Creating the detector once to avoid creating it 1000 times
detector = NaITl(100, 100, 0, 300, 300, 300)

for i in range(100000):
    #Randomly generated angles however, it will be generating a duplicate for 0 & 2pi, I could find a solution for this (go up to 2pi-1x10^-16 for flt16)
    theta = rd.uniform(0, 2*np.pi)
    phi = rd.uniform(-np.pi/2, np.pi/2)
    px = E*np.cos(theta)*np.cos(phi)
    py = E*np.sin(theta)*np.cos(phi)
    pz = E*np.sin(phi)
    #General momentum calculation, energy is based off of Cs 137 gammas. An actual source will not be mono-energetic (other radiation branches)
    Photon = photon(px, py, pz, theta, phi, 0, 0, 0, 0) #replaced i/1000
    #position at 0 0 0 since for now all photons are originating at the origin of the source
    photons.append(Photon)
    #Need to figure out sleep command so that it actually generates particles in a set timeframe. Also, the particles shouldn't be equally spaced out time-wise, they're supposed to represent a randomly disintegrating source
    if detector.detects(Photon):
        detected += 1

print(detected)
