import numpy as np
import scipy.constants as const
import random as rd
import sympy as sy
import math

c = const.c
photons = []
detected = 0
#2darray containing the information about all of the photons. Unsure as to whether I should even be using classes or if I should directly append particles to the array

#Energies are in keV, natural units are assumed, distances are in mm for the time being (for the sake of dimensions of detectors)
E = 662 
#Defining the photon class
class photon:
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
    def energy(self):
        E = math.sqrt(self.px**2 + self.py**2 + self.pz**2)
        return E
    #def position(self):
    #    return self.px, self.py, self.pz, self.theta, self.phi, self.x, self.y, self.z, self.t
    def raytrace(self):
        theta=self.theta
        phi=self.phi
        x=self.x
        y=self.y
        z=self.z
        t=self.t
        px=self.px
        py=self.py
        pz=self.pz
        arraytrace = []
        #Normalise Z
        #dt*c*sin(phi) = 1 --> dt = 1/(c*sin(phi))
        dt = 1/(c*np.sin(phi)) #Normalising my vector and then flipping the sign if negative (it only extends into one direction from the origin, not in both directions, so the sign actually matters)
        if pz<0:                #This is a lot quicker, the inefficient part is the Detection method
            dz = -1
            dx = -np.cos(theta)*np.cos(phi)/np.sin(phi)
            dy = -np.sin(theta)*np.cos(phi)/np.sin(phi)
        else:
            dz = 1
            dx = np.cos(theta)*np.cos(phi)/np.sin(phi)
            dy = np.sin(theta)*np.cos(phi)/np.sin(phi)
        # for i in range(1000):
        #     position = [i*dx+x,i*dy+y,i*dz+z]
        #     arraytrace.append(position)
        vector = [dx, dy, dz]
        origin = sy.Point3D(x, y, z)
        path = sy.Ray(origin, direction_ratio=vector)
        return path
        # while abs(x) < 1000 and abs(y) < 1000 and abs(z) < 1000: #creating an arbitrary box of 2000x2000x2000 mm^3 for out of bounds. Will tweak this later, but it just exists now to avoid particle leakage
        #     x += c*dt*np.cos(theta)*np.cos(phi)
        #     y += c*dt*np.sin(theta)*np.cos(phi)
        #     z += c*dt*np.sin(phi)
        #     position = [x, y, z]
        #     arraytrace.append(position)
        #     #Creating an array for each incrimental position along the particle's track. This is very time expensive
        # return arraytrace
        
class NaITl:
    def __init__(self, x, y, z, X, Y, Z):
        #Lowercase are the minumum bound position and uppercase are the dimensions
        #So far this only allows for cuboidal detectors, I will have to redefine it for more complex geometry
        self.x = x
        self.y = y
        self.z = z
        self.X = X
        self.Y = Y
        self.Z = Z
    def position(self):
        #Going to try using planes to define the geometry of the detector
        x = self.x
        y = self.y
        z = self.z
        X = self.X
        Y = self.Y
        Z = self.Z#This only works for a cuboidal structure
        vrt1 = (x, y, z)
        vrt2 = (x + X, y, z)
        vrt3 = (x, y + Y, z)
        vrt4 = (x + X, y + Y, z)
        vrt5 = (x, y, z+Z)
        vrt6 = (x + X, y, z+Z)
        vrt7 = (x, y + Y, z+Z)
        vrt8 = (x + X, y + Y, z+Z)
        print(vrt1)
        plane1 = sy.geometry.polygon.Polygon(vrt1, vrt2, vrt4, vrt3)#POLYGON DOESNT WORK WITH 3D POINTS OKAYYYY BUT PLANES ARE INFINITE
        plane2 = sy.geometry.polygon.Polygon(vrt1, vrt2, vrt6, vrt5)
        plane3 = sy.geometry.polygon.Polygon(vrt1, vrt3, vrt5, vrt7)
        plane4 = sy.geometry.polygon.Polygon(vrt3, vrt4, vrt7, vrt8)
        plane5 = sy.geometry.polygon.Polygon(vrt2, vrt4, vrt6, vrt8)
        plane6 = sy.geometry.polygon.Polygon(vrt5, vrt6, vrt7, vrt8)
        return plane1, plane2, plane3, plane4, plane5, plane6


def Detection(detector, particle):#The "detector" argument takes a list of planes for a detector
    for i in range(len(detector)):
        intr = detector[i].intersection(particle.raytrace())
        print(intr)
    # for i in range(len(particle.raytrace())):
    #     if detector.position()[0] < particle.raytrace()[i][0] < detector.position()[0] + detector.position()[3] and detector.position()[1] < particle.raytrace()[i][1] < detector.position()[1] + detector.position()[4] and detector.position()[2] < particle.raytrace()[i][2] < detector.position()[2] + detector.position()[5]: 
    #          #cross checking each term in the array arraytrace against the detector's coordinates. This is incredibly inefficient
    #          detection = 1
    #          break
    #     else:
    #         detection = 0
    
    detection = 0
    return detection

#Creating the detector once to avoid creating it 1000 times
detector = NaITl(100, 100, 0, 300, 300, 300)
planes = detector.position()

for i in range(1000):
    #Randomly generated angles however, it will be generating a duplicate for 0 & 2pi, I could find a solution for this (go up to 2pi-1x10^-16 for flt16)
    theta = rd.uniform(0, 2*np.pi)
    phi = rd.uniform(-np.pi/2, np.pi/2)
    px = E*np.cos(theta)*np.cos(phi)
    py = E*np.sin(theta)*np.cos(phi)
    pz = E*np.sin(phi)
    #General momentum calculation, energy is based off of Cs 137 gammas. An actual source will not be mono-energetic (other radiation branches)
    Photon = photon(px, py, pz, theta, phi, 0, 0, 0, i/10000)
    #position at 0 0 0 since for now all photons are originating at the origin of the source
    photons.append(Photon.raytrace())
    #Need to figure out sleep command so that it actually generates particles in a set timeframe. Also, the particles shouldn't be equally spaced out time-wise, they're supposed to represent a randomly disintegrating source
    detected += Detection(planes, Photon)

print(detected)