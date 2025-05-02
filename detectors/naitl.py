import numpy as np
from interfaces import IDetector, IParticle
from particles.photon import photon

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

        
    def detects_batch(self, origins, directions, theta, phi, E):
        xc, yc, zc, xf, yf, zf = self._bounds
        bound_min = np.array([xc, yc, zc])
        bound_max = np.array([xf, yf, zf])

        with np.errstate(divide='ignore', invalid='ignore'):
            t_min = (bound_min - origins) / directions
            t_max = (bound_max - origins) / directions

        t_min = np.where(directions == 0, -np.inf, t_min)
        t_max = np.where(directions == 0, np.inf, t_max)
        tclose = np.max(t_min, axis=1)
        tfar = np.min(t_max, axis=1)
        #Origins is (100000, 3) and I want to pull (100000,x) out of it
        positions = origins + np.dot(tclose, directions)
        positions = positions.reshape(-1,3)
        
        x = np.where(tclose <= tfar, positions[:,0], None)
        y = np.where(tclose <= tfar, positions[:,1], None)
        z = np.where(tclose <= tfar, positions[:,2], None)
        thetas = np.where(tclose <= tfar, theta, None)
        phis = np.where(tclose <= tfar, phi, None)
        batch_detected = np.sum(tclose <= tfar)
        tclose_det = np.where(tclose <= tfar, tclose, None)
        for i in range(100000):
            if tclose[i] <= tfar[i]:
                electron = (photon.photoelectric(thetas[i], phis[i], E, x[i], y[i], z[i], tclose_det[i]))      

        #print(positions)    
        return tclose <= tfar
    
    # def calculate_angles(self):
    #     dO = (self.x, self.y, self.z) #Vector to the technical "origin" of the detector (Closest vertex to the universal origin)
    #     dx = (self.x + self.X, self.y, self.z) 
    #     dy = (self.x, self.y + self.Y, self.z)
    #     dz = (self.x, self.y, self.z + self.Z)
    #     dxy = (self.x + self.X, self.y + self.Y, self.z)
    #     dxz = (self.x + self.X, self.y, self.z + self.Z)
    #     dyz = (self.x, self.y + self.Y, self.z + self.Z)
    #     dxyz = (self.x + self.X, self.y + self.Y, self.z + self.Z)
    #     vectors = np.array([dO, dx, dy, dz, dxy, dxz, dyz, dxyz])
    #     angles = []
    #     #calculate angles to universal origin:
    #     for i in vectors:
    #         angles[i] = (np.arctan(vectors[i][0]/vectors[i][1]), np.arcsin(vectors[i][2])) #array of theta, phi for each vector
    #     #The minimum phi angle will always be the one to dx (or dxy) and the maximum will be the one to dz
    #     #The min theta will be to d0 and the max will be dy
    #     #If minphi<phi<maxphi and mintheta<theta<maxtheta then, for a cuboid, any combination of phi and theta meeting these requirements should lay on the surface of the cuboid
        
        
            
        