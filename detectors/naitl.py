import numpy as np
from interfaces import IDetector, IParticle

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