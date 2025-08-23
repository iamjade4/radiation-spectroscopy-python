#This file is home to all of the different types of scintillators
import numpy as np
import pandas as pd
from interfaces import IDetector, IParticle
from particles.photon import photon
from detectors.scintillators.scintillator import Scintillator
import math

class NaITl(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n = 53*0.846627 + 11*0.153373 #fractional atomic numbers
        self.n_m = 1.474*10**22 #molecular number density
        self.n = 2*self.n_m #atomic number density
        super().__init__()

class CsITl(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n = 53*0.488451 + 55*0.511549 #fractional atomic numbers
        self.n_m = 1.024*10**22 #molecular number density
        self.n = 2*self.n_m #atomic number density
        super().__init__()

class LiIEu(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n = 3*0.051858 + 53*0.948142 #fractional atomic numbers
        self.n_m = 1.8066*10**22 
        self.n = 2*self.n_m
        super().__init__()

class BGO(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n =  8*0.154119 + 32*0.174859 + 83*0.671022#fractional atomic numbers
        self.n_m = 0.3447*10**22
        self.n = 19*self.n_m
        super().__init__()

class CdWO(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n =  8*0.177644 + 48*0.312029 + 74*0.510328#fractional atomic numbers
        self.n_m = 0.1321*10**22
        self.n = 6*self.n_m
        super().__init__()

class CaWO(Scintillator):
    def __init__(self, base, axis, radius, height):
        #this is a cylinder
        self.b = base
        self.a = axis
        self.r = radius
        self.h = height
        self.t = self.b + self.a*self.h #top
        self.Z_n =  8*0.222271 + 20*0.139196 + 74*0.638533#fractional atomic numbers
        self.n_m = 0.1276*10**22
        self.n = 6*self.n_m
        super().__init__()