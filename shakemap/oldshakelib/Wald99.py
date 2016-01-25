#!/usr/bin/env python

#third party imports
import numpy as np

class Wald99(GMICE):
    def __init__(self,magnitude):
        self.magnitude = magnitude

    def convertToIntensity(self,pga,pgv):
        Imm = 2.20 * np.log(pga) + 1.00
        Imm[Imm >= 7] = 2.10 * np.log(pgv) + 3.40
        Imm[Imm < 2.1987] = 3.66 * np.log(pga) - 1.66
        Imm[(Imm >= 2.1987) && (Imm < 7)] = 3.47 * np.log(pgv) + 2.35
        return Imm

    def convertToPGM(self,Imm):
        pga[(Imm >= 5) && (Imm <= 8)] = np.exp((Imm + 1.66)/3.66)
        return (pga,pgv)
                
                

    
