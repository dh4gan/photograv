# Written 23/05/17 by dh4gan
# Class for the Sail Object

import vector
from numpy import cos

# Some useful physical constants

AU = AU = (149.6e6 * 1000)

# Sun
sun_radius = 695700000  # [m]
sun_mass = 1.989 * 10**30  # [kg]
sun_luminosity = 3.86 * 10**26  # [Watt] stellar luminosity
sun_Bfield_1AU = 5.0e-9 # Solar magnetic field strength at 1 AU (Tesla)

# CenA:
L_star_CenA = sun_luminosity * 1.522
R_star_CenA = sun_radius * 1.224
M_star_CenA = sun_mass * 1.105

# CenB:
L_star_CenB = sun_luminosity * 0.503
R_star_CenB = sun_radius * 0.863
M_star_CenB = sun_mass * 0.934

# CenC:
L_star_CenC = sun_luminosity * 138 * 10e-6
R_star_CenC = sun_radius * 0.145
M_star_CenC = sun_mass * 0.123



class Star(object):
    
    
    def __init__(self,m,R,L,B,pos,vel,magmom=vector.Vector3D(0.0,0.0,1.0)):
        '''Initialises star with mass, radius, B-field, position, velocity'''
        self.M = m
        self.R = R
        self.L = L
        self.B1AU = B
        self.position = pos
        self.velocity = vel
        self.magmoment = magmom
        
    def __str__(self):
        s= 'Star: mass %e radius %e luminosity %e\n' % (self.M/sun_mass, self.R/sun_radius, self.L/sun_luminosity)
        s = s+"Position: "+str(self.position)+"\n"
        s = s+"Velocity: "+str(self.velocity)+"\n"
        s = s+"Mag Moment: "+str(self.magmoment)+"\n"
        return s
    
    def get_magnetic_field_dipole(self,position):
        '''Returns a spherically symmetric dipole magnetic field
        NB: Calculated in 3D'''

        sepvector = position.subtract(self.position).scalarmult(1.0/AU)
    
        sep = sepvector.mag()
        sep2 = sep*sep
        sep3 = sep2*sep
    
        sepvector = sepvector.unitVector()
    
        mdotr = self.magmoment.dot(sepvector)
    
        prefac = self.B1AU/sep3
    
        Bfield = sepvector.scalarmult(3.0*prefac*mdotr)
        Bfield = Bfield.subtract(self.magmoment.scalarmult(prefac))

        return Bfield
    
    
    