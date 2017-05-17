# Written 17/5/17 by dh4gan
# Script tests the magnetic field implementation
# Computes the field on a cartesian grid and plots the results

import vector
import numpy as np
import matplotlib.pyplot as plt
from photograv import get_magnetic_field_dipole, sun_Bfield_1AU, AU

starposition = vector.Vector3D(0.0,0.0,0.0)
magmoment = vector.Vector3D(0.0,1.0,0.0) # unit vector describing the stellar dipole
Bfield_1AU = sun_Bfield_1AU

npoints = 300

px = np.linspace(-2.0*AU, 2.0*AU, num=npoints)
py = np.linspace(-2.0*AU, 2.0*AU, num=npoints)
#pz = np.linspace(-10.0*AU, 10.0*AU, num=npoints)

Bx = np.zeros((npoints,npoints))
By = np.zeros((npoints,npoints))
Bmag = np.zeros((npoints,npoints))

for i in range(npoints):
    for j in range(npoints):
            
        position = vector.Vector3D(px[i],py[j],0.0)

        Bfield = get_magnetic_field_dipole(position, Bfield_1AU, magmoment,starposition)
            
        Bx[j,i] = Bfield.x
        By[j,i] = Bfield.y
        Bmag[j,i] = Bfield.mag()
       
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.streamplot(px/AU,py/AU, Bx,By)
ax1.scatter(starposition.x/AU, starposition.y/AU, color='red')

Bmag = Bmag/Bfield_1AU

print np.amin(Bmag), np.amax(Bmag)

contourlevels = np.logspace(-3,2, num=30)
#contourlevels = contourlevels*Bfield_1AU
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
CS = ax2.contour(px/AU,py/AU,Bmag, levels=contourlevels, vmin = 1.0e-4, vmax = 20.0)
plt.clabel(CS, inline=1, fontsize=10)

plt.show()