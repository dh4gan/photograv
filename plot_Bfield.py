# Written 17/5/17 by dh4gan
# Script tests the magnetic field implementation
# Computes the field on a cartesian grid and plots the results

import vector
import numpy as np
import matplotlib.pyplot as plt
import star

starposition = vector.Vector3D(0.0,0.0,0.0)
starvelocity = vector.Vector3D(0.0,0.0,0.0)
magmoment = vector.Vector3D(0.0,0.0,1.0) # unit vector describing the stellar dipole

magmoment.rotateX(0.5*np.pi)

# Create star (Alpha Cen A Parameters)
s = star.Star(star.M_star_CenA,star.R_star_CenA,star.L_star_CenA,star.sun_Bfield_1AU,starposition,starvelocity,magmom=magmoment)

npoints = 100

AU = star.AU

# Define grid for plotting field

px = np.linspace(-2.0*AU, 2.0*AU, num=npoints)
py = np.linspace(-2.0*AU, 2.0*AU, num=npoints)
pz = np.linspace(-2.0*AU, 2.0*AU, num=npoints)

Bx = np.zeros((npoints,npoints))
By = np.zeros((npoints,npoints))
Bz = np.zeros((npoints,npoints))
Bmag = np.zeros((npoints,npoints))

print 'Computing field in z=0 plane'

for i in range(npoints):
    for j in range(npoints):        
            
            position = vector.Vector3D(px[i],py[j],0.0)

            Bfield = s.get_magnetic_field_dipole(position)
            
            Bx[j,i] = Bfield.x
            By[j,i] = Bfield.y            
            Bmag[j,i] = Bfield.mag()
       
       

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.streamplot(px/AU,py/AU, Bx,By)
ax1.scatter(s.position.x/AU, s.position.y/AU, color='red')
ax1.set_xlabel('x (AU)',fontsize=18)
ax1.set_ylabel('y (AU)',fontsize=18)

print 'Computing field in y=0 plane'

for i in range(npoints):
    for j in range(npoints):        
            
            position = vector.Vector3D(px[i],0.0,pz[j])

            Bfield = s.get_magnetic_field_dipole(position)
            
            Bx[j,i] = Bfield.x
            Bz[j,i] = Bfield.z            
            Bmag[j,i] = Bfield.mag()

ax2 = fig1.add_subplot(212)
ax2.streamplot(px/AU,pz/AU, Bx,Bz)
ax2.scatter(s.position.x/AU, s.position.z/AU, color='red')
ax2.set_xlabel('x (AU)', fontsize=18)
ax2.set_ylabel('z (AU)',fontsize=18)
plt.show()

Bmag = Bmag/s.B1AU

print np.amin(Bmag), np.amax(Bmag)

# contourlevels = np.logspace(-3,2, num=30)
# #contourlevels = contourlevels*Bfield_1AU
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# CS = ax2.contour(px/AU,py/AU,Bmag, levels=contourlevels, vmin = 1.0e-4, vmax = 20.0)
# plt.clabel(CS, inline=1, fontsize=10)
# 
# plt.show()