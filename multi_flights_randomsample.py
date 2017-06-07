# Written 6/6/17 by dh4gan
# Carries out a large number of flights,
# Varying sail and stellar parameters slightly


import vector
import star
import sail
from sail import AU
import matplotlib.pyplot as plt
from numpy import pi
from numpy.random import rand

nships = 5
shiparray = []

# Define distributions for each parameter

Bmin = 0.5*star.sun_Bfield_1AU
Bmax = 1.5*star.sun_Bfield_1AU

B_thetamin = -5*pi/180.0
B_thetamax = 5*pi/180.0

B_phimin = -5*pi/180.0
B_phimax = 5*pi/180.0

chargemin = -1.0e-3
chargemax = 0.0
#chargemin = 0.0
#chargemax = 1.0e-3

Lmin = 0.999*star.L_star_CenA
Lmax = 1.001*star.L_star_CenA

xmin = 2.5*star.R_star_CenA
xmax = 3.5*star.R_star_CenA

# Define fixed sail parameters :
nsteps = 5000  # Number of timesteps to compute
timestep = 60 * 10  # 0.1 One timestep every 5 minutes

speed = 1270 # [km/sec], initial speed of spaceship
ship_sail_area = 10  # sail surface im square meters.
ship_mass = .001  # [kg]

afterburner_distance = 10e10  # [R_star]
minimum_distance_from_star = 5 * star.R_star_CenA  # closest distance to star. Only used to check, not to adjust sim!

# Setup output files

outputfile = 'multiflights.output'
f_obj = open(outputfile, 'w')

line = str(nships) + ' \n'
f_obj.write(line)


closest_approach = []
maxforcemag = []

x_final = []
y_final = []
z_final = []

vx_final = []
vy_final = []
vz_final = []
vmag = []

# Loop over total number of flights

for iship in range(nships):

    # Define sail

    #ship_charge = chargemin + rand()*(chargemax-chargemin)
    ship_charge = 0.0
    x = xmin + rand()*(xmax-xmin)

    ship_position = vector.Vector3D(x,10.0*AU,0.0) # start position vertical / distance travelled
    ship_velocity = vector.Vector3D(0.0,-speed*1000,0.0) # unit conversion; sign: fly downwards

    print "sail parameters"
    print "x: ",str(x/star.R_star_CenA)
    print "q: ",str(ship_charge)
    
    ship = sail.Sail(ship_mass,ship_sail_area,ship_charge,ship_position,ship_velocity,nsteps)
    shiparray.append(ship)

    # Now define star
    star_position = vector.Vector3D(0.0,0.0,0.0)
    star_velocity = vector.Vector3D(0.0,0.0,0.0)

    Lstar = Lmin + rand()*(Lmax-Lmin)
    B1AU = Bmin + rand()*(Bmax-Bmin)
    
    B_theta = B_thetamin + rand()*(B_thetamax-B_thetamin)
    B_phi = B_phimin + rand()*(B_phimax-B_phimin)

    print "Star parameters"
    print "L: ", Lstar/star.L_star_CenA
    print "B: ", B1AU/star.sun_Bfield_1AU
    print "Angles: ", B_theta*180.0/pi, B_phi*180.0/pi


    # Create star object
    magmomvector = vector.Vector3D(0.0,0.0,1.0)
    magmomvector.rotateX(B_theta)
    magmomvector.rotateZ(B_phi)
    
    cenA = star.Star(star.M_star_CenA,star.R_star_CenA,Lstar,B1AU,star_position,star_velocity,magmom=magmomvector)
        
    print ship
    print cenA

    ship.fly(cenA,minimum_distance_from_star,afterburner_distance,timestep,return_mission =False)
    filename = 'flight_'+str(iship+1)+'.telemetry'
    ship.write_telemetry_to_file(filename, cenA)
    # Extract data of interest from telemetry

    x_init = ship.telemetry['px'][0]
    y_init = ship.telemetry['py'][0]
    z_init = ship.telemetry['pz'][0]
    
    x_final.append(ship.telemetry['px'][-1])
    y_final.append(ship.telemetry['py'][-1])
    z_final.append(ship.telemetry['pz'][-1])
    
    vmag.append(ship.telemetry['ship_speed'][-1])
    
    vx_final.append(ship.telemetry['vx'][-1]/vmag[-1])
    vy_final.append(ship.telemetry['vy'][-1]/vmag[-1])
    vz_final.append(ship.telemetry['vz'][-1]/vmag[-1])
    

    
    encounter_time, close_step = ship.get_closest_encounter()
    closest_approach.append(ship.telemetry['stellar_distance'][close_step])
    
    maxforce_time, maxforce_step = ship.get_maximum_total_force()
    maxforcemag.append(ship.telemetry['F_total'][maxforce_step])
    maxforce_x = ship.telemetry['F_tot_x'][maxforce_step]
    maxforce_y = ship.telemetry['F_tot_y'][maxforce_step]
    maxforce_z = ship.telemetry['F_tot_z'][maxforce_step]

#     my_plot = make_figure_flight(
#     ship,cenA,
#     scale = 20,  # of plot in [stellar radii]
#     flight_color='black',  # color of flight trajectorie line
#     redness = 0.7,  # 1:red star, 0:yellow; quad limb is darkening hard-coded
#     show_burn_circle = False,  # show dashed circle for minimum distance
#     star_name = r'$\alpha$ Cen A',
#     weight_ratio = 0.1,  # to print 0.1g (we have 10m^2 and 1g)
#     circle_spacing_minutes = 120,  # grid of circle marks every N minutes
#     annotate_cases = False,  # I, II, III, IV, V
#     caption = 'b')  # Figure "suptitle"
#     
# 
#     plt.show()
#     
    # Record basic output data
    line = str(ship_charge)+ ' '+str(x_init)+ ' '+str(Lstar/star.L_star_CenA) + ' '+ str(B1AU)+ ' '+str(B_theta)+' '+str(B_phi)
    line = line + ' ' +str(closest_approach)+ ' '+str(x_final) + ' '+str(y_final)+' '+str(z_final)+ ' ' +str(vx_final) + ' '+str(vy_final)+' '+str(vz_final)
    line = line + ' '+str(maxforcemag) + ' '+str(maxforce_x)+ ' '+str(maxforce_y)+' '+str(maxforce_z) +' \n'
    
    
    f_obj.write(line)
    
f_obj.close()
    

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for iship in range(nships):
    ax1.plot(shiparray[iship].telemetry['px'], shiparray[iship].telemetry['py'])

ax1.quiver(x_final,y_final,vx_final,vy_final)
plt.show()
    

