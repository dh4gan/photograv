from photograv import make_figure_flight, make_figure_flight_xy_xz, make_video_flight
import vector
import star
import sail
from sail import AU
import matplotlib.pyplot as plt

# Define sail:
nsteps = 9000  # Number of timesteps to compute
timestep = 60 * 5  # 0.1 One timestep every 5 minutes

speed = 1270 # [km/sec], initial speed of spaceship
ship_sail_area = 10  # sail surface im square meters.
ship_mass = .001  # [kg]
ship_charge = 1.0e-6 # Charge in Coulomb

ship_position = vector.Vector3D(3.2*star.R_star_CenA,10.0*AU,1.0) # start position vertical / distance travelled
ship_velocity = vector.Vector3D(0.0,-speed*1000,1.0) # unit conversion; sign: fly downwards

# Create ship object
ship = sail.Sail(ship_mass,ship_sail_area,ship_charge,ship_position,ship_velocity,nsteps)

# Now define star
star_position = vector.Vector3D(0.0,0.0,0.0)
star_velocity = vector.Vector3D(0.0,0.0,0.0)

# Create star object
magmomentvector = vector.Vector3D(0.0,0.0,1.0) # unit vector describing the stellar dipole
#magmomentvector.rotateX(0.6283)
cenA = star.Star(star.M_star_CenA,star.R_star_CenA,star.L_star_CenA,star.sun_Bfield_1AU,star_position,star_velocity, magmom = magmomentvector)
        
afterburner_distance = 10e10  # [R_star]
minimum_distance_from_star = 5 * star.R_star_CenA  # closest distance to star. Only used to check, not to adjust sim!


print '-------------------------------------'
print 'Commencing Single Flight: Parameters'
print '-------------------------------------'
print ship
print cenA

ship.fly(cenA,minimum_distance_from_star,afterburner_distance,timestep,return_mission =False)
ship.write_telemetry_to_file('test.telemetry',cenA)

ship.read_from_telemetry_file('test.telemetry')

#print ship.telemetry['px'][-10:]
#print ship.telemetry['py'][-10:]
#print ship.telemetry['pz'][-10:]

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
#ax1.set_xlim(1.1e6,1.2e6)
ax1.plot(ship.telemetry['time'], ship.telemetry['px'], label='px', color='blue')
ax1.plot(ship.telemetry['time'], ship.telemetry['py'], label='py', color = 'green')
ax1.plot(ship.telemetry['time'], ship.telemetry['pz'], label='pz', color='red', linestyle='dashed')
ax1.legend(loc='upper left')
#ax2 = fig1.add_subplot(312)
#ax2.set_xlim(1.1e6,1.2e6)
#ax2.plot(ship.telemetry['time'], ship.telemetry['sail_x'], label='nx', color='black')
#ax2.plot(ship.telemetry['time'],ship.telemetry['ship_speed'], label='speed', color='black')
#ax2.legend(loc='lower left')
#ax3 = fig1.add_subplot(313)
#ax3.set_xlim(1.1e6,1.2e6)
#ax3.plot(ship.telemetry['time'],ship.telemetry['sail_angle'], label = 'alpha', color='red')
#ax3.legend(loc='upper left')
#ax3.set_xlabel('Time')
plt.show()



# Make figure
fig = plt.gcf()

my_plot = make_figure_flight_xy_xz(
    ship,cenA,
    scale = 20,  # of plot in [stellar radii]
    flight_color='black',  # color of flight trajectorie line
    redness = 0.7,  # 1:red star, 0:yellow; quad limb is darkening hard-coded
    show_burn_circle = False,  # show dashed circle for minimum distance
    star_name = r'$\alpha$ Cen A',
    weight_ratio = 0.1,  # to print 0.1g (we have 10m^2 and 1g)
    circle_spacing_minutes = 120,  # grid of circle marks every N minutes
    annotate_cases = False,  # I, II, III, IV, V
    caption = 'b')  # Figure "suptitle"


plt.show()

my_plot = make_video_flight(
    ship,cenA,
    scale = 20,  # of plot in [stellar radii]
    flight_color='black',  # color of flight trajectorie line
    redness = 0.7,  # 1:red star, 0:yellow; quad limb is darkening hard-coded
    show_burn_circle = False,  # show dashed circle for minimum distance
    star_name = r'$\alpha$ Cen A',
    weight_ratio = 0.1,  # to print 0.1g (we have 10m^2 and 1g)
    circle_spacing_minutes = 120,  # grid of circle marks every N minutes
    annotate_cases = False,  # I, II, III, IV, V
    caption = 'b')  # Figure "suptitle")