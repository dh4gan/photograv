from photograv import *

# Define spaceship:
number_of_steps = 100000  # 300000
timestep = 60 * 5  # 0.1 One timestep every 5 minutes
ship_sail_area = 10  # sail surface im square meters. Weight is 1g hard-coded.
afterburner_distance = 10e10  # [R_star]
speed = 1270 # [km/sec], initial speed of spaceship
minimum_distance_from_star = 5 * R_star_CenA  # closest distance to star. Only used to check, not to adjust sim!

ship_mass = .001  # [kg]

ship_position = vector.Vector3D(2.8*R_star_CenA,10.0*AU,0.0) # start position vertical / distance travelled
ship_velocity = vector.Vector3D(0.0,-speed*1000,0.0) # unit conversion; sign: fly downwards

ship_charge = 1
magmom_star_CenA = vector.Vector3D(0.0,0.0,0.0).scalarmult(sun_magnetic_moment)
#magmom_star_CenA.rotateX(0.85)
star_position = vector.Vector3D(0.0,0.0,0.0)


data = fly(
    ship_position,
    ship_velocity,
    ship_mass,
    ship_sail_area,
    ship_charge,
    M_star_CenA,
    R_star_CenA,
    L_star_CenA,
    magmom_star_CenA,
    star_position,
    minimum_distance_from_star,
    afterburner_distance,
    timestep,
    number_of_steps,
    return_mission=False)

fig1 = plt.figure()
ax1 = fig1.add_subplot(311)
ax2 = fig1.add_subplot(312)
ax3 = fig1.add_subplot(313)
ax1.plot(data['time'],data['sail_angle'], label='sailangle')
#ax1.plot(data['time'],data['sail_y'],label='saily')
ax2.plot(data['time'],data['F_photon_x'],label='Fx')
ax2.plot(data['time'],data['F_photon_y'],label='Fy')
ax2.plot(data['time'],data['F_mag_x'], label='Fmagx')
ax2.plot(data['time'],data['F_mag_y'], label='Fmagy')
#ax2.plot(data['time'],data['F_mag_z'], label='Fmagz')
ax3.plot(data['time'],data['px'],label='x')
ax3.plot(data['time'],data['py'],label='y')
ax3.plot(data['time'],data['pz'],label='z')
#ax1.set_xlim(1180000,1220000)
#ax2.set_xlim(1180000,1220000)
#ax3.set_xlim(1180000,1220000)
#ax1.plot(data['time'],data['vx'],label='vx')
#ax1.plot(data['time'],data['vy'],label='vy')
ax1.legend(loc='lower left')
ax2.legend(loc='lower left')
ax3.legend(loc='lower left')

plt.show()


# Make figure
fig = plt.gcf()
# my_plot = make_video_flight(
#         data,
#         stellar_radius=R_star_CenA,
#         scale=20,
#         flight_color='black',
#         redness=0.7,
#         show_burn_circle=False,
#         star_name=r'$\alpha$ Cen A',
#         weight_ratio=0.1,
#         circle_spacing_minutes=120,
#         annotate_cases=False,
#         caption='movie')




my_plot = make_figure_flight(
    data = data,
    stellar_radius = R_star_CenA,
    scale = 20,  # of plot in [stellar radii]
    flight_color='black',  # color of flight trajectorie line
    redness = 0.7,  # 1:red star, 0:yellow; quad limb is darkening hard-coded
    show_burn_circle = False,  # show dashed circle for minimum distance
    star_name = r'$\alpha$ Cen A',
    weight_ratio = 0.1,  # to print 0.1g (we have 10m^2 and 1g)
    circle_spacing_minutes = 120,  # grid of circle marks every N minutes
    annotate_cases = False,  # I, II, III, IV, V
    caption = 'b')  # Figure "suptitle"

#fig.savefig("2b.pdf", bbox_inches = 'tight')
plt.show()