from photograv import make_figure_multiple_stars
import vector
import star
import sail
from sail import AU
import matplotlib.pyplot as plt
from numpy import linspace

nstars = 20
stararray = []

rotation_angle = linspace(0.0,1.0, num=nstars)


# Define sail:
nsteps = 10000  # Number of timesteps to compute
timestep = 60 * 5  # 0.1 One timestep every 5 minutes

speed = 1270 # [km/sec], initial speed of spaceship
ship_sail_area = 10  # sail surface im square meters.
ship_mass = .001  # [kg]
ship_charge = 0.0 # Charge in Coulomb

ship_position = vector.Vector3D(2.5*star.R_star_CenA,10.0*AU,0.0) # start position vertical / distance travelled
ship_velocity = vector.Vector3D(0.0,-speed*1000,0.0) # unit conversion; sign: fly downwards

# Create sail
ship = sail.Sail(ship_mass,ship_sail_area,ship_charge,ship_position,ship_velocity,nsteps)

# Now define star
star_position = vector.Vector3D(0.0,0.0,0.0)
star_velocity = vector.Vector3D(0.0,0.0,0.0)

# Create star object

for istar in range(nstars):
    cenA = star.Star(star.M_star_CenA,star.R_star_CenA,star.L_star_CenA,star.sun_Bfield_1AU,star_position,star_velocity)
    cenA.magmom = vector.Vector3D(0.0,0.0,1.0).rotateX(rotation_angle[istar]) # unit vector describing the stellar dipole

    stararray.append(cenA)

afterburner_distance = 10e10  # [R_star]
minimum_distance_from_star = 5 * star.R_star_CenA  # closest distance to star. Only used to check, not to adjust sim!


my_plot = make_figure_multiple_stars(
    ship,stararray,
    minimum_distance_from_star,
    afterburner_distance,
    timestep,
    return_mission=False,
    scale=20,
    caption='multi',
    colorbar=False)



plt.show()
