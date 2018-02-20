# Written 23/05/17 by dh4gan
# Class for the Sail Object

import vector
import numpy as np
from scipy import optimize
import star
from star import sun_radius
import linecache

G = 6.67428e-11  # the gravitational constant G
c = 299792458  # [m/sec] speed of light
AU = (149.6e6 * 1000)  # [m] 149.6 million km
mu0 = 1e-7*4.0*np.pi # permeability of free space (SI units)
epsilon0 = 8.85418782e-12 # permittivity of free space (SI units)

class Sail(object):
    
    
    def __init__(self,m,A,c,pos,vel,steps):
        '''Initialises sail with mass,area, charge, position, velocity'''
        self.mass = m
        self.area = A
        self.charge = c
        self.position = pos
        self.velocity = vel
        self.sail_normal = vector.Vector3D(0.0,0.0,0.0)
        self.F_grav = vector.Vector3D(0.0,0.0,0.0)
        self.F_photon = vector.Vector3D(0.0,0.0,0.0)
        self.F_mag = vector.Vector3D(0.0,0.0,0.0)
        self.F_tot = vector.Vector3D(0.0,0.0,0.0)
        self.nsteps = steps
        
        
        # Initialise array to store self.telemetry for output
        
        self.telemetry = np.zeros((self.nsteps,), dtype=[
                               ('step', 'int32'),
                               ('time', 'f8'),
                               ('px', 'f8'),
                               ('py', 'f8'),
                               ('pz', 'f8'),
                               ('vx', 'f8'),
                               ('vy', 'f8'),
                               ('vz', 'f8'),
                               ('sail_x', 'f8'),
                               ('sail_y', 'f8'),
                               ('sail_z', 'f8'),
                               ('sail_angle', 'f8'),
                               ('F_gravity', 'f8'),
                               ('F_photon', 'f8'),
                               ('F_magnetic', 'f8'),
                               ('F_total', 'f8'),
                               ('F_photon_x','f8'),
                               ('F_photon_y','f8'),
                               ('F_photon_z','f8'),
                               ('F_mag_x','f8'),
                               ('F_mag_y','f8'),
                               ('F_mag_z','f8'),
                               ('F_grav_x','f8'),
                               ('F_grav_y','f8'),
                               ('F_grav_z','f8'),
                               ('F_tot_x','f8'),
                               ('F_tot_y','f8'),
                               ('F_tot_z','f8'),
                               ('photon_acceleration', 'f8'),
                               ('ship_speed', 'f8'),
                               ('stellar_distance', 'f8'),
                               ('sail_not_parallel', 'bool')])
        self.telemetry[:] = np.NAN
        
    def __str__(self):
        s= 'Sail: mass %e area %e charge %e\n' % (self.mass, self.area, self.charge)
        s = s+"Position: "+str(self.position)+"\n"
        s = s+"Velocity: "+str(self.velocity)+"\n"
        s = s+"Sail Normal: "+str(self.sail_normal)+"\n"
        return s
    
    def clone(self):
        
        return Sail(self.mass,self.area,self.charge,self.position,self.velocity,self.nsteps)

    def get_gravity_force(self,star):
        """Return the gravity force due to a star"""

        relpos = self.position.subtract(star.position)
        distance = relpos.mag()
        force = -G * self.mass * star.M / (distance**3)
        
        self.F_grav = relpos.scalarmult(force)
         

    
    def optimise_sail(self,star):
        """Finds the sail normal (n) that optimises the craft's deceleration
        Optimised when (n.r)(n.v) is minimised"""
    
        def sailfunction(nvalues, position,velocity):
        
            r = position.unitVector()
            v = velocity.unitVector()
            n = vector.Vector3D(nvalues[0],nvalues[1],nvalues[2]).unitVector()
        
            return n.dot(v)*n.dot(r)
    
        firstguess = [1.0,0.0,0.0]

        n_optimal = optimize.minimize(sailfunction, firstguess,args=(self.position,self.velocity), method='SLSQP',bounds=((-1,1),(-1,1),(-1,1)), tol=1.0e-10)
    
        sail_normal = vector.Vector3D(n_optimal.x[0],n_optimal.x[1],n_optimal.x[2])
        sail_normal = sail_normal.unitVector()

        return sail_normal

    def get_photon_force(self,star):
        """Returns the photon force in 3D"""
   
        r = self.position.subtract(star.position).mag() # distance between sail and star in m

        # Radial photon pressure force
        F_photon_r = star.L * self.area / (3 * np.pi * c * star.R**2) * \
            (1 - (1 - (star.R / r)**2)**(1.5))

        # Now need to adjust sail alignment
        # Maximum deceleration when sail normal minimises (n.r)(n.v)
    
        self.sail_normal = self.optimise_sail(star)
    
        rdotn = self.sail_normal.dot(self.position.unitVector())    
        self.F_photon = self.sail_normal.scalarmult(rdotn*F_photon_r)


    def get_magnetic_force(self,star):
        '''Computes the Lorentz force from a stellar magnetic field on the charged sail'''
    
        Bfield = star.get_magnetic_field_dipole(self.position)
    
        self.F_mag = self.velocity.cross(Bfield).scalarmult(self.charge)
        
        
    def integrate(self, timestep):
        '''Integrates the system given a total force F'''
    
        newvelocity = self.velocity.add(self.F_tot.scalarmult(timestep/self.mass))
        newposition = self.position.add(newvelocity.scalarmult(timestep))
    
        self.velocity = newvelocity
        self.position = newposition
        
        
    def fly(self,star,minimum_distance_from_star,afterburner_distance,timestep,return_mission =False):
        """Loops through the simulation, returns result array"""

    
        self.telemetry[:] = np.NAN
        deceleration_phase = True

        print 'Beginning flight'
        # Main loop
        for step in range(self.nsteps):
        
            # Check distance from star
            star_distance = self.position.subtract(star.position).mag()
        
            # Inferno time: If we are inside the star, the simulation ends
            if star_distance < star.R:
                print('Exit due to ship being inside the star.')
                break

            self.F_tot = vector.Vector3D(0.0,0.0,0.0)
    
            #################
            # Compute forces
            #################
        
            # Gravity force
            self.get_gravity_force(star)

            self.F_tot = self.F_tot.add(self.F_grav)
        
            # Magnetic Force is velocity dependent!
            self.get_magnetic_force(star)
            
            self.F_tot = self.F_tot.add(self.F_mag)

            # Now photon pressure force
        
            # Check if we are past closest encounter. If yes, switch sail off
            previous_distance = self.telemetry['stellar_distance'][step-1]*star.R
            if step > 2 and star_distance > previous_distance:
                #if deceleration_phase:print "Past closest approach at time ",str(step),str(self.telemetry['time'][step]),": disengaging sail"
                deceleration_phase = False

            # Check if we are past the star and at afterburner distance
            # If yes, switch sail on again
            if not return_mission and self.position.y < 0 and star_distance > afterburner_distance:
                #print "At Afterburner distance"
                deceleration_phase = True  # Actually acceleration now!

                # In case we are inside the minimum distance, the simulation ends
            if star_distance < minimum_distance_from_star / star.R:
                print('Exit due to ship being inside the minimum distance')
                break

            # Special case return mission
            # This is an ugly hardcoded special case. To be optimized in the future
            if return_mission and self.position.x < 0 and self.position.y / star.R > 4.75:
                deceleration_phase = True  # Actually acceleration now!

            # Photon pressure force
            self.get_photon_force(star)
    
            if deceleration_phase:
                self.F_tot = self.F_tot.add(self.F_photon)
                
            # If we do not decelerate: sail shall be parallel with zero photon force
            if not deceleration_phase:
                self.sail_normal = self.position.cross(self.velocity).unitVector() # sail normal to both position and velocity (zero force)
                self.F_photon = vector.Vector3D(0.0,0.0,0.0)

            # Update positions
            self.integrate(timestep)

            """Calculate interesting values"""
        
            # Write interesting values into return array
            self.telemetry['step'][step] = step
            self.telemetry['time'][step] = step * timestep
            self.telemetry['px'][step] = self.position.x / sun_radius
            self.telemetry['py'][step] = self.position.y / sun_radius
            self.telemetry['pz'][step] = self.position.z / sun_radius
            self.telemetry['vx'][step] = self.velocity.x / 1000
            self.telemetry['vy'][step] = self.velocity.y / 1000
            self.telemetry['vz'][step] = self.velocity.z / 1000
            self.telemetry['F_gravity'][step] = self.F_grav.mag()
            self.telemetry['F_photon'][step] = self.F_photon.mag()
            self.telemetry['F_magnetic'][step] = self.F_mag.mag()
            self.telemetry['F_total'][step] = self.F_tot.mag()
            self.telemetry['F_grav_x'][step] = self.F_grav.x
            self.telemetry['F_grav_y'][step] = self.F_grav.y
            self.telemetry['F_grav_z'][step] = self.F_grav.z
            self.telemetry['F_photon_x'][step] = self.F_photon.x
            self.telemetry['F_photon_y'][step] = self.F_photon.y
            self.telemetry['F_photon_z'][step] = self.F_photon.z
            self.telemetry['F_mag_x'][step] = self.F_mag.x
            self.telemetry['F_mag_y'][step] = self.F_mag.y
            self.telemetry['F_mag_z'][step] = self.F_mag.z
            self.telemetry['F_tot_x'][step] = self.F_tot.x
            self.telemetry['F_tot_y'][step] = self.F_tot.y
            self.telemetry['F_tot_z'][step] = self.F_tot.z
            self.telemetry['sail_x'][step] = self.sail_normal.x
            self.telemetry['sail_y'][step] = self.sail_normal.y
            self.telemetry['sail_z'][step] = self.sail_normal.z
            self.telemetry['sail_angle'][step] = np.arccos(self.sail_normal.dot(self.position.unitVector()))    
            self.telemetry['photon_acceleration'][step] = self.F_photon.mag()/self.mass
            self.telemetry['ship_speed'][step] = self.velocity.mag()/1000.
            self.telemetry['stellar_distance'][step] = star_distance/star.R
            self.telemetry['sail_not_parallel'][step] = deceleration_phase

        print 'Flight complete'

    
    def print_flight_report(self):
        '''Prints a brief report of the flight telemetry'''
        
        print('speed [km/sec] start', '{:10.1f}'.format(self.telemetry['ship_speed'][1]),
        ' - speed end', '{:10.1f}'.format(self.telemetry['ship_speed'][-1]))
        print('Overall closest encounter to star [stellar radii]',
        '{:1.3f}'.format(np.amin(self.telemetry['stellar_distance'])))

        index_min = np.argmin(self.telemetry['stellar_distance'])
        print('Closest encounter at time [min]',
              '{:1.3f}'.format(np.amin(self.telemetry['time'][index_min] / 60)))
        print('Closest encounter at speed [km/s]',
              '{:1.3f}'.format(np.amin(self.telemetry['ship_speed'][index_min])))

        angle = abs(np.arctan2(self.telemetry['py'][-1], self.telemetry['px'][-1]) * (360 / (2 * np.pi))) - 90
        print('Deflection angle [degree]',
              '{:1.1f}'.format(angle), '{:1.1f}'.format(180 - angle))
        print('Total time range is [seconds]',
              self.telemetry['time'][0], 'to', self.telemetry['time'][-1])
    
        has_sail_switched_off = False
        print(len(self.telemetry['time']))
        
        for step in range(len(self.telemetry['time'])):
            speed_change = self.telemetry['ship_speed'][step] - self.telemetry['ship_speed'][step-1]  # km/sec
            time_change = self.telemetry['time'][step] - self.telemetry['time'][step-1]  # sec
            g_force = speed_change / time_change * 1000 / 9.81

            total_force = self.telemetry['F_photon'][step] - self.telemetry['F_gravity'][step]
            gee = total_force / 0.086 / 9.81
            print(step, self.telemetry['stellar_distance'][step], self.telemetry['ship_speed'][step], speed_change, time_change, g_force, gee)


            if self.telemetry['sail_not_parallel'][step] == False:
                has_sail_switched_off = True
                percentage_when_off = self.telemetry['time'][step] / self.telemetry['time'][-1] * 100
                print('Sail switched off at time [seconds]',
                    self.telemetry['time'][step], '{:10.1f}'.format(percentage_when_off), '%')
                break
        if not has_sail_switched_off:
            print('Sail was always on')
   
    
    def get_closest_encounter(self):
        '''Find the time and location of closest approach'''
    
        closest_encounter = float("inf")
        step_of_closest_encounter = 0
     
        for step in range(len(self.telemetry['stellar_distance'])):
            if self.telemetry['stellar_distance'][step] < closest_encounter:
                closest_encounter = self.telemetry['stellar_distance'][step]
                step_of_closest_encounter = step
        
        encounter_time = self.telemetry['time'][step_of_closest_encounter] / 3600  # hours
        
        print('Closest encounter at step', step_of_closest_encounter,
          'with distance [stellar radii]', '{:1.1f}'.format(closest_encounter),
          'at time', '{:1.1f}'.format(encounter_time))

        return encounter_time, step_of_closest_encounter
    
    def get_maximum_total_force(self):
        '''Find the time of maximum force on the sail'''
        
        step_of_maxforce = 0     
        maxforce = 0.0
     
        for step in range(self.nsteps):
            if abs(self.telemetry['F_total'][step] > maxforce):
                maxforce = self.telemetry['F_total'][step]
                step_of_maxforce = step
        
        force_time = self.telemetry['time'][step_of_maxforce] / 3600  # hours
        
        print('Maximum force at step', step_of_maxforce,
          'with magnitude [N]', '{:1.1f}'.format(maxforce),
          'at time', '{:1.1f}'.format(force_time))

        return force_time, step_of_maxforce
        
        
    def write_telemetry_to_file(self,filename,star):
        
        headerstring = 'Recorded telemetry for flight around: \n'
        headerstring = headerstring + str(star)+ '\n'
        headerstring = headerstring + str(self) + '\n'
        headerstring = headerstring + '-------- \n'
        headerstring = headerstring + 'Columns:   step  time  x y z  vx vy vz  F_grav  F_photon  F_mag   F_total  F_grav(x,y,z)  F_photon(x,y,z)  F_magnetic(x,y,z)  sail(x,y,z) sailangle photonacceleration  stellar distance  sail parallel?  \n '
        headerstring = headerstring + '-------- \n'
        
        np.savetxt(filename, self.telemetry, header = headerstring)
        
    def read_from_telemetry_file(self,filename):
        
        # Read header and find sail parameters(mass,area, charge)
        line = linecache.getline(filename, 7)
        header = line.split()
        print header
        self.mass = float(header[3])
        self.area = float(header[5])
        self.charge = float(header[7])
        
        # Read telemetry data and add to sail 
        data = np.genfromtxt(filename,skiprows=15)
        
        self.nsteps = len(data[:,0])
        
        
        self.telemetry['step'] = data[:,0]
        print data[:,0]
        

    