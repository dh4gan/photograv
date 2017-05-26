
import numpy
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from numpy import arctan2, sqrt, pi, sin, cos, radians

import star

from star import sun_radius

def quadratic_limb_darkening(impact, limb1, limb2):
    "Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"
    impact = cos(1 - impact)
    return 1 - limb1 * (1 - impact) - limb2 * (1 - impact) ** 2

def make_figure_flight(
        sail,
        star,
        scale,
        flight_color,
        redness,
        show_burn_circle,
        star_name,
        weight_ratio,
        circle_spacing_minutes,
        annotate_cases,
        caption,
        colorbar=True):
    
    fig = plt.gcf()
    ax = fig.add_subplot(111, aspect='equal')

    # Flight trajectory
    px_stellar_units = sail.telemetry['px'] * sun_radius / star.R
    py_stellar_units = sail.telemetry['py'] * sun_radius / star.R
    
    plt.plot(
        px_stellar_units,
        py_stellar_units,
        color=flight_color,
        linewidth=0.5)

    # G-Forces

    cNorm = matplotlib.colors.Normalize(vmin=-11.5, vmax=2)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='Reds_r')

    last_label = 18.5

    if colorbar:

        for step in range(sail.nsteps):
            # Color bar
            x1 = 15
            x2 = 20
            y1 = sail.telemetry['py'][step] * sun_radius / star.R
            y2 = y1
            if -20 < y1 < 20:
                speed_change = sail.telemetry['ship_speed'][step] - sail.telemetry['ship_speed'][step-1]  # km/sec
                time_change = sail.telemetry['time'][step] - sail.telemetry['time'][step-1]  # sec
                g_force = speed_change / time_change * 1000 / 9.81
                colorVal = scalarMap.to_rgba(g_force)
                plt.plot([x1, x2], [y1, y2], color = colorVal, linewidth = 4.)
            if y1 < last_label + 1 and y1 < 18 and y1 > -18:
                if g_force < -10:
                    textcolor = 'white'
                else:
                    textcolor = 'black'
                text = r'${:2.1f} g$'.format(g_force)
                plt.annotate(
                    text, xy=(17, y1),
                    fontsize=12,
                    horizontalalignment='center',
                    color=textcolor)

                last_label = last_label - 4

    # If desired, show dashed circle for minimum distance
    if show_burn_circle:
        circle2 = plt.Circle(
            (0, 0),
            5,
            color='black',
            fill=False,
            linestyle='dashed',
            linewidth=0.5)  # dashed version
        # shaded version
        # circle2=plt.Circle((0,0), 5, color=(0.9, 0.9, 0.9), fill=True)

        fig.gca().add_artist(circle2)
        

    # Star in the center with limb darkening
    star_quality = 50  # number of shades
    limb1 = 0.4703  # quadratic limb darkening parameters
    limb2 = 0.236
    for i in range(star_quality):
        impact = i / float(star_quality)
        LimbDarkening = quadratic_limb_darkening(impact, limb1, limb2)
        Sun = plt.Circle(
            (0, 0),
            1 - i / float(star_quality),
            color=(LimbDarkening, redness * LimbDarkening, 0))
        plt.gcf().gca().add_artist(Sun)

    # Calculate and print deflection angle (and star name) if desired
    if star_name != '':
        ax.annotate(star_name, xy=(-2.5, 1.75))
        deflection_angle = abs(
            arctan2(sail.telemetry['py'][-1], sail.telemetry['px'][-1]) * (360 / (2 * pi))) - 90

        # print angle
        text_deflection_angle = r'$\delta = {:1.0f} ^\circ$'.format(deflection_angle)

        ax.annotate(
            text_deflection_angle,
            xy=(-scale + 1, scale - 2.5),
            fontsize=16)

        # Add a circle mark at the closest encounter
        index_min = numpy.argmin(sail.telemetry['stellar_distance'])
        min_x_location = sail.telemetry['px'][index_min] * sun_radius / star.R
        min_y_location = sail.telemetry['py'][index_min] * sun_radius / star.R
        marker = plt.Circle(
            (min_x_location, min_y_location), 0.2, color='red', fill=False)
        plt.gcf().gca().add_artist(marker)

        # print weight ratio
        text_weight_ratio = r'$\sigma = {:1.1f}$ g/m$^2$'.format(weight_ratio)
        ax.annotate(text_weight_ratio, xy=(-scale + 1, scale - 5), fontsize=16)

        # print entry speed
        text_entry_speed = r'$\nu_{{\infty}} = {:10.0f}$ km/s'.format(sail.telemetry['ship_speed'][1])
        ax.annotate(text_entry_speed, xy=(-scale + 1, scale - 8), fontsize=16)

        # Add circle marks every [circle_spacing_minutes] and annotate them
        current_marker_number = 0
        it = numpy.nditer(sail.telemetry['time'], flags=['f_index'])
        while not it.finished:
            if (it[0] / 60) % circle_spacing_minutes == 0:
                x_location = sail.telemetry['px'][it.index] * sun_radius / star.R
                y_location = sail.telemetry['py'][it.index] * sun_radius / star.R
                speed = sail.telemetry['ship_speed'][it.index]

                # Check if inside the plotted figure (faster plot generation)
                if abs(x_location) < scale and abs(y_location) < scale:

                    # Yes, mark it
                    marker = plt.Circle(
                        (x_location, y_location),
                        0.1,
                        color='black',
                        fill=True)
                    plt.gcf().gca().add_artist(marker)
                    current_marker_number = current_marker_number + 1

                    # Add sail with angle as marker

                    angle = sail.telemetry['sail_angle'][it.index]
                    if y_location > 0:
                        angle = -angle

                    # Sail is parallel, put it in direction of star:
                    if not sail.telemetry['sail_not_parallel'][it.index]:
                        v1_theta = arctan2(0, 0)
                        v2_theta = arctan2(y_location, x_location)
                        angle = (v2_theta - v1_theta) * (180.0 / pi)
                        angle = radians(angle)

                    # Calculate start and end positions of sail line
                    length = 1
                    endy = y_location + (length * sin(angle))
                    endx = x_location + (length * cos(angle))
                    starty = y_location - (length * sin(angle))
                    startx = x_location - (length * cos(angle))
                    plt.plot([startx, endx], [starty, endy], color='black')

                    # Make arrow between first and second circle mark
                    if current_marker_number == 1:
                        ax.arrow(
                            x_location,
                            y_location - 2, 0.0,
                            -0.01,
                            head_width=0.75,
                            head_length=0.75,
                            fc='k',
                            ec='k',
                            lw=0.01)

                    # To avoid crowding of text labels, set them sparsely
                    if sail.telemetry['ship_speed'][it.index] > 300:
                        n_th_print = 1  # print all labels
                    if 150 < sail.telemetry['ship_speed'][it.index] < 300:
                        n_th_print = 3  # print every 5th label
                    if sail.telemetry['ship_speed'][it.index] < 150:
                        n_th_print = 5  # print every 5th label

                    if -19 < y_location < 19 and (it[0] / 60) % (circle_spacing_minutes * n_th_print) == 0:
                        
                        text_speed = r'{:10.0f} km/s'.format(speed)
                        ax.annotate(
                            text_speed,
                            xy=(x_location + 1, y_location - 0.5),
                            fontsize=12)

            it.iternext()

    # Format the figure
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    minor_ticks = numpy.arange(-scale, scale, 1)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(minor_ticks, minor=True)
    plt.xlim(-scale, scale)
    plt.ylim(-scale, scale)
    plt.xlabel('Distance [Stellar Radii]', fontweight='bold')
    plt.ylabel('Distance [Stellar Radii]', fontweight='bold')

    if caption != '':
        fig.suptitle(
            caption,
            fontsize=16,
            fontweight='bold',
            weight='heavy',
            x=0.22,
            y=0.95)

    if annotate_cases:
        # Annotate cases
        ax.annotate('I, II', xy=(-1, 3))
        ax.annotate('III', xy=(-16, -5))
        ax.annotate('IV', xy=(6, -18))

    return plt


def make_figure_speed(sail, scale, caption):

    encounter_time, step_of_closest_encounter = sail.get_closest_encounter()

    fig = plt.figure(figsize=(6, 6))
    ax = plt.gca()
    minor_xticks = numpy.arange(-scale, scale, scale/10)
    ax.set_xticks(minor_xticks, minor=True)
    minor_yticks = numpy.arange(0, 14000, 400)
    ax.set_yticks(minor_yticks, minor=True)
    ax.get_yaxis().set_tick_params(direction='in')
    ax.get_xaxis().set_tick_params(direction='in')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    fig.suptitle(caption, fontsize=14, fontweight='bold', x=0.131, y=0.95)
    if caption == 'a':
        ax.annotate(r'$\alpha$ Cen A', xy=(2, 12000))
    if caption == 'b':
        ax.annotate(r'$\alpha$ Cen B', xy=(2, 12000))
    if caption == 'c':
        ax.annotate(r'$\alpha$ Cen C', xy=(2, 12000))

    time = sail.telemetry['time'] / 3600 - encounter_time
    speed = sail.telemetry['ship_speed']
    plt.plot(time, speed, color='black', linewidth=0.5)

    # Vertical line
    plt.plot(
        [0, 0],
        [0, 14000],
        linestyle='dotted',
        color='black',
        linewidth=0.5)

    plt.xlim(-scale, +scale)
    plt.ylim(0, 14000)
    plt.xlabel('Time Around Closest Approach [Hours]', fontweight='bold')
    plt.ylabel('Speed [km/s]', fontweight='bold')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    return plt





def make_figure_forces(sail, star,scale):

    # Find closest encounter
    encounter_time, step_of_closest_encounter = sail.get_closest_encounter()

    # select sail.telemetry from start to closest encounter
    px_stellar_units = sail.telemetry['px'] * sun_radius / star.R
    py_stellar_units = sail.telemetry['py'] * sun_radius / star.R

    # distance between sail and origin in [m]
    distance = sqrt((px_stellar_units**2) + (py_stellar_units**2))
    distance[:step_of_closest_encounter] = -distance[:step_of_closest_encounter]

    total_force = sail.telemetry['F_photon'] - sail.telemetry['F_gravity']
    gee = total_force / 0.086 / 9.81
    print('max total_force', numpy.max(total_force))
    print('max gee', numpy.max(gee))

    # Figure
    fig = plt.figure(figsize=(6, 6))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    ax = plt.gca()

    minor_xticks = numpy.arange(-scale, scale, 5)
    ax.set_xticks(minor_xticks, minor=True)
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    fig.suptitle('a', fontsize=14, fontweight='bold', x=0.131, y=0.95)
    time = sail.telemetry['time'] / 3600 - encounter_time
    total_force = sail.telemetry['F_photon'][:step_of_closest_encounter]
    - sail.telemetry['F_gravity'][:step_of_closest_encounter]
    ax.plot(distance[:step_of_closest_encounter],
        sail.telemetry['F_gravity'][:step_of_closest_encounter],
        color='black', linewidth=1.)
    ax.plot(distance[step_of_closest_encounter:],
        sail.telemetry['F_gravity'][step_of_closest_encounter:],
        color='black', linewidth=1.)
    ax.plot(distance[:step_of_closest_encounter],
        sail.telemetry['F_photon'][:step_of_closest_encounter],
        color='black', linestyle = 'dashed', linewidth=1.)
    offset = 0.02  # to show photon and total lines separatley
    ax.plot(distance[:step_of_closest_encounter],
        total_force-(offset*total_force), color='red', linewidth=1.)

    # Vertical line
    ax.plot(
        [0, 0],
        [0, numpy.amax(total_force * 1.1)],
        linestyle='dotted',
        color='black',
        linewidth=0.5)

    ax.annotate('Photon', xy=(-90, 12), rotation=20, color='black')
    ax.annotate('Total', xy=(-90, 4), rotation=20, color='red')
    ax.annotate('Gravity', xy=(-90, 2.1e-3), rotation=20, color='black')
    ax.yaxis.tick_left()
    ax.set_xlabel('Stellar Distance [Stellar Radii]', fontweight='bold')
    ax.set_ylabel('Force Acting on the Sail [Newton]', fontweight='bold')
    ax.set_xlim(-scale, + scale)
    ax.set_ylim(10e-4, numpy.amax(total_force * 1.2))
    ax.set_yscale('log')

    # Right axis: g-force
    ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.set_xlim(-scale, + scale)
    ax2.set_ylim(5, numpy.max(gee) * 1.2)
    ax2.set_yscale('log')
    ax2.set_ylabel(
        r'Deceleration of the Sail [$g=9.81$m/s$^2$]', fontweight='bold')
    plt.yscale('log')

    # Inset with zoom
    ax1 = fig.add_subplot(111)
    axins = inset_axes(
        ax,
        1.75,
        1.75,
        loc=2,
        bbox_to_anchor=(0.56, 0.89),
        bbox_transform=ax.figure.transFigure)  # no zoom
    offset = 0.001  # to show photon and total lines separatley

    axins.plot(
        distance[:step_of_closest_encounter],
        sail.telemetry['F_photon'][:step_of_closest_encounter],
        color='black',
        linewidth=1,
        linestyle='dashed')
    axins.plot(
        distance[:step_of_closest_encounter],
        total_force - (offset * total_force),
        color='red',
        linewidth=1)

    # Vertical line
    axins.plot(
        [0, 0],
        [0, numpy.amax(total_force * 1.1)],
        linestyle='dotted',
        color='black',
        linewidth=0.5)
    axins.set_xlim(-5.5, -5)
    axins.set_ylim(1200, 1400)

    return plt


def make_figure_distance_pitch_angle(sail, scale, star):

    encounter_time, step_of_closest_encounter = sail.get_closest_encounter()

    # select sail.telemetry from start to closest encounter
    stellar_radii = sun_radius / star.R
    px_stellar_units = sail.telemetry['px'][:step_of_closest_encounter] * stellar_radii
    py_stellar_units = sail.telemetry['py'][:step_of_closest_encounter] * stellar_radii
    pitch_angle = sail.telemetry['alpha'][:step_of_closest_encounter]
    pitch_angle = -numpy.degrees(pitch_angle)

    # distance between sail and origin in [m]
    distance = sqrt((px_stellar_units**2) + (py_stellar_units**2))

    # make figure
    fig = plt.figure(figsize=(6, 6))
    ax = plt.gca()
    fig.suptitle('b', fontsize=14, fontweight='bold', x=0.131, y=0.95)
    plt.plot(distance, pitch_angle, linewidth=0.5, color='black')
    plt.xscale('log')
    ax.set_xticks([5, 10, 100])
    ax.set_xticklabels(["50", "10", "100"])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim(4, 100)
    plt.ylim(-45, 0)
    plt.xlabel('Stellar Distance [Stellar Radii]', fontweight='bold')
    plt.ylabel(
        r'Sail Pitch Angle $\alpha$ During Approach [Degrees]',
        fontweight='bold')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    return plt



def make_figure_multiple_sails(
    sailarray,
    star,
    minimum_distance_from_star,
    afterburner_distance,
    timestep,
    return_mission,
    scale,
    caption,
    colorbar=False):
    
    '''Creates a figure with multiple flights, each with a different sail'''

    color = 'black'
    redness = 0.7
    star_name = ''
    circle_spacing_minutes = 2e10  # never
    show_burn_circle = True
    
    annotate_cases = False  # Avoid printing it many times

    niterations = len(sailarray)
    iter_counter = 0
    
    for ship in sailarray:
        iter_counter += 1
        weight_ratio = 1/ship.area
        
        print(
            'Now running iteration', iter_counter,
            'of', int(niterations))
        
        print ship
        ship.fly(star,minimum_distance_from_star,afterburner_distance,timestep,return_mission)
        
        
        color_shade = float(iter_counter) / float(niterations)
        
        if color_shade > 1:
            color_shade = 1
        if color_shade < 0:
            color_shade = 0
        
        color = [1 - color_shade, 0, color_shade]
        print color

        my_figure = make_figure_flight(
            ship,
            star,
            scale,
            color,
            redness,
            show_burn_circle,
            star_name,
            weight_ratio,
            circle_spacing_minutes,
            annotate_cases,
            caption='',
            colorbar=False)

    fig = plt.gcf()
    ax = fig.add_subplot(111, aspect='equal')
    ax.annotate('I, II', xy=(-1, 3))
    ax.annotate('III', xy=(-16, -5))
    ax.annotate('IV', xy=(6, -18))
    fig.suptitle('a', fontsize=16, fontweight='bold', weight='heavy', x=0.22, y=0.95)

    return my_figure

def make_figure_multiple_stars(
    sail1,
    stararray,
    minimum_distance_from_star,
    afterburner_distance,
    timestep,
    return_mission,
    scale,
    caption,
    colorbar=False):
    
    '''Creates a figure with multiple flights, each with a different star'''

    color = 'black'
    redness = 0.7
    star_name = ''
    circle_spacing_minutes = 2e10  # never
    weight_ratio = 1/sail1.area
    show_burn_circle = True
    
    annotate_cases = False  # Avoid printing it many times

    niterations = len(stararray)
    iter_counter = 0
    
    for s in stararray:
        iter_counter += 1
        
        
        print(
            'Now running iteration', iter_counter,
            'of', int(niterations))
        print s
        sail = sail1.clone()
        print sail
        
        sail.fly(s,minimum_distance_from_star,afterburner_distance,timestep,return_mission =False)
        #sail.print_flight_report()
        encounter_time, step_of_closest_encounter = sail.get_closest_encounter()
        print encounter_time, step_of_closest_encounter
        color_shade = float(iter_counter) / float(niterations)
        
        if color_shade > 1:
            color_shade = 1
        if color_shade < 0:
            color_shade = 0
        
        color = [1 - color_shade, 0, color_shade]

        my_figure = make_figure_flight(
            sail,
            s,
            scale,
            color,
            redness,
            show_burn_circle,
            star_name,
            weight_ratio,
            circle_spacing_minutes,
            annotate_cases,
            caption='',
            colorbar=False)

    fig = plt.gcf()
    ax = fig.add_subplot(111, aspect='equal')
    ax.annotate('I, II', xy=(-1, 3))
    ax.annotate('III', xy=(-16, -5))
    ax.annotate('IV', xy=(6, -18))
    fig.suptitle('a', fontsize=16, fontweight='bold', weight='heavy', x=0.22, y=0.95)

    return my_figure


def make_video_flight(
        sail,
        star,
        scale,
        flight_color,
        redness,
        show_burn_circle,
        star_name,
        weight_ratio,
        circle_spacing_minutes,
        annotate_cases,
        caption):
    fig = plt.gcf()
    ax = fig.add_subplot(111, aspect='equal')

    # Flight trajectorie
    px_stellar_units = sail.telemetry['px'] * sun_radius / star.R
    py_stellar_units = sail.telemetry['py'] * sun_radius / star.R

    # If desired, show dashed circle for minimum distance
    if show_burn_circle:
        circle2 = plt.Circle(
            (0, 0),
            5,
            color='black',
            fill=False,
            linestyle='dashed',
            linewidth=0.5)  # dashed version
        # shaded version
        # circle2=plt.Circle((0,0), 5, color=(0.9, 0.9, 0.9), fill=True)
        fig.gca().add_artist(circle2)

    # Star in the center with limb darkening
    star_quality = 50  # number of shades
    limb1 = 0.4703  # quadratic limb darkening parameters
    limb2 = 0.236
    time_start = 0
    cNorm = matplotlib.colors.Normalize(vmin=-11.5, vmax=2)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='Reds_r')

    # Calculate and print deflection angle (and star name) if desired
    it = numpy.nditer(sail.telemetry['time'], flags=['f_index'])
    while not it.finished:

        x_location = sail.telemetry['px'][it.index] * sun_radius / star.R
        y_location = sail.telemetry['py'][it.index] * sun_radius / star.R
        speed = sail.telemetry['ship_speed'][it.index]
        time = sail.telemetry['time'][it.index]

        # Check if inside the plotted figure (faster plot generation)
        if abs(x_location) < scale and abs(y_location) < scale + 2:
            if time_start == 0:
                time_start = sail.telemetry['time'][it.index]
            plt.cla()

            # Colorbar
            last_label = 18.5
            for step in range(it.index):
                # Color bar
                x1 = 15
                x2 = 20
                y1 = sail.telemetry['py'][step] * sun_radius / star.R
                y2 = y1
                if -20 < y1 < 20:
                    speed_change = sail.telemetry['ship_speed'][step] - sail.telemetry['ship_speed'][step-1]  # km/sec
                    time_change = sail.telemetry['time'][step] - sail.telemetry['time'][step-1]  # sec
                    g_force = speed_change / time_change * 1000 / 9.81
                    colorVal = scalarMap.to_rgba(g_force)
                    plt.plot([x1, x2], [y1, y2], color=colorVal, linewidth=4)
                if y1 < last_label + 1 and y1 < 18 and y1 > -18:
                    if g_force < -10:
                        textcolor = 'white'
                    else:
                        textcolor = 'black'
                    text = r'${:2.1f} g$'.format(g_force)
                    plt.annotate(
                        text,
                        xy=(17, y1),
                        fontsize=12,
                        horizontalalignment='center',
                        color=textcolor)

                    last_label = last_label - 4

            # Flight line
            plt.plot(
                px_stellar_units[:it.index + 1],
                py_stellar_units[:it.index + 1],
                color=flight_color,
                linewidth=0.5)

            # Star
            for i in range(star_quality):
                impact = i / float(star_quality)
                LimbDarkening = quadratic_limb_darkening(impact, limb1, limb2)
                Sun = plt.Circle(
                    (0, 0),
                    1 - i / float(star_quality),
                    color=(LimbDarkening, redness * LimbDarkening, 0))
                plt.gcf().gca().add_artist(Sun)

            ax.annotate(star_name, xy=(-2.5, 1.75))
            deflection_angle = abs(arctan2(sail.telemetry['py'][-1], sail.telemetry['px'][-1]) * (360 / (2 * pi))) - 90

            # print angle
            text_deflection_angle = r'$\delta = {:1.0f} ^\circ$'.format(deflection_angle)
            ax.annotate(
                text_deflection_angle,
                xy=(-scale + 1, scale - 2.5),
                fontsize=16)

            # Add a circle mark at the closest encounter
            index_min = numpy.argmin(sail.telemetry['stellar_distance'])
            if it.index > index_min:  # past the closest encounter
                stellar_radii = sun_radius / star.R
                min_x_location = sail.telemetry['px'][index_min] * stellar_radii
                min_y_location = sail.telemetry['py'][index_min] * stellar_radii
                # print('drawing closest encounter at (x,y)=', min_x_location, min_y_location)
                marker = plt.Circle(
                    (min_x_location, min_y_location),
                    0.2,
                    color='red',
                    fill=False)
                plt.gcf().gca().add_artist(marker)

            # print weight ratio
            text_weight_ratio = r'$\sigma = {:1.1f}$ g/m$^2$'.format(weight_ratio)
            ax.annotate(
                text_weight_ratio,
                xy=(-scale + 1, scale - 5),
                fontsize=16)

            # print entry speed
            text_entry_speed = r'$\nu_{{\infty}} = {:10.0f}$ km/s'.format(sail.telemetry['ship_speed'][1])
            ax.annotate(
                text_entry_speed,
                xy=(-scale + 1, scale - 8),
                fontsize=16)

            # Add circle marks every [circle_spacing_minutes] and annotate them
            current_marker_number = 0

            # Yes, mark it
            marker = plt.Circle(
                (x_location, y_location), 0.1, color='black', fill=True)
            plt.gcf().gca().add_artist(marker)
            current_marker_number = current_marker_number + 1

            # Show sail angle as line
            angle = sail.telemetry['sail_angle'][it.index]
            if y_location > 0:
                angle = -angle

            # Sail is parallel, put it in direction of star:
            if not sail.telemetry['sail_not_parallel'][it.index]:
                v1_theta = arctan2(0, 0)
                v2_theta = arctan2(y_location, x_location)
                angle = (v2_theta - v1_theta) * (180.0 / 3.141)
                angle = radians(angle)

            # Calculate start and end positions of sail line
            length = 1
            endy = y_location + (length * sin(angle))
            endx = x_location + (length * cos(angle))
            starty = y_location - (length * sin(angle))
            startx = x_location - (length * cos(angle))
            plt.plot([startx, endx], [starty, endy], color='black')

            """
            # Add circle marks every [circle_spacing_minutes] and annotate them
            current_marker_number = 0
            inner_loop = numpy.nditer(sail.telemetry['time'], flags=['f_index'])
            while not inner_loop.finished and inner_loop.index < it.index:
                if (inner_loop[0] / 60) % circle_spacing_minutes == 0:
                    x_location = sail.telemetry['px'][inner_loop.index] * sun_radius / star.R
                    y_location = sail.telemetry['py'][inner_loop.index] * sun_radius / star.R
                    speed = sail.telemetry['ship_speed'][inner_loop.index]

                    # Check if inside the plotted figure (faster plot generation)
                    if abs(x_location) < scale and abs(y_location) < scale:

                        # Yes, mark it
                        marker = plt.Circle(
                            (x_location, y_location),
                            0.1,
                            color='black',
                            fill=True)
                        plt.gcf().gca().add_artist(marker)
                        current_marker_number = current_marker_number + 1

                        # To avoid crowding of text labels, set them sparsely
                        if sail.telemetry['ship_speed'][inner_loop.index] > 300:
                            n_th_print = 1  # print all labels
                        if 150 < sail.telemetry['ship_speed'][inner_loop.index] < 300:
                            n_th_print = 3  # print every 5th label
                        if sail.telemetry['ship_speed'][inner_loop.index] < 150:
                            n_th_print = 5  # print every 5th label

                        if -19 < y_location < 19 and (inner_loop[0] / 60) % (circle_spacing_minutes * n_th_print) == 0:
                            test_speed = ''
                            text_speed = r'{:10.0f} km/s'.format(speed)
                            ax.annotate(
                                text_speed,
                                xy=(x_location + 1, y_location - 0.5),
                                fontsize=12)
                inner_loop.iternext()
            """

            # Annotate current speed and time in lower left corner
            text_speed = r'Speed: {:10.0f} km/s'.format(speed)
            text_time = r'Time: {:10.0f} hrs'.format((time - time_start) / 60 / 60)
            ax.annotate(text_speed, xy=(-19, -17), fontsize=12)
            ax.annotate(text_time, xy=(-19, -19), fontsize=12)

            # Format the figure
            #plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            minor_ticks = numpy.arange(-scale, scale, 1)
            ax.set_xticks(minor_ticks, minor=True)
            ax.set_yticks(minor_ticks, minor=True)
            plt.xlim(-scale, scale)
            plt.ylim(-scale, scale)
            plt.xlabel('Distance [Stellar Radii]', fontweight='bold')
            plt.ylabel('Distance [Stellar Radii]', fontweight='bold')
            fig.savefig(str(it.index) + ".png", bbox_inches='tight', dpi=400)
            print('it.index', it.index)

        it.iternext()

    return plt




def make_video_flight_black(
        sail,
        star,
        scale,
        flight_color,
        redness,
        show_burn_circle,
        star_name,
        weight_ratio,
        circle_spacing_minutes,
        annotate_cases,
        caption):
    fig = plt.gcf()
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_bgcolor('none')
    ax.set_axis_bgcolor((0, 0, 0))
    ax.patch.set_facecolor('black')
    fig.patch.set_facecolor('black')
    plt.axis('off')

    # Flight trajectorie
    px_stellar_units = sail.telemetry['px'] * sun_radius / star.R
    py_stellar_units = sail.telemetry['py'] * sun_radius / star.R

    # If desired, show dashed circle for minimum distance
    if show_burn_circle:
        circle2 = plt.Circle(
            (0, 0),
            5,
            color='black',
            fill=False,
            linestyle='dashed',
            linewidth=0.5)  # dashed version
        # shaded version
        # circle2=plt.Circle((0,0), 5, color=(0.9, 0.9, 0.9), fill=True)
        fig.gca().add_artist(circle2)

    # Star in the center with limb darkening
    star_quality = 50  # number of shades
    limb1 = 0.4703  # quadratic limb darkening parameters
    limb2 = 0.236
    time_start = 0
    cNorm = matplotlib.colors.Normalize(vmin=-11.5, vmax=2)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap='Reds_r')

    # Calculate and print deflection angle (and star name) if desired
    it = numpy.nditer(sail.telemetry['time'], flags=['f_index'])
    while not it.finished:
        """
        fig = plt.gcf()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_axis_bgcolor('none')
        ax.set_axis_bgcolor((0, 0, 0))
        ax.patch.set_facecolor('black')
        fig.patch.set_facecolor('black')
        plt.axis('off')
        ax.set_axis_off()
        """
        x_location = sail.telemetry['px'][it.index] * sun_radius / star.R
        y_location = sail.telemetry['py'][it.index] * sun_radius / star.R
        speed = sail.telemetry['ship_speed'][it.index]
        time = sail.telemetry['time'][it.index]

        # Check if inside the plotted figure (faster plot generation)
        if abs(x_location) < scale and abs(y_location) < scale + 2:
            if time_start == 0:
                time_start = sail.telemetry['time'][it.index]
            plt.cla()

            # Colorbar
            last_label = 18.5
            for step in range(it.index):
                # Color bar
                x1 = 15
                x2 = 20
                y1 = sail.telemetry['py'][step] * sun_radius / star.R
                y2 = y1
                if -20 < y1 < 20:
                    speed_change = sail.telemetry['ship_speed'][step] - sail.telemetry['ship_speed'][step-1]  # km/sec
                    time_change = sail.telemetry['time'][step] - sail.telemetry['time'][step-1]  # sec
                    g_force = speed_change / time_change * 1000 / 9.81
                    colorVal = scalarMap.to_rgba(g_force)
                    plt.plot([x1, x2], [y1, y2], color=colorVal, linewidth=4)
                if y1 < last_label + 1 and y1 < 18 and y1 > -18:
                    if g_force < -10:
                        textcolor = 'white'
                    else:
                        textcolor = 'black'
                    text = r'${:2.1f} g$'.format(g_force)
                    plt.annotate(
                        text,
                        xy=(17, y1),
                        fontsize=12,
                        horizontalalignment='center',
                        color=textcolor)

                    last_label = last_label - 4

            # Flight line
            plt.plot(
                px_stellar_units[:it.index + 1],
                py_stellar_units[:it.index + 1],
                color=flight_color,
                linewidth=0.5)

            # Star
            for i in range(star_quality):
                impact = i / float(star_quality)
                LimbDarkening = quadratic_limb_darkening(impact, limb1, limb2)
                Sun = plt.Circle(
                    (0, 0),
                    1 - i / float(star_quality),
                    color=(LimbDarkening, redness * LimbDarkening, 0))
                plt.gcf().gca().add_artist(Sun)

            ax.annotate(star_name, xy=(-2.5, 1.75))
            deflection_angle = abs(arctan2(sail.telemetry['py'][-1], sail.telemetry['px'][-1]) * (360 / (2 * pi))) - 90

            # print angle
            text_deflection_angle = r'$\delta = {:1.0f} ^\circ$'.format(deflection_angle)
            ax.annotate(
                text_deflection_angle,
                xy=(-scale + 1, scale - 2.5),
                fontsize=16,
                color='white')

            # Add a circle mark at the closest encounter
            index_min = numpy.argmin(sail.telemetry['stellar_distance'])
            if it.index > index_min:  # past the closest encounter
                stellar_radii = sun_radius / star.R
                min_x_location = sail.telemetry['px'][index_min] * stellar_radii
                min_y_location = sail.telemetry['py'][index_min] * stellar_radii
                # print('drawing closest encounter at (x,y)=', min_x_location, min_y_location)
                marker = plt.Circle(
                    (min_x_location, min_y_location),
                    0.2,
                    color='red',
                    fill=False)
                plt.gcf().gca().add_artist(marker)

            # print weight ratio
            text_weight_ratio = r'$\sigma = {:1.1f}$ g/m$^2$'.format(weight_ratio)
            ax.annotate(
                text_weight_ratio,
                xy=(-scale + 1, scale - 5),
                fontsize=16,
                color='white')

            # print entry speed
            text_entry_speed = r'$\nu_{{\infty}} = {:10.0f}$ km/s'.format(sail.telemetry['ship_speed'][1])
            ax.annotate(
                text_entry_speed,
                xy=(-scale + 1, scale - 8),
                fontsize=16,
                color='white')

            """
            # Add circle marks every [circle_spacing_minutes] and annotate them
            current_marker_number = 0

            # Yes, mark it
            marker = plt.Circle(
                (x_location, y_location), 0.1, color='white', fill=True)
            plt.gcf().gca().add_artist(marker)
            current_marker_number = current_marker_number + 1
            """

            # Show sail angle as line
            angle = sail.telemetry['alpha'][it.index]
            if y_location > 0:
                angle = -angle

            # Sail is parallel, put it in direction of star:
            if not sail.telemetry['sail_not_parallel'][it.index]:
                v1_theta = arctan2(0, 0)
                v2_theta = arctan2(y_location, x_location)
                angle = (v2_theta - v1_theta) * (180.0 / 3.141)
                angle = radians(angle)

            # Calculate start and end positions of sail line
            length = 1
            endy = y_location + (length * sin(angle))
            endx = x_location + (length * cos(angle))
            starty = y_location - (length * sin(angle))
            startx = x_location - (length * cos(angle))
            plt.plot([startx, endx], [starty, endy], color='white')

            """
            # Add circle marks every [circle_spacing_minutes] and annotate them
            current_marker_number = 0
            inner_loop = numpy.nditer(sail.telemetry['time'], flags=['f_index'])
            while not inner_loop.finished and inner_loop.index < it.index:
                if (inner_loop[0] / 60) % circle_spacing_minutes == 0:
                    x_location = sail.telemetry['px'][inner_loop.index] * sun_radius / star.R
                    y_location = sail.telemetry['py'][inner_loop.index] * sun_radius / star.R
                    speed = sail.telemetry['ship_speed'][inner_loop.index]

                    # Check if inside the plotted figure (faster plot generation)
                    if abs(x_location) < scale and abs(y_location) < scale:

                        # Yes, mark it
                        marker = plt.Circle(
                            (x_location, y_location),
                            0.1,
                            color='black',
                            fill=True)
                        plt.gcf().gca().add_artist(marker)
                        current_marker_number = current_marker_number + 1

                        # To avoid crowding of text labels, set them sparsely
                        if sail.telemetry['ship_speed'][inner_loop.index] > 300:
                            n_th_print = 1  # print all labels
                        if 150 < sail.telemetry['ship_speed'][inner_loop.index] < 300:
                            n_th_print = 3  # print every 5th label
                        if sail.telemetry['ship_speed'][inner_loop.index] < 150:
                            n_th_print = 5  # print every 5th label

                        if -19 < y_location < 19 and (inner_loop[0] / 60) % (circle_spacing_minutes * n_th_print) == 0:
                            test_speed = ''
                            text_speed = r'{:10.0f} km/s'.format(speed)
                            ax.annotate(
                                text_speed,
                                xy=(x_location + 1, y_location - 0.5),
                                fontsize=12)
                inner_loop.iternext()
            """

            # Annotate current speed and time in lower left corner
            text_speed = r'Speed: {:10.0f} km/s'.format(speed)
            text_time = r'Time: {:10.0f} hrs'.format((time - time_start) / 60 / 60)
            ax.annotate(text_speed, xy=(-19, -17), fontsize=12, color='white')
            ax.annotate(text_time, xy=(-19, -19), fontsize=12, color='white')

            # Format the figure
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            #minor_ticks = numpy.arange(-scale, scale, 1)
            #ax.set_xticks(minor_ticks, minor=True)
            #ax.set_yticks(minor_ticks, minor=True)
            plt.xlim(-scale, scale)
            plt.ylim(-scale, scale)
            #plt.xlabel('Distance [Stellar Radii]', fontweight='bold')
            #plt.ylabel('Distance [Stellar Radii]', fontweight='bold')
            ax.set_axis_bgcolor('none')
            #ax.set_axis_bgcolor((0, 0, 0))

            plt.axis('off')
            ax.set_axis_off()
            #fig.axes.get_xaxis().set_visible(False)
            #fig.axes.get_yaxis().set_visible(False)
            ax.patch.set_facecolor('black')
            fig.patch.set_facecolor('black')
            #circ=
            ax.add_patch(plt.Circle((0,0), radius=100, color='black', fill=True))
            fig.savefig(str(it.index) + ".png", bbox_inches='tight', dpi=400, pad_inches = 0,  facecolor=fig.get_facecolor())
            print('it.index', it.index)

        it.iternext()

    return plt
