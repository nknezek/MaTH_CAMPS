"""
build_planet.py

Framework to self-consistently calculate pressure, temperature, and density, among other
physical parameters. Takes in a list of homogenous layers and iterates to a consistent
solution. Built using the burnman framework: http://www.burnman.org/.
Allows insertion of heFESTo values for mantle.

Created by Nicholas Knezek based on code by Ian Rose
CIDER July 2014
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import RectBivariateSpline
import re
import csv

import burnman
import burnman.minerals as minerals
import burnman.composite as composite
import burnman.geotherm as gt
import solidus

# constants
G = 6.67e-11


class Planet:
    def __init__(self, boundaries, compositions, T_bounds, methods=None, heFESTo_filename=None, masses=None):
        """
        args:
            boundaries: list of increasing outer boundaries of layers (core to surface)
            compositions: list of burnman.Composite of burnman.Material
            methods: list of EOS fitting method
            T_bounds: list of temperatures for boundaries
            heFESTo_filename: csv file containing pressure, temperature, and density values
            masses: list containing mass of each layer
        """
        for c in compositions:
            assert ( isinstance(c, burnman.Material) )

        for i, b in enumerate(boundaries):
            assert ( i == 0 or b >= boundaries[i - 1])

        assert ( len(boundaries) == len(compositions) )

        assert (T_bounds[-1] > 0)

        self.boundaries = boundaries
        self.compositions = compositions
        self.Nlayer = len(boundaries)
        self.T_bounds = T_bounds
        self.E_internal = np.empty_like(boundaries)

        if masses is None:
            self.masses = self.compute_masses
        else:
            self.masses = masses

        if methods is None:
            meths = ['slb3'] * self.Nlayer
        else:
            meths = methods

        for m, comp in zip(meths, self.compositions):
            comp.set_method(m)

        # Calculate Freezing Point for Fe/FeS mixture
        self.solidus_solver = solidus.Solver()

        # Import heFESTo mantle values for P,rho,T
        if not (heFESTo_filename == None):
            self.inv_pressure_mantle_func, self.pressure_mantle_func, self.temperature_mantle_func, self.density_mantle_func = self.import_heFESTo_mantle(
                heFESTo_filename)

    def import_heFESTo_mantle(self, inputfilename):
        '''
        Imports heFESTo pressure, temperature, and density values for use in mantle
        '''
        density_mantle = []
        pressure_mantle = []
        temperature_mantle = []

        with open(inputfilename, 'rU') as f:
            lines = (line.rstrip() for line in f)  # All lines including the blank ones, ignore header
            lines = (line for line in lines if line)  # Non-blank lines
            next(lines)  # Ignore header
            for line in lines:
                line = re.split(',', line)
                pressure_mantle.append(float(line[0]))
                temperature_mantle.append(float(line[1]))
                density_mantle.append(float(line[-1]))
        assert (len(density_mantle) == len(pressure_mantle) and len(pressure_mantle) == len(temperature_mantle))

        pressure_mantle = np.array(pressure_mantle) * 1.e9
        density_mantle = np.array(density_mantle) * 1.e3
        temperature_mantle = np.array(temperature_mantle)
        N = len(pressure_mantle)
        xx = np.linspace(0., 1., N)

        inv_pressure_mantle_func = interp1d(pressure_mantle, xx)
        pressure_mantle_func = interp1d(xx, pressure_mantle)
        temperature_mantle_func = interp1d(xx, temperature_mantle)
        density_mantle_func = interp1d(xx, density_mantle)

        return inv_pressure_mantle_func, pressure_mantle_func, temperature_mantle_func, density_mantle_func

    def set_boundary_temperatures(self, T_bounds):
        '''
        Sets the boundary temperatures
        args:
            T_bounds: array of temperatures at the top of each boundary [K]
        '''
        self.T_bounds = T_bounds
        assert (self.T_bounds[-1] > 0)

    def set_boundaries(self, R_bounds):
        '''
        Sets the radius of the boundaries between layers
        args:
            R_bounds: array of radii of boundaries [m]
        '''
        for i, b in enumerate(R_bounds):
            assert ( i == 0 or b >= R_bounds[i - 1])
        self.boundaries = R_bounds

    def set_compositions(self, compositions):
        '''
        Sets the composition of each boundary
        args:
            compositions: array of compositions, matching the number of boundaries [burnman.mineral]
        '''
        for c in compositions:
            assert ( isinstance(c, burnman.Material) )
        self.compositions = compositions

    def set_masses(self, masses):
        ''' Sets the mass of each layer
            args:
                masses: array of masses matching number of boundaries
        '''
        self.masses = masses

    def evaluate_eos(self, pressure, temperature, radii):
        '''
        Find densities for a given set of pressure and temperature
        args:
            pressure: array of starting pressure in Pa
            temperature: array of constant temperature in K
            radii: array of radii in m
        '''
        assert (radii.max() <= self.boundaries[-1] and radii.min() >= 0. )

        densities = np.empty_like(radii)
        vp = np.empty_like(radii)
        vs = np.empty_like(radii)
        vphi = np.empty_like(radii)
        K = np.empty_like(radii)
        G = np.empty_like(radii)
        # iterate over layers
        last = -1.
        for bound, comp, T in zip(self.boundaries, self.compositions, self.T_bounds):
            layer = (radii > last) & ( radii <= bound)
            rrange = radii[layer]
            drange = np.empty_like(rrange)
            vptemp = np.empty_like(rrange)
            vstemp = np.empty_like(rrange)
            vphitemp = np.empty_like(rrange)
            Ktemp = np.empty_like(rrange)
            Gtemp = np.empty_like(rrange)
            prange = pressure[layer]
            trange = temperature[layer]
            for i in range(len(rrange)):
                density, vpi, vsi, vphii, Ki, Gi = burnman.velocities_from_rock(comp, np.array([prange[i]]),
                                                                                np.array([trange[i]]))
                drange[i] = density
                vptemp[i] = vpi
                vstemp[i] = vsi
                vphitemp[i] = vphii
                Ktemp[i] = Ki
                Gtemp[i] = Gi
            densities[layer] = drange
            vp[layer] = vptemp
            vs[layer] = vstemp
            vphi[layer] = vphitemp
            K[layer] = Ktemp
            G[layer] = Gtemp
            last = bound  # update last boundary
        physical_params = {'vp': vp, 'vs': vs, 'vphi': vphi, 'K': K, 'G': G}
        return densities, physical_params

    def compute_adiabat(self, pressure, radii):
        '''
        Find temperature adiabat for a given set of pressures and radii and self.T_bounds
        args:
            pressure: array of pressures in Pa
            radii: array of radii in m
        '''
        Ta = np.empty_like(pressure)
        bound = list(self.boundaries)
        bound.insert(0, -1.0)
        T_bottom_of_previous = self.T_bounds[-1]
        for i, (comp, T_in) in reversed(list(enumerate(zip(self.compositions, self.T_bounds)))):
            if (T_in < 0):
                T = T_bottom_of_previous
            elif (T_in > 0):
                T = T_in
            layer = (radii > bound[i]) & ( radii <= bound[i + 1])
            rrange = radii[layer]
            prange = pressure[layer]
            temp = gt.adiabatic(prange[::-1], T, comp)
            Ta[layer] = temp[::-1]
            T_bottom_of_previous = temp[-1]  # store value of temp at bottom of boundary
        return Ta

    def compute_gravity(self, density, radii):
        '''
        Find gravity given density and radii
        args:
            density: array of densities [kg/m^3]
            radii: array of radii [m]
        '''
        rhofunc = UnivariateSpline(radii, density)
        poisson = lambda p, x: 4.0 * np.pi * G * rhofunc(x) * x * x
        grav = np.ravel(odeint(poisson, 0.0, radii))
        grav[1:] = grav[1:] / radii[1:] / radii[1:]
        grav[0] = 0.0
        return grav

    def compute_pressure(self, density, gravity, radii):
        '''
        Find pressure given density, gravity, and radii
        args:
            density: array of densities [kg/m^3]
            gravity: array of local gravities [m/s^2]
            radii: array of radii [m]
        '''
        depth = radii[-1] - radii
        rhofunc = UnivariateSpline(depth[::-1], density[::-1])
        gfunc = UnivariateSpline(depth[::-1], gravity[::-1])
        pressure = np.ravel(odeint((lambda p, x: gfunc(x) * rhofunc(x)), 0.0, depth[::-1]))
        return pressure[::-1]

    def compute_masses(self, radii, density):
        """
        computes mass of each layer
        :param radii: array of radius values
        :param density: array of density values
        :return: masses array for each layer
        """
        dr = radii[1]-radii[0]
        bounds = list(self.boundaries)
        bounds.insert(0,0)
        masses_tmp = []
        for i in range(len(self.boundaries)):
            layer = (radii > bounds[i]) & ( radii <= bounds[i + 1])
            rrange = radii[layer]
            masses_tmp.append(4 * np.pi * dr * sum(rrange*rrange*density[layer]))
        return masses_tmp

    def compute_planet_mass(self, radius, density):
        '''
        calculates the total planet mass.
        args:
            radius: array of radii
            density: array of densities
        '''
        dr = radius[1] - radius[0]
        return 4 * np.pi * dr * sum(radius * radius * density)

    def compute_thermal_energy(self, radii, density, temperature):
        '''
        Finds the internal energy due to temperature based on stated heat capacity
        of compositions
        returns thermal energy for each layer
        if 'Cv' is not specified in compositions parameters, Debye model is used to
        calculate
        args:
            density: array of densities [kg/m^3]
            temperature: array of temperatures [K]
            radii: array of radii [m]
        '''
        dr = radii[1] - radii[0]
        last = -1
        Eth_out = []
        for i, (bound, comp) in enumerate(zip(self.boundaries, self.compositions)):
            if 'Cv' in comp.params:
                Cv = comp.params['Cv']
            else:
                Cv = comp.heat_capacity_v()
            layer = (radii > last) & ( radii <= bound)
            rrange = radii[layer]
            drange = density[layer]
            trange = temperature[layer]
            Eth_out.append(Cv * 4.0 * np.pi * dr * sum(drange * trange * rrange * rrange))
            last = bound  # update last boundary
        return Eth_out

    def compute_formation_energy(self, radius, density, pressure, solid, Lh=750e3, Lh_prime=0):
        '''
        Computes energy of formation for all solid points
        Returns the energy of formation
        args:
            radius: array of radii [m]
            density: array of densities [kg/m^3]
            pressure: array of pressures [Pa]
            solid: array of T/F values for solid or not
        '''
        dr = radius[1] - radius[0]
        Lh_at_P = pressure[solid] * Lh_prime + np.ones_like(pressure[solid]) * Lh
        return 4 * np.pi * dr * sum(Lh_at_P * density[solid] * radius[solid] ** 2.0)

    def compute_gravitational_energy(self, radius, density, pressure, temperature, solid):
        '''
        computes energy released to gravitational potential of solid inner core
        Returns the gravitational energy freed since last timestep
        args:
            radius: array of radii [m]
            density: array of densities [kg/m^3]
            temperature: array of temperatures [K]
            pressure: array of pressures [Pa]
            solid: arrayzx of T/F values for solid or not
        '''
        pass

    def adjust_boundaries(self,radius,density,target_masses,n_iter=5):
        layer_masses = self.compute_masses(radius,density)
        boundaries=[]
        for i in range(len(self.boundaries)):
            boundaries.append(((target_masses[i]/layer_masses[i]-1)/n_iter+1)*self.boundaries[i])
        self.set_boundaries(boundaries)
        return boundaries

    def check_solid(self, pressure, temperature, S_percent):
        '''
        Checks whether any point in P-T space is in the solidus.
        Returns a vector of true or false for solid phase at each point
        args:
            pressure: array of pressures [Pa]
            temperature: array of temperatures [K]
            radii: array of radii [m]
            S_percent: wt % of Sulfur in liquid
        '''
        solid = []
        for i, (P_in, T_in) in enumerate(zip(pressure, temperature)):
            solid.append(self.solidus_solver.check_solid(S_percent, P_in, T_in))
        return solid

    def integrate(self, n_slices, P0, T0, n_iter=5):
        '''
        Iteratively determine the pressure, temperature and gravity profiles for the
        planet.
        Usage:
            density, gravity, pressure, temperature, physical_params = integrate(n_slices,P0,T0,n_iter=5)
        args:
            n_slices: number of radial slices
            P0: initial central pressure in Pa
            T0: exterior temperature in K
            n_iter: number of iterations (default: 5)
            tol: not implemented
        '''

        radius = np.linspace(0.e3, self.boundaries[-1], n_slices)
        pressure = np.linspace(P0, 0.0, n_slices)  # initial guess at pressure profile
        temperature = np.ones_like(pressure) * T0
        gravity = np.empty_like(radius)
        for i in range(n_iter):
            density, physical_params = self.evaluate_eos(pressure, temperature, radius)
            gravity = self.compute_gravity(density, radius)
            pressure = self.compute_pressure(density, gravity, radius)
            temperature = self.compute_adiabat(pressure, radius)
        return radius, density, gravity, pressure, temperature, physical_params

    def integrate_heFESTo_mantle(self, n_slices, P0, T0, n_iter=5, calc_T_cmb_from_mantle=True):
        '''
        Iteratively determine the pressure, temperature and gravity profiles for the
        planet. Replaces mantle with heFESTo values using T(P) and density(P).
        Usage:
            density, gravity, pressure, temperature, physical_params = integrate(n_slices,P0,T0,n_iter=5)
        args:
            n_slices: number of radial slices
            P0: initial central pressure in Pa
            T0: exterior temperature in K
            n_iter: number of iterations (default: 5)
            calc_T_cmb_from_mantle: boolean, whether to use specified T_cmb or replace with
                T at bottom of mantle
            tol: not implemented
        '''

        radius = np.linspace(0.e3, self.boundaries[-1], n_slices)
        pressure = np.linspace(P0, 0.0, n_slices)  # initial guess at pressure profile
        temperature = np.ones_like(pressure) * T0
        gravity = np.empty_like(radius)
        boolean_mantle = radius > self.boundaries[-2]
        r_mantle = radius[boolean_mantle]
        r_cmb = self.boundaries[-2]
        r_surf = self.boundaries[-1]
        for i in range(n_iter):
            density, physical_params = self.evaluate_eos(pressure, temperature, radius)
            assert (pressure[boolean_mantle].max() < 24.0e9)
            index_vector = []
            for i, P in enumerate(pressure[boolean_mantle]):
                index_vector.append(self.inv_pressure_mantle_func(P))
            density[boolean_mantle] = self.density_mantle_func(index_vector)

            gravity = self.compute_gravity(density, radius)

            pressure = self.compute_pressure(density, gravity, radius)
            index_vector = []
            for i, P in enumerate(pressure[boolean_mantle]):
                index_vector.append(self.inv_pressure_mantle_func(P))
            if calc_T_cmb_from_mantle:
                temperature[boolean_mantle] = self.temperature_mantle_func(index_vector)
                self.T_bounds[-2] = temperature[boolean_mantle][0]
            temperature = self.compute_adiabat(pressure, radius)
            temperature[boolean_mantle] = self.temperature_mantle_func(index_vector)
            density[boolean_mantle] = self.density_mantle_func(index_vector)
        return radius, density, gravity, pressure, temperature, physical_params

    def integrate_conserve_mass(self, n_slices, P0, T0, target_masses, n_iter=5):
        '''
        Iteratively determine the pressure, temperature and gravity profiles for the
        planet by varying layer radius to conserve layer mass.
        Usage:
            density, gravity, pressure, temperature, physical_params = integrate(n_slices,P0,T0,n_iter=5)
        args:
            n_slices: number of radial slices
            P0: initial central pressure in Pa
            T0: exterior temperature in K
            n_iter: number of iterations (default: 5)
            tol: not implemented
        '''

        radius = np.linspace(0.e3, self.boundaries[-1], n_slices)
        pressure = np.linspace(P0, 0.0, n_slices)  # initial guess at pressure profile
        temperature = np.ones_like(pressure) * T0
        gravity = np.empty_like(radius)
        for i in range(n_iter):
            density, physical_params = self.evaluate_eos(pressure, temperature, radius)
            gravity = self.compute_gravity(density, radius)
            pressure = self.compute_pressure(density, gravity, radius)
            temperature = self.compute_adiabat(pressure, radius)
            self.adjust_boundaries(radius,density,target_masses, n_iter)
        return radius, density, gravity, pressure, temperature, physical_params

    def plot_all(self, radius, density, gravity, pressure, temperature, show_plot=True):
        '''
        Plots pressure, temperature, gravity, and density vs radius.
        '''
        f = plt.figure()
        ax1 = plt.subplot(221)
        ax2 = plt.subplot(222)
        ax3 = plt.subplot(223)
        ax4 = plt.subplot(224)
        plt.hold(True)

        # Plots !
        ax1.plot(radius / 1000., density)
        ax1.set_title('density', fontsize=10)
        ax1.set_xlabel('radius (km)', fontsize=10)
        ax1.set_ylabel('density (kg/m^3)', fontsize=10)

        ax2.plot(radius / 1000., gravity)
        ax2.set_title('gravity', fontsize=10)
        ax2.set_xlabel('radius (km)', fontsize=10)
        ax2.set_ylabel('gravity (m/s^2)', fontsize=10)

        ax3.plot(radius / 1000., pressure / 1e9)
        ax3.set_title('pressure', fontsize=10)
        ax3.set_xlabel('radius (km)', fontsize=10)
        ax3.set_ylabel('pressure (GPa)', fontsize=10)

        ax4.plot(radius / 1000., temperature)
        ax4.set_title('Temperature', fontsize=10)
        ax4.set_xlabel('radius (km)', fontsize=10)
        ax4.set_ylabel('temperature (K)', fontsize=10)
        if show_plot:
            plt.show()
        return f

    def plot_heFESTo(self, pressure, density, temperature, show_plot=True):
        '''
        Plots given heFESTo values over the model values to check implementation
        '''
        f = plt.figure()
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)

        plt.hold(True)
        xx = np.linspace(0., 1., 100)

        ax1.plot(pressure, density, label='Calculated')
        ax1.plot(self.pressure_mantle_func(xx), self.density_mantle_func(xx), label='heFESTo')
        ax1.legend()
        ax1.set_xlabel('Pressure (Pa)')
        ax1.set_ylabel('Density (kg/m^3)')
        ax1.set_title('Compare Model to heFESTo')

        ax2.plot(pressure, temperature, label='Calculated')
        ax2.plot(self.pressure_mantle_func(xx), self.temperature_mantle_func(xx), label='heFESTo')
        ax2.legend()
        ax2.set_xlabel('Pressure (Pa)')
        ax2.set_ylabel('Temperature (K)')
        ax2.set_title('Compare Model to heFESTo')

        if show_plot:
            plt.show()
        return f


