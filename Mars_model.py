'''
Mars_model.py

A python script to calculate self-consistent pressure, temperature, density, and gravity
for Mars. Built using the burnman framework: http://www.burnman.org/. Utilizes
build_planet.py originally created by Ian Rose, modified by Nicholas Knezek. Also
utilizes solidus.py to calculate freezing points of Fe-FeS liqiuds in core.
Allows use of mantle parameters from heFESTo for integration instead of burnman minerals.

by Nicholas Knezek
CIDER, July 2014
'''


import numpy as np
import matplotlib.pyplot as plt

import FeS
import burnman
import burnman.minerals as minerals
import burnman.mineral_helpers as bmb
import burnman.composite as composite
from build_planet import Planet
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import RectBivariateSpline


# Material parameters - These are unrealistic
class iron_liquid(burnman.Mineral):
	'''
	Parameters for iron liquid. Most of these are guesses and approximations
	'''
	def __init__(self):
		self.params = {
			'equation_of_state':'mgd2',
			'V_0': 7.978e-6,  #Molar volume [m^3/(mole molecules)]
				#at reference pressure/temperature
			'K_0': 210.0e9, #Reference bulk modulus [Pa]
				#at reference pressure/temperature
			'Kprime_0': 3.5, #pressure derivative of bulk modulus
			'G_0': 1.e-8, #reference shear modulus
				#at reference pressure/temperature
			'Gprime_0': 0., #pressure derivative of shear modulus
			'molar_mass': 0.055845, #molar mass in units of [kg/mol]
			'n': 1, #number of atoms per formula unit
			'Debye_0': 10., #Debye temperature for material.
				# Set very low to approximate E~3R
			'grueneisen_0': 1.85,  #Gruneisen parameter for material.
			'q_0': 0.35, #isotropic strain derivative of gruneisen
			'eta_s_0': 0., #full strain derivative of gruneisen parameter
			'T_0': 1600., # Reference temperature
			'k': 86., # Thermal Conductivity [W/(K*m)]
			'sigma': 1.2e6,  # Electrical Conductivity, [S/m]
			'Cv': 815.0, # Heat Capacity at Constant Volume [J/(kg*K)]
			'Lh':750.e3, # Enthalpy of Formation at 0GPa [J/kg]
			'Lh_prime': 0., # Pressure Derivative of Enthalpy of Formation [J/kg/Pa]
			'S_frac': 0. # Weight fraction of Sulfur in material
			}

class ironsulfide_liquid(burnman.Mineral):
	'''
	Parameters for iron sulfide liquid. Most of these are guesses and approximations
	'''
	def __init__(self):
		self.params = {
			# Uses the same
			'equation_of_state':'mgd2',
			'V_0': 1.915e-05,  #Molar volume [m^3/(mole molecules)]
				#at reference pressure/temperature
			'K_0': 160.0e9, #Reference bulk modulus [Pa]
				#at reference pressure/temperature
			'Kprime_0': 3.5, #pressure derivative of bulk modulus
			'G_0': 1.e-8, #reference shear modulus
				#at reference pressure/temperature
			'Gprime_0': 0., #pressure derivative of shear modulus
			'molar_mass': 0.08791, #molar mass in units of [kg/mol]
			'n': 2, #number of atoms per formula unit
			'Debye_0': 10., #Debye temperature for material.
				# Set very low to approximate E~3R
			'grueneisen_0': 1.6,  #Gruneisen parameter for material.
			'q_0': 0.35, #isotropic strain derivative of gruneisen
			'eta_s_0': 0., #full strain derivative of gruneisen parameter
			'T_0': 1600., # Reference temperature
			'k': 86., # Thermal Conductivity [W/(K*m)]
			'sigma': 1.2e6,  # Electrical Conductivity, [S/m]
			'Cv': 815.0, # Heat Capacity at Constant Volume [J/(kg*K)]
			'Lh':750.e3, # Latent Heat of Formation at 0GPa [J/kg]
			'Lh_prime': 0, # Pressure Derivative of Enthalpy of Formation [J/kg/Pa]
			'S_frac': 0.3647 # Weight fraction of Sulfur in material
			}

class iron_ironsulfide_solidsolution(bmb.HelperSolidSolution):
	'''
	Fe-FeS liquid mixture created by simple solid solution averaging by burnman.
	'''
	def __init__(self, FeS_molfrac):
		base_materials = [iron_liquid(), ironsulfide_liquid()]
		molar_fraction = [1. - FeS_molfrac, 0.0 + FeS_molfrac] # keep the 0.0 +, otherwise it is an array sometimes
		bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class olivine(burnman.Mineral):
	'''
	Starting mineral for mantle. Isn't actually used if heFESTo values are inserted
	instead.
	'''
	def __init__(self):
		self.params = {
			'equation_of_state':'slb3',
			'V_0': 11.24e-6,
			'K_0': 161.0e9,
			'Kprime_0': 3.9,
			'G_0': 130.9e9,
			'Gprime_0': 1.92,
			'molar_mass': .0403,
			'n': 2,
			'Debye_0': 773.,
			'grueneisen_0': 1.5,
			'q_0': 1.5,
			'eta_s_0': 2.3
			}



'''
Define Parameters:
'''
# FeS Weight Percent Sulfur
wtS = 0.10

Mars_mass = 6.4185e23
Mars_g = 3.711

# Structural Parameters
r_cmb = 1650.e3 # core-mantle boundary radius in m (~1600e3 - 1850e3)
r_surf = 3390.e3 # surface radius in m

# Boundary temperatures for calculating adiabat.
# When using heFESTo parameters, T_surf is ignored
T_surf = 1100;
T_cmb = 2544;

# integration parameters
n_slices = 200
# Initial guess at central pressure and temperature (doesn't matter so much)
P0 = 40.0e9
T0 = 1500.0

# Create compositions for each layer
ol = olivine()
mol_FeS = FeS.wpS_to_mf(wtS)
fe_fes_ss = iron_ironsulfide_solidsolution(mol_FeS)

# Which heFESTo Mantle to use
# heFESTo_filename = None
heFESTo_filename = './input/DWHot.csv'


'''
Perform Calculations
'''

#### build planet!
# Mars = Planet([r_cmb,r_surf],[fe_fes_ss,ol],[T_cmb,T_surf],methods=['mgd2','slb3'])
Mars = Planet([r_cmb,r_surf],[fe_fes_ss,ol],[T_cmb,T_surf],methods=['mgd2','slb3'],heFESTo_filename=heFESTo_filename)

#### Integrate!
# radius, density, gravity, pressure, temperature, physical_params = Mars.integrate(n_slices,P0,T0,n_iter=3)
radius, density, gravity, pressure, temperature, physical_params = Mars.integrate_heFESTo_mantle(n_slices,P0,T0,n_iter=3,calc_T_cmb_from_mantle=True)

#### Check solidus in core
Rcore = radius[radius<r_cmb]
RHOcore = density[radius<r_cmb]
Pcore = pressure[radius<r_cmb]
Tcore = temperature[radius<r_cmb]
solid = Mars.check_solid(Pcore,Tcore,wtS)

#### Find energy of core
thermal_energy = Mars.compute_thermal_energy(Rcore,RHOcore,Tcore)
formation_energy = Mars.compute_formation_energy(Rcore,RHOcore,Pcore,solid)
# Gravitational buoyancy energy function not completed yet.
# gravitational_energy = Mars.compute_gravitational_energy(radius, density,pressure,temperature
Mass_core = Mars.calculate_planet_mass(Rcore,RHOcore)
print Mass_core
#### Calculate model gravity and mass error
Mars_g_model = gravity[-1]
g_err = float(Mars_g_model - Mars_g)/float(Mars_g)
Mars_mass_model = Mars.calculate_planet_mass(radius,density)
Mass_err = float(Mars_mass_model - Mars_mass)/float(Mars_mass)
print '\n\nsurface gravity error = '+str(g_err)
print 'Mass Error = '+str(Mass_err)


'''
Output and plot calculated values
'''

#### Plot
f = Mars.plot_all(radius,density,gravity,pressure,temperature,show_plot=True)
f2 = Mars.plot_heFESTo(pressure,density,temperature,show_plot=False)
f2.savefig('heFESTo.png')
f.savefig('MarsPTgrho.png')


#### Store Model Values in csv file
filename = "MarsValues_"+str(r_cmb/1.e3)+'km_'+str(wtS)+'wtS_'+str(T_cmb)+'Kcmb_DWhot'
fout = open(filename+".csv",'w')
fout.write("Mass error,"+str(Mass_err)+", Surface gravity error,"+str(g_err)+"\n")
fout.write("radius (m),density (kg/m^3),pressure (Pa),temperature (K),gravity (m/s^2), solid?, vp (m/s),vs (m/s),vphi (m/s),K (Pa),G (Pa)\n")
for r,rho,P,T,g,sol,vp,vs,vphi,K,G in zip(radius,density,pressure,temperature,gravity,solid,physical_params['vp'],physical_params['vs'],physical_params['vphi'],physical_params['K'],physical_params['G']):
	fout.write(str(r)+','+str(rho)+','+str(P)+','+str(T)+','+str(g)+','+str(sol)+','+str(vp)+','+str(vs)+','+str(vphi)+','+str(K)+','+str(G)+'\n')
fout.close





