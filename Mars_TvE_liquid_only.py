'''
Mars_TvE_liquid_only.py
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


# Material parameters - These are mainly guesses
class iron_liquid(burnman.Mineral):
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
			'DEth':270., # Enthalpy of Formation at 0GPa [J/kg]
			'DEth_prime': 1e-6, # Pressure Derivative of Enthalpy of Formation [J/kg/Pa]
			'S_frac': 0. # Weight fraction of Sulfur in material
			}

class ironsulfide_liquid(burnman.Mineral):
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
			'DEth':270., # Enthalpy of Formation at 0GPa [J/kg]
			'DEth_prime': 1e-6, # Pressure Derivative of Enthalpy of Formation [J/kg/Pa]
			'S_frac': 0.3647 # Weight fraction of Sulfur in material
			}

class iron_ironsulfide_solidsolution(bmb.HelperSolidSolution):
	def __init__(self, FeS_molfrac):
		base_materials = [iron_liquid(), ironsulfide_liquid()]
		molar_fraction = [1. - FeS_molfrac, 0.0 + FeS_molfrac] # keep the 0.0 +, otherwise it is an array sometimes
		bmb.HelperSolidSolution.__init__(self, base_materials, molar_fraction)

class olivine(burnman.Mineral):
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


fe_l = iron_liquid()
fes_l = ironsulfide_liquid()

# FeS Molar fraction
wtS = 0.14
# mol_FeS = FeS.wpS_to_mf(0.04) 0.073 # 4.% wt fraction S (small core)
# mol_FeS = FeS.wpS_to_mf(0.2)  # 20.% wt fraction S (large core)
mol_FeS = FeS.wpS_to_mf(wtS)


fe_fes_ss = iron_ironsulfide_solidsolution(mol_FeS)
ol = olivine()
mg_pv = minerals.SLB_2011.mg_perovskite()

Mars_mass = 6.4185e23
Mars_g = 3.711

# Structural Parameters
r_cmb = 1500e3
r_surf = 3390e3

# These don't matter, T_surf set by heFESTo and T_cmb is iterated
T_surf = 1100;
T_cmb = 2050;

# integration parameters
n_slices = 100
P0 = 40.0e9
T0 = 0.




# build planet!
# Mars = Planet([r_cmb,r_surf],[fe_fes_ss,ol],[T_cmb,T_surf],methods=['mgd2','slb3'])
Mars = Planet([r_cmb,r_surf],[fe_fes_ss,ol],[T_cmb,T_surf],methods=['mgd2','slb3'],heFESTo_filename='DWHot.csv')

# # Integrate!
# radius, density, gravity, pressure, temperature, physical_params = Mars.integrate(n_slices,P0,T0,n_iter=3)
radius, density, gravity, pressure, temperature, physical_params = Mars.integrate_heFESTo_mantle(n_slices,P0,T0,n_iter=3)
thermal_energy = Mars.compute_thermal_energy(density,temperature,radius)

# Find Error in gravity and Mass
Mars_g_model = gravity[-1]
g_err = float(Mars_g_model - Mars_g)/float(Mars_g)
Mars_mass_model = Mars.calculate_planet_mass(radius,density)
Mass_err = float(Mars_mass_model - Mars_mass)/float(Mars_mass)
print '\n\nsurface gravity error = '+str(g_err)
print 'Mass Error = '+str(Mass_err)



# Plot
f = Mars.plot_all(radius,density,gravity,pressure,temperature,show_plot=False)
f2 = Mars.plot_heFESTo(pressure,density,temperature,show_plot=False)
f2.savefig('heFESTo.png')
f.savefig('MarsPTgrho.png')
#
# Store Model Values in csv file
filename = "MarsValues_"+str(r_cmb/1.e3)+'km_'+str(wtS)+'wtS_DWhot'
fout = open(filename+".csv",'w')
fout.write("Total Mass = "+str(Mars_mass)+", Surface Gravity = "+str(Mars_g)+"\n")
fout.write("radius (m),density (kg/m^3),pressure (Pa),temperature (K),gravity (m/s^2),vp (m/s),vs (m/s),vphi (m/s),K (Pa),G (Pa)\n")
for r,rho,P,T,g,vp,vs,vphi,K,G in zip(radius,density,pressure,temperature,gravity,physical_params['vp'],physical_params['vs'],physical_params['vphi'],physical_params['K'],physical_params['G']):
	fout.write(str(r)+','+str(rho)+','+str(P)+','+str(T)+','+str(g)+','+str(vp)+','+str(vs)+','+str(vphi)+','+str(K)+','+str(G)+'\n')
fout.close



####### Calculate T vs E values #######
Etherm = []
g_err_vec = []
Mass_err_vec = []
Tcenter = []
Pcenter = []
P_cmb = []
#
Npoints = 101
Temps = np.linspace(3000.0,1500.0,Npoints)
#
for T in Temps:
	Mars.set_boundary_temperatures([T,T_surf])
	radius, density, gravity, pressure, temperature, physical_params = Mars.integrate_heFESTo_mantle(n_slices,pressure,temperature,n_iter=3,calc_T_cmb_from_mantle=False)
	Etherm.append(Mars.compute_thermal_energy(density,temperature,radius)[0])
	Tcenter.append(temperature[0])
	Pcenter.append(pressure[0])
	Mars_g_model = gravity[-1]
	g_err_vec.append(float(Mars_g_model - Mars_g)/float(Mars_g))
	P_core = pressure[radius < r_cmb]
	P_cmb.append(P_core[-1])

	Mars_mass_model = Mars.calculate_planet_mass(radius,density)
	Mass_err_vec.append(float(Mars_mass_model - Mars_mass)/float(Mars_mass))

	print '\nCompleted Temperature: ' + str(T)

# Store Eth in csv file
filename = "TvE"
fout = open(filename+".csv",'w')
fout.write("Temperature at CMB adiabat vs core internal energy. 14wt% Sulfur \n")
fout.write("Core Internal Energy, Temperature at CMB, Pressure at CMB, Temperature at Center, Pressure at Center,Mass Error, Surface Gravity Error\n")
for E,T,Pcm,Tc,Pc,dm,dg in zip(Etherm,Temps,P_cmb,Tcenter,Pcenter,Mass_err_vec,g_err_vec):
	fout.write(str(E)+','+str(T)+','+str(Pcm)+','+str(Tc)+','+str(Pc)+','+str(dm)+','+str(dg)+'\n')
fout.close



####### Calculate Inner Core Freezing #########
# S_percent = 0.072
# solid = Mars.check_solid(pressure,temperature,radius,S_percent)
# print solid





