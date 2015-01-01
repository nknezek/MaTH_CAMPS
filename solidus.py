import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

import burnman
import burnman.minerals as minerals
import burnman.mineral_helpers as bmb
import burnman.composite as composite
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import UnivariateSpline


class Solver:
	def __init__(self):
		'''
		Class to calculate Solidus for Fe/FeS mixtures
		'''
		# Low Sulfur Freeze P-T-S Values
		S_l = [0.,3.,6.,9.,12.]
		P_l = [0.0,10.,14.,23.,40.]

		S_14GPa_l = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,18.27,]
		T_14GPa_l = [2219,2071,1954,1863,1793,1743,1708,1685,1669,1658,1648,1635,1615,1585,1542,1482,1400,1294,1160,1119,]
		S_21GPa_l = [0.,3.,6.,9.,12.,15.4]
		T_21GPa_l = [2271.,2185.,2081.,1949.,1758.,1348.]
		S_23GPa_l = [0.,3.,6.,9.,12.,16.]
		T_23GPa_l = [2308.,2211.,2106.,1985.,1819.,1435.]
		S_40GPa_l = [0.,3.,6.,9.,12.]
		T_40GPa_l = [2498.,2372.,2215.,1997.,1565.]

		# Find T freeze for 30 GPa
		x = [23.,40.]
		S_30GPa_l = [0.,3.,6.,9.,12.]
		T_30GPa_l = np.empty_like(S_40GPa_l)
		T23_for_loop = [2308.,2211.,2106.,1985.,1819.]
		for i,(T40,T23) in enumerate(zip(T_40GPa_l,T23_for_loop)):
			T_30GPa_l[i] = np.interp(30.,x,[T23, T40])


		# High Sulfur P-T-S Values
		S_14GPa_mid = [19,20,20.8]
		T_14GPa_mid = [1144,1160,1173]
		S_14GPa_h = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36]
		T_14GPa_h = [1236,1394,1471,1528,1576,1618,1656,1690,1722,1752,1781,1808,1834,1858,1882,1905]
		S_21GPa_h = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,]
		T_21GPa_h = [1370,1405,1428,1443,1454,1461,1467,1472,1477,1480,1532,1586,1637,1685,1731,1775,1818,1859,1898,1937,1974,]
		S_23GPa_h = [16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,]
		T_23GPa_h = [1435,1450,1461,1470,1485,1492,1555,1614,1669,1720,1768,1812,1852,1888,1921,1950,1976,1997,2015,2029,2040,]
		S_30GPa_h = [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,]
		T_30GPa_h = [1500,1575,1646,1707,1759,1804,1842,1876,1906,1934,1958,1982,2003,2023,2042,2060,2077,2093,2109,2124,2138,2151,2164,]
		S_40GPa_h = [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,]
		T_40GPa_h = [1647,1730,1787,1833,1873,1908,1941,1971,1999,2026,2051,2075,2098,2120,2142,2162,2182,2202,2220,2239,2256,2274,2291,2308,]

		S_21GPa_h1 = S_21GPa_h[0:10]
		S_21GPa_h2 = S_21GPa_h[10:-1]
		T_21GPa_h1 = T_21GPa_h[0:10]
		T_21GPa_h2 = T_21GPa_h[10:-1]
		S_23GPa_h1 = S_23GPa_h[0:5]
		S_23GPa_h2 = S_23GPa_h[5:-1]
		T_23GPa_h1 = T_23GPa_h[0:5]
		T_23GPa_h2 = T_23GPa_h[5:-1]

		self.p14lsp = UnivariateSpline(S_14GPa_l,T_14GPa_l,k=2)
		self.p21lsp = UnivariateSpline(S_21GPa_l,T_21GPa_l,k=2)
		self.p23lsp = UnivariateSpline(S_23GPa_l,T_23GPa_l,k=2)
		self.p30lsp = UnivariateSpline(S_30GPa_l,T_30GPa_l,k=2)
		self.p40lsp = UnivariateSpline(S_40GPa_l,T_40GPa_l,k=2)
		self.p14msp = UnivariateSpline(S_14GPa_mid,T_14GPa_mid,k=2)
		self.p14hsp = UnivariateSpline(S_14GPa_h,T_14GPa_h,k=2)
		self.p21hsp1 = UnivariateSpline(S_21GPa_h1,T_21GPa_h1,k=2)
		self.p21hsp2 = UnivariateSpline(S_21GPa_h2,T_21GPa_h2,k=2)
		self.p23hsp1 = UnivariateSpline(S_23GPa_h1,T_23GPa_h1,k=2)
		self.p23hsp2 = UnivariateSpline(S_23GPa_h2,T_23GPa_h2,k=2)
		self.p30hsp = UnivariateSpline(S_30GPa_h,T_30GPa_h,k=2)
		self.p40hsp = UnivariateSpline(S_40GPa_h,T_40GPa_h,k=2)

	'''
	Define Functions to find temperature given sulfur at each pressure 14, 21, 23, 30, 40
	'''
	## Low S side of eutectic
	def T_p14l(self,S):
		return self.p14lsp(S)

	def T_p21l(self,S):
		return self.p21lsp(S)

	def T_p23l(self,S):
		return self.p23lsp(S)

	def T_p30l(self,S):
		return self.p30lsp(S)

	def T_p40l(self,S):
		return self.p40lsp(S)

	## High S side of eutectic
	def T_p14h(self,S):
		if S<21.:
			return self.p14msp(S)
		else:
			return self.p14hsp(S)

	def T_p21h(self,S):
		if S<25.:
			return self.p21hsp1(S)
		else:
			return self.p21hsp2(S)

	def T_p23h(self,S):
		if S<21.:
			return self.p23hsp1(S)
		else:
			return self.p23hsp2(S)

	def T_p30h(self,S):
		return self.p30hsp(S)

	def T_p40h(self,S):
		return self.p40hsp(S)

	def T_SPl(self,S,P):
		Pfit = [14.,21.,23.,30.,40.]
		Tfit = [self.T_p14l(S),self.T_p21l(S),self.T_p23l(S),self.T_p30l(S),self.T_p40l(S)]
		return float(UnivariateSpline(Pfit,Tfit,k=1,s=0)(P))

	def T_SPh(self,S,P):
		Pfit = [14.,21.,23.,30.,40.]
		Tfit = [self.T_p14h(S),self.T_p21h(S),self.T_p23h(S),self.T_p30h(S),self.T_p40h(S)]
		return float(UnivariateSpline(Pfit,Tfit,k=1,s=0)(P))

	def T_SP(self,S,P):
		'''
		Finds solidus temperature, given a sulfur wt% and pressure
		args:
			S: sulfur weight percent (0.0 to 1.0)
			P: pressure (Pa)
		'''
		P = P/(1.e9)
		S = S*100.
		assert(S>=0. and S<=35.)
		assert(P>=0. and P<=45.)
		if S<12.:
			return self.T_SPl(S,P)
		elif S>18.27:
			return self.T_SPh(S,P)
		else:
			return max(self.T_SPl(S,P),self.T_SPh(S,P))

	def check_solid(self,S,P,T):
		'''
		checks whether point in S,P,T space is above or below solidus
		args:
			S: sulfur weight percent (0.0 to 1.0)
			P: pressure (Pa)
			T: temperature (K)
		'''
		return bool(self.T_SP(S,P) > T)

