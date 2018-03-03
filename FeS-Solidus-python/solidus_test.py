import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

import solidus as s


Sinterp = np.linspace(0.,35.0,100)/100.
Pinterp = np.array([14.,20.,21.,23.,27.,30.,33.,37.,40.])*1.e9
Pall = np.linspace(14.,40.,27)*1.e9

Porig = np.array([14.,21.,23.,30.,40.])*1.e9

sol = s.Solver()
Tplot = dict()

for P in Pall:
	Ttemp = []
	for S in Sinterp:
		Ttemp.append(sol.T_SP(S,P))
	Tplot[str(P)] = Ttemp
	plt.plot(Sinterp,Ttemp,'k--')

for P in Porig:
	Ttemp = []
	for S in Sinterp:
		Ttemp.append(sol.T_SP(S,P))
	Tplot[str(P)] = Ttemp
	plt.plot(Sinterp,Ttemp,'-',label=str(P/1.e9)+' GPa')

plt.xlabel('Sulfur (wt%)')
plt.ylabel('Temperature (K)')
plt.title('Fe-FeS Liquidus')
plt.legend(loc=3,fontsize=10)
plt.grid(True)
# plt.show()
plt.savefig('SolidusTest14vals.png')