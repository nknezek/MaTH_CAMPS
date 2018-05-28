# Core Thermal Section

## Notes

Past-tense, active, "we"



We created a way to track the thermal evolution of the core of Martian history. Previous authors have created thermal models of planetary cores to examine the evolution of Earth and Mars (e.g. Stevenson 1983, Nimmo 2015, Labrosse 2015). However, these authors all assume that the thermal profile within the core is adiabatic, which is no longer true if the heat flow out of the top of the core falls below the conductive gradient at the CMB for a significant period of time. A dense basal mantle layer as we proposed in our model insulates the core and suppresses heat flow for a significant period of Martian history. Therefore, the thermal profile departed from an adiabat significantly which required tracking the evolution of the full thermal profile within the Martian core. 

​	We computed the evolution of the thermal profile numerically by splitting the core into 150 spherically-symetric radial layers indexed by $i$. We then tracked the temperature $T$ evolution of each layer due to both conduction and convection. In each timestep, we conducted heat $Q^{cond}$ between neighboring layers based on their relative temperatures, using a finite volume formulation

$Q_{i}^{cond} = \frac{-A_{i} \,k\,(T_{i+1}-T_{i})}{\Delta r} $

where $A_i = 4 \pi r_i^2$ is the outer surface area of the layer with index $i$, and $r_i$ is the radial location of the outer boundary of the layer.  Then the change in temperature of each layer is computed using

$\Delta T_i = \frac{(Q_{i}^{cond}-Q_{i-1}^{cond}) \Delta t}{\rho_i C_p V_i}$

where $\Delta t$ is the timestep size, $\rho$ is density, $C_p$ is heat capacity, $V_i = 4/3 \pi (r_i^3-r_{i-1}^3)$ is the volume, and $T_i$ is the average temperature of layer $i$. At the outermost layer, $Q_{i}^{cond}$ is simply replaced by the heat transferred into the mantle, $Q_{CMB}$. On the other hand, at the centermost layer $A_{i-1}=0$ , so that $Q_{i-1}^{cond}=0$ and thus the centermost layer only loses heat by conduction to layers above. 

Then, we accounted for the comparatively rapid transfer of heat by convection. We computed the adiabatic gradient downward from each layer using an exponential gradient

$ T_{a,i}= T_je^{(r_{c,j}^2-r_{c,i}^2)/D^2}$

where $T_j$ and $r_{c,j}$ are the temperature and center radial location of layer $j$, $D$ is the adiabatic length scale in the core, and $i$ is the new layer to calculate the relative adiabaitic temperature. 

If the next lower layer had a temperature greater than or equal than the local adiabatic temperature extrapolated to that depth

$T_i >=T_{a,i}|_{T_{i+1}}$

they two layers were defined to be in a convecting region. This metric was computed for all layers, and multiple consecutive layers satisfying this criteria were grouped together into one convective region. This method allowed for mutliple convecting and conducting regions within the core, which could arise due to highly variable CMB heat flow. 

In a convective region, heat was transported instantaneously upward each timestep until the temperature of all layers in the convecting region matched the local adiabat. The amount of heat transferred from the lower layers to the top layer is chosen to maintain the total energy of the convecting region, computed by

$ E = \int_r \rho(r) c_p T(r) dV = \sum_i \rho_i c_p T_i V_i$

Over multiple timesteps, convective regions will cool off due to conductive heat transfer out of their top layer. If the CMB heat flow is large enough so that the entire core convects, this method gives the same result as the approach used by previous authors. The core maintains an adiabatic gradient while losing a total internal energy equivalent to the heat lost to the mantle across the CMB, which slowly decreases the temperature over time. 

However, this method also allowed the temperature profile within the core to depart significantly from an adiabat while conserving energy properly. We found that in many runs a stratified layer would grow to encompass the majority of the core, and in some runs, the core became nearly isothermal due to lower-mantle layer insulation and radiogenic heating. 

