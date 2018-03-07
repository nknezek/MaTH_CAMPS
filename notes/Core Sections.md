MaTH-CAMPS paper outline

## Core Section

Key points

- 11-16wt% S proposed
  - depresses liquidus and keeps core molten. 
  - allows for potentially snowing cores in some S compositions
  - allows for completely molten cores in our model
- Layered model predicts strongly sub-adiabatic heatflow through most of martian history
  - leads to strongly stratified core, almost isothermal for some periods of Martian history
  - requires development of code to treat both convection and conduction through core, dependent on heat flow across CMB.

#### Core Thermal Model

- 1D spherical parameterized convection-conduction code. 
- Parameterized convection, finite-volume conduction
- Euler stepping with timestep dt
  - dt ~ 1 Myr found to be numerically stable for reasonable heat flows
  - dt ~10 kyr numerically stable even for implausibly large heatflows out of or into the core, and solutions within 1% of results from 100kyr runs, so errors from euler method assumed to be insignificant compared to other systemic erros in model.
- **Conduction**
  - conducts heat through layers, maintaining energy balance exactly. 
  - 150 layers found to achieve good accuracy.
  - finite-volume method in radius, spherically symmetric
  - parameters
    - $k $  thermal conductivity =constant130 W/m-k from Nimmo 2015
    - $\rho$ density = parameterized after the method of Nimmo 2015
      - varies from central rho to smaller rho at CMB. Low density chosen for low P and high S content
      - TODO: maybe calculate densities of solid-solution hard-sphere Fe-S mixtures at various P-T and correct for temperatures?
    - $P $ Pressure - varies following Nimmo 2015
      - $P_{CMB}$ given from mantle code based on core size
    - $\alpha$ thermal expansivity - 1.25e-5 /K from Nimmo 2015
    - $T_{adiabat}$ - follows Nimmo 2015 for radial variation 
    - $\alpha_c$ Compositional expansivity 1.1 from Nimmo 2015 (TODO: maybe need better value for FeS?)
    - $L_H$ Latent heat of solidification 750 kJ/kg Nimmo 2015 (for Earth, but cite Stevenson 1983 too?)
- **Convection**
  - at a particular timestep, a region of the core is defined as "convecting" if the radial temperature is greater than or equal to the local adiabatic gradient.
    ​	$ \nabla_r T(r) \ge \nabla_r T_a(r)$
  - For convecting regions, the change in temperature of the entire region is equal to the heatflow out of the top layer into either the overlying mantle or an overlying conducting stratified layer. 
    ​	$q_{region} = \nabla_rT(r_{top})$
    - At the CMB, the heatflow out of the top of the core is set by the mantle. For a convective region beneath a stratified layer, the heatflow is set by the conduction across the convective-stratified boundary. 
  - This heatflow is equal to the change in total internal energy of the convective region, as the shape of the adiabat dictates that the greatest conductive heatflow occurs at the largest radius, and the rest of the region convects energy up the adiabat instantaneously on geologic timescales. 
  - Internal energy is calculated
    - $ E = \int_r \rho(r) c_p T(r) dV = \sum_i \rho_i c_p T(r_i) dV_i$
  - The new temperature gradient after a timestep is set by choosing a temperature of the upper layer and tracing the local adiabat throughout the convective region. This T is chosen such that the heat flow equals the change in internal energy.
    - $ E = \int_r \rho(r) c_p T(r) dV = \sum_i \rho_i c_p T(r_i) dV_i$

#### Fe-S Phase Diagram

- Limited data availability in region of interest, need good phase diagram.
  - Potential for solidificaiton at low wt% S (5-8% ) with our model. 
    - could allow for re-start of Martian dynamo, affects thermal evolution of coupled core-mantle.
  - potentially snowing cores /mostly solid cores in canonical Mars model that does not include layer
    - introduces all sorts of interesting effects (TODO:CITE)
- New data available (TODO:CITE)
  - fit curves to data points from all previous authors
- Phase diagram construction 
  - parameterize phase diagram with 2D splines in pressure and sulfur space. 
  - functional form of eutectic between $Fe-L$ and $Fe_{(3+x)}S_2-L$ regions

