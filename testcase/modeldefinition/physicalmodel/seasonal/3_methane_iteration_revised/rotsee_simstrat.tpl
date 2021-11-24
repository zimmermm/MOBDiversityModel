ptf #
*** Files *************************************************
Rotsee_Input/rotsee_initialconditions.dat
Rotsee_Input/rotsee_grid.dat
Rotsee_Input/rotsee_morphology.dat
Rotsee_Input/rotsee_forcing.dat
Rotsee_Input/rotsee_absorption.dat
Rotsee_Output/
Rotsee_Input/rotsee_outputdepths.dat
Rotsee_Input/rotsee_outputtime.dat
Rotsee_Input/rotsee_Qinp.dat
Rotsee_Input/rotsee_Qout.dat
Rotsee_Input/rotsee_Tinp.dat
Rotsee_Input/rotsee_Sinp.dat
*** Model set up: time steps and grid resolution ***
300.0   	Timestep dt [s]
16147.0 	Start time [d]
17198.0 	End time [d]
*** Model, conditions and output selection ****************
1.0     	Turbulence model (1:k-epsilon, 2:MY)
1.0     	Stability function (1:constant, 2:quasi-equilibrium)
1.0     	Flux condition (0:Dirichlet condition, 1:no-flux)
3.0     	Forcing (1:Wind+Temp+SolRad, 2:(1)+Vap, 3:(2)+Cloud, 4:Wind+HeatFlux+SolRad)
0.0     	Use filtered wind to compute seiche energy (0/default:off, 1:on) (if 1:on, one more column is needed in forcing file)
2.0     	Seiche normalization (1:max N^2, 2:integral)
3.0     	Wind drag model (1/default:lazy (constant), 2:ocean (increasing), 3:lake (Wüest and Lorke 2003))
0.0     	Inflow placement (0/default:manual, 1:density-driven)
0.0     	Pressure gradients (0:off, 1:Svensson 1978, 2:?)
1.0     	Enable salinity transport (0:off, 1/default:on)
0.0     	Display simulation (0:off, 1:when data is saved, 2:at each iteration, 3:extra display)
1.0     	Display diagnose (0:off, 1:standard display, 2:extra display)
10.0    	Averaging data (not implemented, set to 10)
*** Model parameters **************************************
#lat       #	Lat [°]         Latitude for Coriolis parameter
#p_air     #	p_air [mbar]	Air pressure
#a_seiche  #	a_seiche [-]	Fraction of wind energy to seiche energy
#q_NN      #	q_NN			Fit parameter for distribution of seiche energy
#f_wind    #	f_wind	[-]		Fraction of forcing wind to wind at 10m (W10/Wf)
#C10       #	C10 [-]			Wind drag coefficient (used if wind drag model is 1:lazy)
#CD        #	CD [-]			Bottom friction coefficient
#fgeo      #	fgeo [W/m2]		Geothermal heat flux
#k_min     #	k_min [J/kg]	Minimal value for TKE
#p1        #	p1				Fit parameter for absorption of IR radiation from sky
#p2        #	p2				Fit parameter for convective and latent heat fluxes
#beta      #	beta [-]		Fraction of short-wave radiation directly absorbed as heat
#albsw     #	albsw [-]		Albedo for reflection of short-wave radiation