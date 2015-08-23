Planetary Code Collection

by Norbert Schorghofer


This program collection contains
* a semi-implicit one-dimensional thermal model for planetary surfaces,
* an explicit subsurface vapor diffusion and deposition model,
* models for the equilibrium ice table on Mars,
* fast method for subsurface-atmosphere vapor exchange for Mars,
* a Monte-Carlo model for ballistic hops in the lunar exosphere, and
* a 3D model of shadowing, illumination, and scattering.


==Basic Models for Planetary Surfaces==

*General Purpose*

flux_noatm.f90: 
	Insolation on horizontal or sloped surface without atmosphere

psv.f: 
	vapor pressure of H2O

psvco2.f: 
       	vapor pressure of CO2

generalorbit.f: 
	distance, longitude, and declination of the sun from orbital elements

conductionQ.f:
	1D thermal conduction with heterogeneous thermal properties and 
	flux/radiation boundary condition, semi-implicit solver

julday.for: 
	from Numerical Recipes (C), but without the pause statement

tridag.f: 
	from Numerical Recipes (C), without stop and pause statements, NMAX=1000

derivs.f90: 
	first and second derivatives on irregular grid

grids.f: 
	creates appropriate 1D grids, calculates thermal properties of soil with ice



*Basic Thermal Model*

asteroid_thermal.f90: 
	main program

oneasteroid_thermal1d.f90: 
	1D diffusion of temperature for asteroid

insolonly.f90:
	insolation only



==Mars Subsurface Ice Model (M-SIM)==

*Mars Thermal Model*

flux.f: 
	Insolation on horizontal or sloped surface with a poor man's Mars atmosphere

marsorbit.f:
	position of the sun as seen from Mars; data from Allison & McEwen (2000)

conductionT.f:
	1D thermal conduction with heterogeneous thermal properties and 
	temperature boundary condition, semi-implicit solver

tprofile.m: 
	Matlab script that compares temperature profile with  analytic solution

modeldescription.pdf:
	Notes on Numerics (see Part 1)

mars_thermal1d.f: 
	1D diffusion of temperature for Mars; 
	prototype example of how to call conductionQ/T with seasonal CO2 frost cover

input.par: 
	Example input file for mars_thermal1d.f



*Vapor Diffusion Model*

vapordiffusioni.f:
	Diffusion of water vapor with phase transitions on irregular grid, explicit solver

adsorption.f:
	amount of adsorbed H2O and more

modeldescription.pdf:
	Notes on Numerics (see Part 2)

exper_thermal1d.f:
	1D diffusion of temperature and optionally also water vapor with prescribed 
	surface temperature

exper.par:
	Example input file for exper_thermal1d.f



*Equilibrium Ice Table (on Mars)*

jsub.f:
	net flux between ice table and surface, includes thermal model

jsubv.f90:
	vectorized version of jsub, includes emission from one surface to another

mars_mapi.f:
	determines equilibrium ice table depth for a list of locations (such as the entire 
	globe); contains leftover Numerical Recipes code (C)

mapgrid.dat:
	Example input file for mars_mapi.f

mars_mapt2.f:
	calls jsub for a list of locations (such as the entire globe)

mars_mapiqh2v.f90:
	version of mars_mapi that uses jsubv and slope coupling, configured for cluster



*Fast Method for Subsurface Ice Dynamics (on Mars)*

fast_modules.f90:
	numerically accelerated routines for growth and depletion of subsurface ice, 
	Fortran modules

fast_subs_univ.f90:
	numerically accelerated routines for growth and depletion of subsurface ice, 
	general subroutines

fast_subs_exper.f90:
	numerically accelerated routines for growth and depletion of subsurface ice

fast_subs_mars.f90:
	numerically accelerated routines for growth and depletion of subsurface ice

stabgrow_fast.f90:
	numerically accelerated growth of pore ice

exper_fast.f90:
	numerically accelerated growth and depletion of subsurface ice

input_fast.par:
	Example input file for stabgrow_fast.f90 and exper_fast.f90

mars_fast.f90:
	numerically accelerated growth and depletion of subsurface ice

lats.ph:
	Example input file for mars_fast.f90



==Lunar Models==

*Monte-Carlo Model for Surface-bounded Exosphere*

Exosphere/moon4.f90:
	main program; event driven Monte Carlo model for ballistic hops of 
	water molecules on the lunar surface

Exosphere/montecarlo.f90:
	ballistic hops, event scheduler

Exosphere/geogrid.f90:
	everything specific to the geographic grid

Exosphere/geogrid_D.f90:
	a different geographic grid

Common/subl_subs.f90:
	miscellaneous physical parametrizations

Common/gasdev.for:
	Gaussian probability distribution, Numerical Recipes(C)

Common/ran2.for:
	random number generator, Numerical Recipes(C)



*Shadowing and illumination*

see Craters/README.txt




ACKNOWLEDGMENTS:

2006: Troy Hudson discovered a grid-point offset in conductionT and conductionQ, which has been corrected.

2005: Thanks to Mischa Kreslavsky for providing correct formulas for energy balance on a slope.

* Many Thanks to Andy Vaught for developing an open-source Fortran 95 compiler (www.g95.org).  The early versions of this code were developed with this compiler.

2001: Samar Khatiwala invented an elegant implementation of the upper radiation boundary condition for the Crank-Nicolson method.

SUPPORT: This code development was supported by NASA, Caltech, and the University of Hawaii. Undoubtedly, some parts were written without support.



