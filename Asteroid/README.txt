Asteroid/

Basic Thermal Model for Astroid Surfaces
========================================

asteroid_thermal.f90;
	main program	

oneasteroid_thermal1d.f90:
	simple 1D diffusion of temperture for asteroid

insolonly.f90:
	Insolation only, very simple


Asynchronously Coupled Model of Ice Evolution, Temperature, and Impact Stirring for Astroid Surfaces
====================================================================================================

asteroid_fast2.f90: 
       	main program

impactstirring.f90:
	Model of ice homegenization due to impact stirring,
	for use with "fast" model

fast_modules_asteroid.f90 
	modules for fast method

fast_subs_asteroid.f90 
	subroutines and functions for fast method



