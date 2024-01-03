Common/Misc/

Subroutines for experimental purposes
=====================================

conductionQ+iter.f90:
	1D thermal conduction with heterogeneous thermal properties and 
	Stefan-Boltzmann Law boundary condition, semi-implicit solver, with
	built-in predictor-corrector that iteratively adjusts linearization 
	temperature. (This was the default version of conductionQ.f90 for a 
	number of years.)

conductionQ+smooth.f90:
	1D thermal conduction with heterogeneous thermal properties and 
	Stefan-Boltzmann Law boundary condition, semi-implicit solver, plus
	wrapper subroutine that peforms artificial flux smoothing for extra
	stability
	(also provides subroutine cranknQ)

conductionQ+volt.f90:
	1D thermal conduction with heterogeneous thermal properties and 
	Stefan-Boltzmann Law boundary condition, semi-implicit solver, with 
	Volterra predictor for linearization temperature
	
makefile:
	shows file dependencies

