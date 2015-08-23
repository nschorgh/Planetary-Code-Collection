Craters/

Shadowing, illumination, and scattering on 3D topography
========================================================

shadows.f90:
	main program; calculates horizons

fieldofviews.f90:
	main program; calculates horizons and field of views

cratersQ_snapshot.f90:
	main program; instantaneous surface temperature with 3D shadowing and reflections

cratersQ_moon.f90:
	main program; surface temperature with 3D shadowing and reflections

topos.f90:
	input topography information

crater_modules.f90: 
	interface definitions

crater_common.f90: 
	common routines

shadow_subs.f90: 
	subroutines for shadows

fieldofview_subs.f90: 
	subroutines for fieldofviews

model_subs.f90: 
	subroutines for cratersQ_*

topo40.xyz: 
	example input topography

makefile:
	shows file dependencies




