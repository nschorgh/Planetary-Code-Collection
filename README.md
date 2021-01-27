Planetary Code Collection
=========================

*by Norbert Schorghofer*


This program collection contains:

* semi-implicit one-dimensional thermal model for planetary surfaces
* explicit subsurface vapor diffusion and ice deposition model
* long-term ice evolution:  
  * equilibrium ice table on Mars
  * fast (asynchronously coupled) method for subsurface-atmosphere vapor exchange on Mars
  * fast (asynchronously coupled) method for ice retreat on asteroids
  * lunar ice pump (boundary-value formulation)   
* 3D model of shadowing, illumination, and scattering for airless bodies, Mars, and Mauna Kea
* Monte-Carlo model for ballistic hops in the exospheres of the Moon and Ceres



Mars Subsurface Ice Model (M-SIM)
---------------------------------

### Mars Thermal Model

Mars/mars_thermal1d.f: (main program)  
Mars/flux_mars.f90  
Mars/marsorbit.f90  
Mars/soilthprop_mars.f90  
Common/conductionQ.f90  
Common/conductionT.f90  
Common/grids.f90  
Common/julday.for  
Common/psvco2.f  
Common/tridag.for  
*Documentation: User Guide Part 1*  


### Vapor Diffusion Model

Mars/exper_thermal1d.f: (main program)  
Mars/vapordiffusioni.f  
Mars/adsorption.f  
*Documentation: User Guide Part 2  
Documentation: Schorghofer, N. & Aharonson, O. (2005) J. Geophys. Res. 110, E05003, Appendices*  
Mars/Misc/  an animation for illustration  


### Equilibrium Ice Table

Mars/mars_mapi.f: (main program)  
Mars/mars_mapii.f90: (main program)  
Mars/mars_mapt.f: (main program)  
Mars/mars_mapip2.f90: (main program)  
Mars/jsub.f90    
Mars/jsubv.f90  
Mars/marsorbit.f90  
Mars/flux_mars.f90  
Mars/soilthprop_mars.f90  
Common/conductionQ.f90  
Common/conductionT.f90  
Common/grids.f90  
Common/julday.for  
Common/psv.f  
Common/psvco2.f  
Common/tridag.for  
*Documentation: User Guide Section 3.1  
Documentation: Schorghofer, N. & Aharonson, O. (2005) J. Geophys. Res. 110, E05003*  


### Mars Long-Term Thermal Model

Mars/insol_driver.f90 (main program)  
Mars/tempr_driver.f90 (main program)  
Mars/flux_mars.f90  
Mars/generalorbit.f90  
Mars/soilthprop_mars.f90  
Common/conductionQ.f90  
Common/conductionT.f90  
Common/grids.f90  
Common/psvco2.f  
Common/tridag.for  
Mars/MilankOutput/  surface temperatures from last 21Myr  


### Fast (asynchronously-coupled) Method for Subsurface Ice Dynamics

Mars/stabgrow_fast.f90: (main program)  
Mars/exper_fast.f90: (main program)  
Mars/mars_fast.f90: (main program)  
Mars/fast_modules.f90  
Mars/fast_subs_univ.f90  
Mars/fast_subs_exper.f90  
Mars/fast_subs_mars.f90  
Mars/soilthprop_mars.f90  
Common/conductionQ.f90  
Common/conductionT.f90  
Common/derivs.f90  
Common/generalorbit.f  
Common/grids.f90  
Common/psv.f  
Common/psvco2.f  
Common/tridag.for  
*Documentation: Schorghofer, N. (2010) Icarus 208, 598-607*  


Other Models for Planetary Surfaces
-----------------------------------

### Basic Thermal Model for Airless Bodies

Asteroids/asteroid_thermal.f90: (main program)  
Asteroids/oneasteroid_thermal1d.f90  
Asteroids/insolonly.f90   
Common/conductionQ.f90  
Common/flux_noatm.f90  
Common/generalorbit.f
Common/grids.f90
Common/tridag.for  
*Documentation: User Guide Part 1*  


### Asynchronous Models for Temperature, Impact Stirring, and Ice Loss on Asteroids

Asteroids/asteroid_fast1.f90: (main program)  
Asteroids/asteroid_fast2.f90: (main program)  
Asteroids/fast_modules_asteroid.f90  
Asteroids/fast_subs_asteroid1.f90  
Asteroids/fast_subs_asteroid2.f90  
Asteroids/impactstirring.f90  
Asteroids/common_subs.f90  
Asteroids/sphere1d_implicit.f90: (main program)  
Common/conductionQ.f90  
Common/derivs.f90  
Common/flux_noatm.f90  
Common/generalorbit.f  
Common/grids.f90  
Common/psv.f  
Common/ran2.for  
Common/tridag.for  
*Documentation: Schorghofer, N. (2016) Icarus 276, 88-95*  


### Lunar Ice Pump

Lunar/oscidea1.f90: (main program)  
Lunar/oscidea_subs.f90  
Lunar/subl_h2o.f90  
*Documentation: User Guide Section 3.3  
Documentation: Schorghofer, N. & Aharonson, O. (2014) ApJ 788, 169*  


### Irradiance Model for Terrestrial Analog

EarthAnalogs/insol3d_earth.f90: (main program)  
EarthAnalogs/insol_flat.f90: (simple main program)  
EarthAnalogs/mk_atmosphere.f90  
EarthAnalogs/sunpos.f90  
EarthAnalogs/filemanager.f90  
Common/flux_noatm.f90  
Topo3D/topo3d_modules.f90  
Topo3D/topo3d_common.f90  
Topo3D/topo3d_subs.f90  
*Documentation: User Guide Sections 4.1, 5.1, 5.3-5.7*  


### Terrain Shadowing and Illumination

Topo3D/shadows.f90: (main program)  
Topo3D/fieldofviews.f90: (main program)  
Topo3D/cratersQ_equilbr.f90: (main program)  
Topo3D/cratersQ_moon.f90: (main program)  
Topo3D/cratersQ_mars.f90: (main program)  
Topo3D/cratersQ_mars_parallel.f90: (main program)  
Topo3D/cratersQ_mars_full.f90: (main program)  
Topo3D/insol3d_mars.f90: (main program)  
Topo3D/filemanager.f90  
Topo3D/topo3d_modules.f90   
Topo3D/topo3d_common.f90  
Topo3D/topo3d_geometry.f90  
Topo3D/shadow_subs.f90  
Topo3D/fieldofview_subs.f90  
Topo3D/topo3d_subs.f90  
Topo3D/topo3d_subs_mars.f90  
Topo3D/multigrid.f90  
Topo3D/fieldproperties.f90: (main program)  
Topo3D/hpsort.for  
Topo3D/makegaussian.c   
Common/conductionQ.f90    
Common/conductionQ2.f90  
Common/conductionT2.f90  
Common/flux_noatm.f90  
Common/grids.f90  
Common/julday.for    
Common/tridag.for  
Mars/flux_mars.f90  
Mars/marsorbit.f90  
*Documentation: User Guide Part 5*  


### Monte-Carlo Model for Surface-bounded Exospheres

Exosphere/moon_exo.f90: (main program)  
Exosphere/ceres_exo.f90: (main program)  
Exosphere/body.f90  
Exosphere/montecarlo.f90  
Exosphere/geogrid.f90  
Exosphere/geogrid_D.f90  
Exosphere/subl_subs.f90  
Exosphere/gasdev.for  
Common/ran2.for  
*Documentation: User Guide Part 6*  


---

NOTE: Third party source code from Numerical Recipes is covered by a separate copyright. These are files ending with .for.  A few code snippets from other sources are also used, as documented in the source code.


### ACKNOWLEDGMENTS

2019: Thanks to Sam Potter for comments that helped me speed up the view factor calculations  

Mar 2016: Thanks to Lior Rubanenko for a bug report in cratersQ_*

2010, 2005: Thanks to Oded Aharonson for improvements on mars_mapi2p and a better treatment of the frost/no-frost surface boundary condition.

2006: Troy Hudson discovered a grid-point offset in conductionT and conductionQ, which has been corrected.

2005: Thanks to Mischa Kreslavsky for providing formulas for energy balance on a planar slope.

Many Thanks to Andy Vaught for developing an open-source Fortran 95 compiler (www.g95.org).  The early versions of this code were developed with this compiler.

2001: Samar Khatiwala invented an elegant implementation of the upper radiation boundary condition for the Crank-Nicolson method.

SUPPORT: This code development was supported mainly by NASA, and in smaller parts by Caltech and the University of Hawaii. Undoubtedly, some parts were written in my spare time.

