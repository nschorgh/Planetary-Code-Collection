Planetary Code Collection
=========================

*by Norbert Schorghofer*


This program collection contains:

* Semi-implicit one-dimensional thermal model for planetary surfaces (Crank-Nicolson solver with Stefan-Boltzmann Law boundary condition)  
* 3D model of direct insolation, terrain shadowing, and terrain irradiance for airless bodies, Mars, and Mauna Kea  
* Monte-Carlo model for ballistic hops in the exospheres of the Moon and Ceres
* Miscellaneous
  * insolation model for Mauna Kea
  * asynchronously coupled method for ice retreat on asteroids    
  * lunar ice pump (in boundary value formation)  
  * vapor diffusion at lunar conditions (microphysical model)  


The theory behind the numerical methods is described in [UserGuide.pdf](./UserGuide.pdf) or in individual journal articles cited below.


Mars Subsurface Ice Model (MSIM)
--------------------------------

__NOTE:__ MSIM now has its own repository at https://github.com/nschorgh/MSIM and is no longer part of the Planetary-Code-Collection. 

In brief, MSIM contains:

* Mars Thermal Model
* Vapor Diffusion Model for Mars
* Equilibrium Ice Table
* Fast (asynchronously-coupled) Method for Subsurface Ice Dynamics


Other Models for Planetary Surfaces
-----------------------------------

### Basic Thermal Model for Airless Bodies

Standard thermal model for asteroidal surfaces. The one-dimensional heat equation is solved, based on solar energy input and thermal radiation from the surface to space.  The solver for the one-dimensional heat equation is semi-implicit, with the Stefan-Boltzmann radiation law as upper boundary condition.  The finite-difference method is flux-conservative even on an irregularly spaced vertical grid and the thermal properties of the soil can vary spatially and with time.  (This is a simplified version of the Mars Thermal Model. They both use the same flux-conservative spatial discretization and the same Crank-Nicolson solver with nonlinear boundary condition.)  

Directory: `Asteroids/`  
The core subroutines for the thermal model are located in `Common/` and available in Fortran, C, Matlab, and Python.  
*Documentation: User Guide Part 1*  


### Asynchronous Models for Temperature, Impact Stirring, and Ice Loss on Asteroids

This set of models for asteroidal surfaces combines diurnally-resolved temperatures, the long-term loss of near-surface ice to space, and probabilistic impact stirring (only one-dimensional). It can be used to estimate long-term ice loss from near the surface due to sublimation as the orbit is changing, the Sun brightens, or ice is mixed due to impacts. 
A significant complexity in this model arises from partially ice-filled pore spaces (necessary to incorporate the consequences of impact stirring), because the re-distribution of ice within the pores due to vapor diffusion and deposition adds another partial differential equation. A simpler two-layer version, where pore spaces are either empty or full, is also implemented.  

Directory: `Asteroids/`  
*Documentation: [Schorghofer (2016)](https://doi.org/10.1016/j.icarus.2016.04.037)*  


### Lunar Ice Pump

Two types of models for water vapor migration in the shallow lunar subsurface are implemented. One is a particle-based microphysical model that follows water molecules that undergo a random walk. The other uses the time-averaged boundary value formulation for the same process. Both models are one-dimensional, and primarily intended for the low temperature regime, when molecular adsorption times are long.

Directory: `Lunar/`  
*Documentation for boundary-value model: [Schorghofer & Aharonson (2014)](https://doi.org/10.1088/0004-637X/788/2/169)*  
No documentation has yet been written for the random walk model.  


### Irradiance Model for Terrestrial Analog

Clear-sky direct and indirect short-wave irradiance on Mauna Kea summit, based on optical path length but an otherwise 0-dimensional atmospheric model.
The incoming irradiance can be calculated for a flat unobstructed surface, but also for a tilted suface with horizons, i.e., 3D sky irradiance.  

Directory: `EarthAnalogs/`  
*Documentation: User Guide Sections 3.1 and 2.5*  


### Terrain Shadowing and 3D Surface Energy Balance 

This model of the three-dimensional surface energy balance calculates horizons from a digital elevation model for terrain shadowing calculations and, optionally, also view factors for use in terrain irradiance calculations. The surface energy balance model can then be coupled to the model of subsurface heat conduction introduced above. This then provides a complete thermal model for  rugged terrain on Mars or airless bodies. The model is still at prototype stage, but has been used in several research studies.

Directory: `Topo3D/`  
*Documentation: User Guide Part 2*  


### Monte-Carlo Model for Surface-bounded Exospheres

Ballistic trajectories of neutral molecules or atoms are modeled for large airless bodies (the Moon and Ceres).  Individual molecules are launched with probabilistically distributed velocities. The model then computes the molecule's impact location and time analytically, using a closed-form solution for the intersection of an ellipse with a sphere, i.e., without numerical integration of the particle trajectory.
An event-driven algorithm processes landing and launching events in time-order.
Molecules may be lost or destroyed by photo-dissociation before they land.
Surface temperatures are based on the thermal model for airless bodies.  

Directory: `Exospheres/`  
*Documentation: User Guide Part 4*  


---

### Notes

Third party source code from Numerical Recipes is covered by a separate copyright. These are files ending with .for.  A few code snippets from other sources are also used, as documented in the source code.


### Acknowledgments

2019: Thanks to Sam Potter for comments that helped me speed up the view factor calculations  

2016: Thanks to Lior Rubanenko for a bug report in cratersQ_*

2006: Troy Hudson discovered a grid-point offset in conductionT and conductionQ, which has been corrected.

Many Thanks to Andy Vaught for developing an open-source Fortran 95 compiler (G95).  The early versions of this code were developed with this compiler.

2001: Samar Khatiwala invented an elegant implementation of the upper radiation boundary condition for the Crank-Nicolson method.

SUPPORT: This code development was supported mainly by NASA, and in smaller parts by Caltech and the University of Hawaii. Undoubtedly, some parts were written in my spare time.

