Mars surface temperatures over time
===================================

Mean annual surface temperature on Mars over the past 21 Myr, calculated by tempr_driver.f90

These model outputs correspond to Figure 3 in Schorghofer, Geophys. Res. Lett. 35, L18201 (2008), and are here provided in electronic form for a longer time period.

Inputs are lats.21Myr (latitude, albedo, thermal inertia) and orbital elements according to Laskar et al., Icarus 170, 343 (2004).

Outputs:

   info_geo.21Myr shows general model parameters for each latitude.

   Surface temperatures are stored in out_subsurf.21Myr
   Columns contain:
      time before present (Earth years)
      latitude (deg)
      mean annual surface temperature (K)
      temperature at bottom of the computational domain (K)

   out_geo.21Myr repeats the orbital elements (obliquity, eccentricity, and longitude of perihelion from moving equinox)


History
-------
Original calculations over 5 Myr, Jan 2008
Recalculated over 21 Myr and archived in electronic form, May 2020

Norbert Schorghofer
