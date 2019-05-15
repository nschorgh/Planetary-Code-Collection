Validation and Example Result
=============================

Bowl-shaped crater with depth-to-diameter ratio of 1:5
Domain with 81x81 pixels, 5m horizontal resolution, sphere diameter 50m
(created with program maketopo.f90)


PREP STEP: PRE-CALCULATE GEOMETRIC INFORMATION

Two variants were run:
* shadows.f90 with 180 azimuth rays, takes 22 seconds with serial implementation
* fieldofviews.f90, takes 1.5 minutes

Input: topo81.xyz
Output: horizons.topo81  (180 azimuth rays)
Output: fieldofviews.topo81 


RUN 1: COMPARISON WITH INGERSOLL SOLUTION

cratersQ_snaphost.f90 calculates instantaneous solution, takes 30 seconds
albedo 0.1, emissivity 0.95
sun elevation 10 degree, solar constant 1365 W/m^2

Output: qinst.topo81 (meaning of columns is in write(21,... )
Used qinst.m to make figure qinst.png (insolation and temperature)
Used qtest_ingersoll.m to make figure qtest_ingersoll.eps


RUN 2: ONE SOLAR DAY

cratersQ_moon.f90 calculates one solar day, takes 5.5 minutes
latitude 80 degree, declination 0 degree (hence 10 degree maximum solar elevation)
time step 1/100th of solar day

Output: qmean.topo81 (meaning of columns is in write(22,... )
Used qmean1.m to make figure qmean.png


RUN 3: WITH SUBSURFACE CONDUCTION

cratersQ_moon.f90 as Run 2, but with subsurface conduction, takes 1 hour
domain depth = 5 skin depths, 10+1 solar days, thermal inertia 100 tiu

Output: qmean.topo81_wsurfcond (meaning of columns is in write(22,... )
Used qmean_compare.m to make figure qmean_compare.eps




