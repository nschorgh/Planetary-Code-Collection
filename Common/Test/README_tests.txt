Validation and Example Results
==============================

1. Test of a semi-implicit heat equation solver with prescribed surface temperature
conductionT.f90

make testcrankT
a.out

- Tprofile_testcrankT contains 12 temperature profiles from the last sol
- z.testcrankT contains the coordinates of the vertical grid
- test_Tprofile.m reads these profiles and compares them to the analytical solution

Note: The analytical solution is for an infinitely thick domain, so the small deviations toward the bottom boundary are okay.



2. Example solution of semi-implicit heat equation solver with thermal emission boundary condition
conductionQ.f90 or conductionQ2.f90
Crank-Nicoloson solver with nonlinear upper boundary condition

make testcrankQ
a.out 

- Tprofile_testcrankQ contains 12 temperature profiles from the last sol
- z.testcrankQ contains the coordinates of the vertical grid
- stdout_of_testcrankQ.txt contains the time-averaged temperature profile and heat flux
After the temperature has equilibrated, the time-averaged heat flux must be constant with depth and must equal the heat flux specified at the bottom boundary.
