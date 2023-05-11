Validation and Example Results
==============================

1. Test of a semi-implicit heat equation solver with prescribed surface temperature
(conductionT.f90)

make testcrankT
a.out

- Tprofile_testcrankT contains 12 temperature profiles from the last sol
- z.testcrankT contains the coordinates of the vertical grid
- test_Tprofile.m reads these profiles and compares them to the analytical solution

The result is shown in Figure 1.1 of the UserGuide.pdf

Note: The analytical solution is for an infinitely thick domain, so the small
      deviations toward the bottom boundary are okay.



2. Example solution of semi-implicit heat equation solver with Stefan-Boltzmann Law
boundary condition (conductionQ.f90 or conductionQ2.f90)

make testcrankQ
a.out 

- Tprofile_testcrankQ contains 12 temperature profiles from the last sol
- z.testcrankQ contains the coordinates of the vertical grid
- stdout_of_testcrankQ.txt contains the time-averaged temperature profile and heat flux
After the temperature has equilibrated, the time-averaged heat flux must be
constant with depth and must equal the heat flux specified at the bottom boundary.



3. Test implementation of nonlinear boundary condition in Crank-Nicolson solver
(conductionQ.f90)

make testcrankQ_asymp
a.out

The numerical solution is in the 2nd column of the output file 'Tsurface' and
the analytical solution for small times is in the 4th column. The result is shown
in Figure 1.4 of the UserGuide.pdf.



4. Test rate of convergence of conductionQ.f90 with time step

make testcrankQ_conv
a.out

The output file 'Tprofiles', contains temperature profiles for eight different
values for the time step. Errors can be defined by differences between these profiles.
