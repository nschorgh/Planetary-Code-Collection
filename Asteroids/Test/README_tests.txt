Validation and Example Results
==============================

1. Run basic 1D thermal model for airless body

make asteroid_thermal
./asteroid_thermal  

compare output with content of file 'out'


2. Run asteroid_fast1

make asteroid_fast1
./fast 0
(This might take 1 hour to finish, because of the long integration time.)

compare with depths_fast1.0


3. Run asteroid_fast2

make asteroid_fast2
./fast 0
(This might take 1.5 hours to finish, because of the long integration time.)

compare with depths.0


4. Run sphere1d_implicit

make sphere1d
./a.out

compare output with depths_sphere.dat
