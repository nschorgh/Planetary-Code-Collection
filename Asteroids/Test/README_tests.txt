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

compare with depths_fast2.0


4. Run sphere1d_implicit

make sphere1d
./a.out

compare output with depths_sphere.dat


5. Run trojans_fast3

make trojans_fast3
./a.out 3

compare output depths.3 with depths_trojan.3
(takes about 10 minutes to finish)
