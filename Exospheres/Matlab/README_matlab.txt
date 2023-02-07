
This directory contains a simplified Matlab/Octave version of the exosphere model.
The model follows ballistic trajectories with an event-driven algorithm.
The driving program is 'ceres_exo.m'. All other '.m' files contain functions.

Compared to the Fortran version, the surface temperature model is much simplified and photodissociation is limited to on/off (no night-time shadowing), which suffices for Ceres, where gravitational escape dominates.

It works with Matlab and Octave.
