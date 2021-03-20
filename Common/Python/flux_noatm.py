import numpy as np
from math import sin, cos, sqrt, acos, pi


def flux_noatm(R,decl,latitude,HA,SlopeAngle,azFac):
#**********************************************************************
#   flux_noatm: calculates incoming solar flux without atmosphere
#     R: distance from sun (AU)
#     decl: planetocentric solar declination (radians)
#     latitude: (radians)
#     HA: hour angle (radians from noon, clockwise)
#     SlopeAngle: >0, (radians) 
#     azFac: azimuth of topographic gradient (radians east of north)
#            azFac=0 is south-facing  
#**********************************************************************
    So=1365.  # solar constant [W/m^2]
  
    c1 = cos(latitude)*cos(decl)
    s1 = sin(latitude)*sin(decl)
    # beta = 90 minus incidence angle for horizontal surface
    # beta = elevation of sun above (horizontal) horizon 
    sinbeta = c1*cos(HA) + s1
  
    cosbeta = sqrt(1-sinbeta**2)
    # ha -> az (option 2)
    buf = (sin(decl)-sin(latitude)*sinbeta)/(cos(latitude)*cosbeta)
    # buf can be NaN if cosbeta = 0
    if buf>+1.:
        buf=+1.0  # roundoff
    if buf<-1.:
        buf=-1.0  # roundoff
    azSun = acos(buf)
    if sin(HA)>=0:
        azSun=2*pi-azSun

    # theta = 90 minus incidence angle for sloped surface
    sintheta = cos(SlopeAngle)*sinbeta - \
        sin(SlopeAngle)*cosbeta*cos(azSun-azFac)
    if cosbeta==0.:
        sintheta = cos(SlopeAngle)*sinbeta
    if sintheta<0.:
        sintheta = 0. # horizon
    if sinbeta<0.:
        sintheta=0.  # horizontal horizon at infinity
  
    flux_noatm = sintheta*So/(R**2)
  
    return flux_noatm

