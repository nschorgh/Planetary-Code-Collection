function z = setgrid(nz,zmax,zfac)
  % construct regularly or geometrically spaced 1D grid
  % z(n)-z(1) = 2*z(1)*(zfac**(n-1)-1)/(zfac-1)
  % choice of z(1) and z(2) is compatible with conductionQ

  dz = zmax/nz;
  z = [0.5:nz]' * dz;
  
  if (zfac>1.), 
    dz = zmax/(3.+2.*zfac*(zfac^(nz-2)-1.)/(zfac-1.));
    z(1) = dz;
    z(2) = 3*z(1);
    for i=3:nz
      z(i) = (1.+zfac)*z(i-1) - zfac*z(i-2);
      % z(i) = z(i-1) + zfac*(z(i-1)-z(i-2)); ! equivalent
    end
  end
  
  return


