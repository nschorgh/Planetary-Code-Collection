module allinterfaces
  
  ! Common/*.f90
  interface
     pure function flux_noatm(R,decl,latitude,HA,SlopeAngle,azFac)
       implicit none
       real(8) flux_noatm
       real(8), intent(IN) :: R,decl,latitude,HA,SlopeAngle,azFac
     end function flux_noatm
  end interface

  interface
     function colint(y,z,nz,i1,i2)
       implicit none
       integer, intent(IN) :: nz, i1, i2
       real(8), intent(IN) :: y(nz), z(nz)
       real(8) colint
     end function colint
  end interface

  ! common_subs3.f90
  interface
     pure function flux2T(Q,albedo,emiss)
       implicit none
       real(8) flux2T
       real(8), intent(IN) :: Q, emiss, albedo
     end function flux2T
  end interface

  interface
     pure function a2Torb(semia)
       implicit none
       real(8), intent(IN) :: semia  ! semimajor axis [AU]
       real(8) a2Torb  
     end function a2Torb
  end interface

  interface
     pure function sols_per_orbit(semia,solarDay)
       implicit none
       real(8) sols_per_orbit
       real(8), intent(IN) :: semia, solarDay
     end function sols_per_orbit
  end interface

  interface
     pure function faintsun(t)
       implicit none
       real(8) faintsun
       real(8), intent(IN) :: t   ! time before present [years]
     end function faintsun
  end interface
    
  ! soilthprop_trojan.f90
  interface
     pure function heatcapacity(T)
       implicit none
       real(8), intent(IN) :: T  ! [K]
       real(8) heatcapacity   ! [J/(kg K)]
     end function heatcapacity
  end interface

  interface
     elemental function meanfreepathinsoil(diam,porosity)
       implicit none
       real(8) meanfreepathinsoil ! ell
       real(8), intent(IN) :: diam, porosity
     end function meanfreepathinsoil
  end interface

  interface
     subroutine assignthermalproperties3(nz,z,T,porosity,ti,rhocv,icefrac,zdepthT)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), T(nz), porosity(nz)
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
       real(8), intent(IN), optional :: icefrac(nz), zdepthT
       real(8), external :: heatcapacity
     end subroutine assignthermalproperties3
  end interface

  ! fast_subs_airlessbody3.f90
  interface
     subroutine icechanges3(nz,z,avSice,elleff,bigstep,zdepthT,porosity,icefrac)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), avSice, elleff, bigstep
       real(8), intent(IN) :: porosity(nz), icefrac(nz)
       real(8), intent(INOUT) :: zdepthT
     end subroutine icechanges3
  end interface

  interface
     pure function getfirst(zdepth,nz,z)
       implicit none
       integer getfirst
       integer, intent(IN) :: nz
       real(8), intent(IN) :: zdepth, z(nz)
     end function getfirst
  end interface

  ! other
  interface
     real(8) function psv(T)
       implicit none
       real(8), intent(IN) :: T
     end function psv
  end interface
end module allinterfaces


