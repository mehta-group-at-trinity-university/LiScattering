
module units
  ! The following constants give one unit in units of another.
  ! For example KPerInvcm means "Kelvin in units of inverse cm".
  ! The numerical value on the RHS is therefore the value of 1K per 1 cm**(-1)
  ! The values were taken mostly from the NIST CODATA tables circa 2020.
  
  real*8, parameter :: KPerInvcm = 1.4387862961655296d0
  real*8, parameter :: mKPerInvcm = 1.4387862961655296d3
  real*8, parameter :: HartreePerInvcm = 0.0000045563352812122295d0
  real*8, parameter :: eVPerInvcm = 0.00012398425731484318d0
  real*8, parameter :: BohrPerAngstrom = 0.529177210903d0
  real*8, parameter :: AngstromPerBohr = 1d0/BohrPerAngstrom
  real*8, parameter :: THzPerJ = 6.626070150d-22
  real*8, parameter :: THzPerAU = 1.5198298460570d-4
  real*8, parameter :: GaussPerTesla = 1d-4
  real*8, parameter :: hSI = 6.62607015d-34
  real*8, parameter :: HartreePerJoule = 4.3597447222071d-18
  

!!$  real*8, parameter ::  AnstperBohr,HartreeperJoule,JouleperHartree,hSI,hHartreesec,hbarcMeVfm,hbarceVnm
!!$  
!!$  Kperinvcm = 
!!$  mKperinvcm = Kperinvcm*1d3
!!$  Hartreeperinvcm = 0.0000045563352812122295d0
!!$  eVperinvcm = 0.00012398425731484318d0
!!$  BohrperAngstrom = 0.529177210903d0
!!$  Angstromperbohr = 1d0/bohrperAngstrom
!!$  HartreeperJoule = 2.2937122783963d17
!!$  JouleperHartree = 1d0/HartreeperJoule
!!$  
!!$  hHartreeSec = hSI*JouleperHartree

end module units
