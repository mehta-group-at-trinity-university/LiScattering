
module units
  ! The following constants give one unit in units of another.
  ! For example KPerInvcm means "Kelvin in units of inverse cm".
  ! The numerical value on the RHS is therefore the value of 1K per 1 cm**(-1)
  ! The values were taken mostly from the NIST CODATA tables circa 2020.

  real*8, parameter :: Pi = 3.14159265358979323846264338328
  real*8, parameter :: KPerInvcm = 1.4387862961655296d0
  real*8, parameter :: mKPerInvcm = 1.4387862961655296d3
  real*8, parameter :: InvcmPerHartree = 4.556335252912088d-6
  real*8, parameter :: HartreePerInvcm = 219474.6313632043
  real*8, parameter :: InvcmPereV = 0.00012398425731484318d0
  real*8, parameter :: BohrPerAngstrom = 0.529177210903d0
  real*8, parameter :: AngstromPerBohr = 1d0/BohrPerAngstrom
  real*8, parameter :: THzPerJ = 6.626070150d-22
  real*8, parameter :: THzPerHartree = 1.5198298460570d-4
  real*8, parameter :: GaussPerTesla = 1d-4
  real*8, parameter :: hSI = 6.62607015d-34
  real*8, parameter :: HartreePerJoule = 4.3597447222071d-18
  real*8, parameter :: HartreePerMHz = 6.579683920502d9
  real*8, parameter :: MHzPerHartree = 1d0/HartreePerMHz
  real*8, parameter :: HartreePerTHz = 6.579683920502d3
  real*8, parameter :: amuKg = 1.66053906660d-27
  real*8, parameter :: ElectronMassKg = 9.109383701528d-31
  real*8, parameter :: amuAU = amuKg/ElectronMassKg
  real*8, parameter :: vdwbetaLi = 65.2049d0
  real*8, parameter :: nKPerHartree = 3.1668115634556D-15
  real*8, parameter :: HartreePernK = 1/nKPerHartree
  real*8, parameter :: HartreePermK = HartreePernK*1d-6
  real*8, parameter :: Bohrpercm = 5.2917724900001d-9

  

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
