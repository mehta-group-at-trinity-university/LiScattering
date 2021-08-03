module EREdata
  implicit none
  integer MatrixDim,LegPoints,xDim,xNumPoints,Order,Left,Right
  integer PotInterpPoints,PotInterpOrder
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision, allocatable :: VS(:),VT(:),R(:),Vbcoef(:),Potknots(:)
  double precision mu,muref,xMin,xMax,lam,a2binv,r2b,DD,x0,alpha,leff
  double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),xPoints(:)
  double precision, allocatable :: S(:,:),L(:,:),H0(:,:),G0(:,:),V0(:,:),delta(:) 
  integer, allocatable :: xBounds(:)   
  character(LEN=30) :: Format1
     
  
end module EREdata
