module TwoBodydata
  implicit none
  integer MatrixDim,LegPoints,xDim,xNumPoints,Order,Left,Right
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision mu,xMin,xMax,lam,r2b,DD,x0,alpha
  double precision, allocatable :: u(:,:,:),ux(:,:,:),uxx(:,:,:),xPoints(:),xIntPoints(:),VR(:)
  double precision, allocatable :: S(:,:),L(:,:),H0(:,:),G0(:,:),V0(:,:),delta(:) 
  integer, allocatable :: xBounds(:)   
  character(LEN=30) :: Format1
     
  
end module TwoBodydata
