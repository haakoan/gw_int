module data_types

!This is the main construct used
! to store the data that we read from
!files. Most subrutines
!take this type as input.
TYPE fruit
 REAL(KIND=8) :: time
 REAL(KIND=8), dimension(:),ALLOCATABLE :: xzn(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: xzr(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: xzl(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: yzn(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: yzl(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: yzr(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: zzn(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: zzr(:)
 REAL(KIND=8), dimension(:),ALLOCATABLE :: zzl(:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: den(:,:,:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: gpo(:,:,:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: vex(:,:,:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: vey(:,:,:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: vez(:,:,:)
 REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: yyw(:,:)
END TYPE fruit

TYPE cell
   REAL(KIND=8) :: ct(2,2,2)
   REAL(KIND=8) :: cn(2)
   INTEGER :: ny,nz
END TYPE cell
TYPE bnd
   REAL(KIND=8) :: ct(2,2)
   INTEGER :: ny,nz
END TYPE bnd


integer lmax
parameter(lmax=2)
end module data_types
