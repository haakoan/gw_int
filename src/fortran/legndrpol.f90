!This module calculate Spherical harmonics and
!Associated Legendre polynomials. It is build heavly
!on num.rec rutines and should be used with some care 
!if licences is an issue.
!TODO: Rewrite code and get rid of num.rec. Replace with
!own code.
!Haakon Andresen 17.11.2016

module sph_leg_pol
!Code interface  
  interface sphleg_ifr
     !module procedure spharm
     module procedure plgndr
     module procedure spharm_r
  end interface sphleg_ifr
CONTAINS

! FUNCTION spharm(l,m,th,ph)    
! !Rreturns the spherical harmonics Y(l,m,theta,phi)
! !Calls the functions plgndr and frc.
!     INTEGER :: l,m
!     real(kind=8) :: th,ph,tp
!     complex :: spharm
!     real(kind=8) :: pc_pi
!     parameter(pc_pi = dacos(-1d0))
!     tp = plgndr(l,m,cos(th)) &
!          *dsqrt( (2d0*dble(l)+1d0)/(4d0*pc_pi) &
!          * frc(l-m)/frc(l+m))
!     spharm = CMPLX(dcos(DBLE(m)*ph)*tp,dsin(DBLE(m)*ph)*tp)
!     RETURN
! END FUNCTION

FUNCTION spharm_r(l,m,th,ph)    
!Rreturns the spherical harmonics RE(Y(l,m,theta,phi))
!Calls the functions plgndr and frc.
    INTEGER :: l,m
    real(kind=8) :: th,ph,tp
    real(kind=8) :: spharm_r
    real(kind=8) :: pc_pi
    parameter(pc_pi = dacos(-1d0))
    tp = plgndr(l,m,dcos(th))*sqrt(dble(2*l+1)/(4.0d0*pc_pi))*&
                              sqrt(frc(l-m)/frc(l+m)) 
    spharm_r = cos(DBLE(m)*ph)*tp
    RETURN
END FUNCTION



recursive FUNCTION frc(l) result(a)    
!Recursive function to calculate factorials (n!)
    INTEGER :: l,m
    real(kind=8) :: a
    if(l .eq. 0) then
       a = 1d0
       RETURN
    endif
    a = dble(l)*frc(l-1)
    return
END FUNCTION

recursive FUNCTION plgndr(l,m,x) result(plx)
!From num.rec. Calculates Associated Legendre polynomials.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l,m
  real(kind=8), INTENT(IN) :: x
  REAL(kind=8) :: plx
  INTEGER :: ll
  REAL(kind=8) :: pll,pmm,pmmp1,somx2
  pmm=1.0
  if(m < 0) then
     plx = ((-1.0d0)**abs(m))*plgndr(l,abs(m),x)*frc(l-abs(m))/frc(l+abs(m))
     return
  end if
  if (m > 0) then
     somx2=dsqrt((1.0d0-x)*(1.0d0+x))
     pmm=product(arth(1.0d0,2.0d0,m))*somx2**m
     if (mod(m,2) == 1) pmm=-pmm
  end if
  if (l == m) then
     plx=pmm
  else
     pmmp1=x*(2*m+1)*pmm
     if (l == m+1) then
        plx=pmmp1
     else
        do ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        end do
        plx=pll
     end if
  end if
END FUNCTION plgndr


FUNCTION arth(first,increment,n)
!From num.rec. Calculates something.
  REAL(kind=8), INTENT(IN) :: first,increment
  INTEGER, INTENT(IN) :: n
  REAL(kind=8), DIMENSION(n) :: arth
  INTEGER :: k,k2
  REAL(kind=8) :: temp
  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  if (n > 0) arth(1)=first
  if (n <= NPAR_ARTH) then
     do k=2,n
        arth(k)=arth(k-1)+increment
     end do
  else
     do k=2,NPAR2_ARTH
        arth(k)=arth(k-1)+increment
     end do
     temp=increment*NPAR2_ARTH
     k=NPAR2_ARTH
     do
        if (k >= n) exit
        k2=k+k
        arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
        temp=temp+temp
        k=k2
     end do
  end if
END FUNCTION arth
end module sph_leg_pol
