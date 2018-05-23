module sphexp
  use sph_leg_pol
  use data_types
  use quads
  use dataread
  use tblinf
contains
 
 subroutine run_int_test
   LOGICAL :: loop
   REAL(KIND=8),dimension(:,:,:,:),ALLOCATABLE :: data(:,:,:,:)
   CHARACTER(LEN=90) :: filename,filein
   INTEGER(hsize_t) ::ij,s,ik,sts
   TYPE(fruit) :: modlya,modlyi
   REAL(KIND=8), dimension(6) :: am1,am2
   REAL(KIND=8) :: akp
   real :: r
   integer :: i,j,reason,NstationFiles
   character(LEN=40), dimension(:), allocatable :: stationFileNames
   2000 format(7(1E20.10))


   ik =1
   !ALLOCATE(cutm, STAT = sts)
   CALL getarg(1, filein)
   filename = 'sph_an'//trim(filein)//'.dat' 	
   open(332,FILE=trim(filename),action="write",position='append')
   filename = 'gw_rot'//trim(filein)//'.dat'
   open(333,FILE=trim(filename),action="write",position='append')

   ! do i = 1,1!NstationFiles

   filename = './'//trim(filein)//'/'//trim(filein) 
   ij = 0
   loop = .TRUE.
   write(*,*) filename
   do while(loop)
      CALL setdata(filename,ij,modlya,s,0)
      !CALL setdata(filename,ij,modlyi,s,1)
      ! CALL setdata(filename,ij,modlyi,s,0)
      ij = ij+1
      !if( ij .EQ. s) then
      loop = .FALSE. 
      !endif
      CALL integral_test(modlya,am1,am2)      
      write(332,2000) modlya%time-0.24892d0,am2
      write(333,2000) modlya%time-0.24892d0,am1
      write(*,2000) 1.0d0,am2
      write(*,2000) 1.0d0,am1
    
      !loop = .FALSE. 
   enddo
   !end do
   close(332)
   close(333)
 end subroutine run_int_test


 subroutine integral_test(modlya,amp1,amp2,modlyi)
   USE sph_leg_pol
   implicit none
   TYPE(fruit), INTENT(IN) :: modlya
   TYPE(fruit),OPTIONAL, INTENT(IN) :: modlyi 
   TYPE(fruit) :: yat,yit
   REAL(KIND=8), DIMENSION(6), INTENT(OUT) :: amp1,amp2 
   REAL(KIND=8) :: pc_gc
   parameter(pc_gc = 6.67259d-08)
   REAL(KIND=8) :: pc_pi
   parameter(pc_pi = dacos(-1d0))
   REAL(KIND=8) :: pc_cl
   parameter(pc_cl = 2.99792458d+10)
   REAL(KIND=8) :: geofac
   parameter(geofac = 1.d0/3.d0)

   REAL(KIND=8) :: x,y,z,rsin,dmas,ph,th,ph1,th1,sph1,sph2,c2,gc2vf
   INTEGER :: nx,ny,nz,ny2,nz2,j,k,ast,l,m,lm,i,js,jks
   REAL(KIND=8) :: qxx,qyy,qzz,qxy,qxz,qyz
   REAL(KIND=8) :: axx,ayy,azz,axy,axz,ayz,dummy

   REAL(KIND=8), dimension(:),ALLOCATABLE :: dph(:),dth(:),dr(:),cth,sth,cph,sph
   REAL(KIND=8) :: corr_sph,corr_cth,corr_sth,deltath,deltaph
   REAL(KIND=8), dimension(:,:),ALLOCATABLE :: vxcofs,vycofs,vzcofs,dencofs
   REAL(KIND=8), dimension(:,:,:),ALLOCATABLE,target :: ud(:,:,:)
   REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: dphidr(:,:,:) 
   REAL(KIND=8), DIMENSION(6) :: amplya,amplyi 
   gc2vf = pc_gc/pc_cl**2d0
   c2 = pc_cl**2d0

   nx=size(modlya%den(:,1,1))
   ny=size(modlya%den(1,:,1))
   nz=size(modlya%den(1,1,:))

   ALLOCATE (dth(ny), STAT = ast)
   ALLOCATE (dph(nz), STAT = ast)
   dth = dcos(modlya%yzl(:))-dcos(modlya%yzr(:))
   dph = modlya%zzr(:)-modlya%zzl(:)

   ! corr_cth = 0.99995140318460529d0! dsqrt(2d0/3d0/corr_cth)
   ! corr_sth = 1.000024301064834d0!dsqrt(4d0/3d0/corr_sth)
   ! corr_sph = 0.99999978633434816d0 !Ewalds factor



   corr_sph = SUM(dsin(modlya%zzn(:))**2d0*dph(:))
   corr_cth = SUM(dcos(modlya%yzn(:))**2d0*dth(:))
   corr_sth = SUM(dsin(modlya%yzn(:))**2d0*dth(:))
   corr_sph = dsqrt(pc_pi/corr_sph)
   corr_cth = dsqrt(2d0/3d0/corr_cth)
   corr_sth = dsqrt(4d0/3d0/corr_sth) 

   ALLOCATE (sph(nz), STAT = ast)
   ALLOCATE (cph(nz), STAT = ast)

   ALLOCATE (cth(ny), STAT = ast)
   ALLOCATE (sth(ny), STAT = ast)

   do j = 1,ny
      cth(j) = dcos(modlya%yzn(j))*corr_cth
      sth(j) = dsin(modlya%yzn(j))*corr_sth
   end do
   do k =1,nz
      cph(k) = dcos(modlya%zzn(k))
      sph(k) = dsin(modlya%zzn(k))*corr_sph
   end do
   ALLOCATE(dr(nx), STAT = ast)
   dr = modlya%xzr(:)**3-modlya%xzl(:)**3
   ALLOCATE(vxcofs(0:lmax,-lmax:lmax), STAT = ast)
   vxcofs = 0.0d0
   ALLOCATE(vycofs(0:lmax,-lmax:lmax), STAT = ast)
   vycofs = 0.0d0
   ALLOCATE(vzcofs(0:lmax,-lmax:lmax), STAT = ast)
   vzcofs = 0.0d0
   ALLOCATE(dencofs(0:lmax,-lmax:lmax), STAT = ast)
   dencofs = 0.0d0

   ALLOCATE (dphidr(nx,ny,nz), STAT = ast)
   do i=1,nx-1
      dphidr(i,:,:) = (modlya%gpo(i+1,1:ny,:)-modlya%gpo(i,1:ny,:))/(modlya%xzr(i)-modlya%xzl(i))
   enddo
   dphidr(nx,:,:) = (modlya%gpo(nx,1:ny,:)-modlya%gpo(nx-1,1:ny,:))/(modlya%xzr(nx)-modlya%xzl(nx))
   amp2 = 0.0d0
   axx = 0.0d0; ayy = 0.0d0; azz = 0.0d0;
   axy = 0.0d0; axz = 0.0d0; ayz = 0.0d0;
   qxx = 0.0d0; qyy = 0.0d0; qzz = 0.0d0;
   qxy = 0.0d0; qxz = 0.0d0; qyz = 0.0d0;

   if(PRESENT(modlyi)) then
      yat = modlya
      yit = modlyi
      !$OMP PARALLEL DO &
      !$OMP & PRIVATE(i,l,m,k,j,qxx,qyy,qzz,qxy,qxz,qyz,z,x,y,th,ph,th1,ph1,dmas,rsin,sph1,sph2) &
      !$OMP & SHARED(modlya,modlyi,dphidr, dr,cph,sph,dth,dph) &
      !$OMP & REDUCTION(+:axx,ayy,azz,axy,axz,ayz,vxcofs,vycofs,vzcofs,dencofs)
      do i = 1,nx
         vxcofs = 0.0d0; vzcofs = 0.0d0; vzcofs = 0.0d0
         dencofs = 0.0d0
         qxx = 0.0d0; qyy = 0.0d0; qzz = 0.0d0;
         qxy = 0.0d0; qxz = 0.0d0; qyz = 0.0d0;
         do k=1,nz
            do j=1,ny
               rsin = sth(j)
               z = cth(j)
               x = rsin * cph(k)
               y = rsin * sph(k)
               th = acos(z)
               ph = atan2(y,x)
               th1 = acos(y)
               ph1 = atan2(z,-x)
               dmas = dph(k)*dth(j)*modlya%yyw(j,k)
               do l = 0,lmax
                  do m = -l,l
                     sph1 = (spharm_r(l,m,th,ph))*(-1.d0**dble(abs(m)))/(sqrt(4.d0*pc_pi*dble(2*l+1)))
                     sph2 = (spharm_r(l,m,th1,ph1))*(-1.d0**dble(abs(m)))/(sqrt(4.d0*pc_pi*dble(2*l+1)))
                     vxcofs(l,m) = vxcofs(l,m) + (modlya%vex(i,j,k)*sph1 + modlyi%vex(i,j,k)*sph2)*dmas
                     vycofs(l,m) = vycofs(l,m) + (modlya%vey(i,j,k)*sph1 + modlyi%vey(i,j,k)*sph2)*dmas
                     vzcofs(l,m) = vzcofs(l,m) + (modlya%vez(i,j,k)*sph1 + modlyi%vez(i,j,k)*sph2)*dmas
                     dencofs(l,m) = dencofs(l,m) + (modlya%den(i,j,k)*sph1 + modlyi%den(i,j,k)*sph2)*dmas
                  enddo
               enddo
            end do
         end do
         do k=1,nz
            do j=1,ny
               rsin = sth(j)
               z = cth(j)
               x = rsin * cph(k)
               y = rsin * sph(k)
               th = acos(z)
               ph = atan2(y,x)
               th1 = acos(y)
               ph1 = atan2(z,-x)

               yat%vex(i,j,k) = epx(th,ph,vxcofs)
               yit%vex(i,j,k) = epx(th1,ph1,vxcofs)

               yat%vey(i,j,k) = epx(th,ph,vycofs)
               yit%vey(i,j,k) = epx(th1,ph1,vycofs)

               yat%vez(i,j,k) = epx(th,ph,vzcofs)
               yit%vez(i,j,k) = epx(th1,ph1,vzcofs)

               yat%den(i,j,k) = epx(th,ph,dencofs)
               yit%den(i,j,k) = epx(th1,ph1,dencofs)
            end do
         end do
         call calc_quads_new_loop(qxx,qyy,qzz,qxy,qxz,qyz,vxcofs,vycofs,vzcofs,dencofs,modlya%xzn(i),dphidr(i,1,1))
         axx = axx + qxx*dr(i); ayy = ayy + qyy*dr(i); azz = azz + qzz*dr(i);
         axy = axy + qxy*dr(i); axz = axz + qxz*dr(i); ayz = ayz + qyz*dr(i);
      end do
      !$OMP END PARALLEL DO
   else
      !dphidr = 0.0d0
      yat = modlya
      !$OMP PARALLEL DO &
      !$OMP & PRIVATE(i,l,m,k,j,qxx,qyy,qzz,qxy,qxz,qyz,z,x,y,th,ph,th1,ph1,dmas,rsin,sph1,sph2) &
      !$OMP & SHARED(modlya,modlyi,dphidr, dr,cph,sph,dth,dph) &
      !$OMP & REDUCTION(+:axx,ayy,azz,axy,axz,ayz,vxcofs,vycofs,vzcofs,dencofs)
      do i = 100,100!nx
         vxcofs = 0.0d0; vzcofs = 0.0d0; vzcofs = 0.0d0
         dencofs = 0.0d0
         qxx = 0.0d0; qyy = 0.0d0; qzz = 0.0d0;
         qxy = 0.0d0; qxz = 0.0d0; qyz = 0.0d0;
         do k=1,nz
            do j=1,ny
               rsin = sth(j)
               z = cth(j)
               x = rsin * cph(k)
               y = rsin * sph(k)
               th = acos(z)
               ph = atan2(y,x)
               dmas = dph(k)*dth(j)!*modlya%yyw(j,k)
               do l = 0,lmax
                  do m = -l,l
                     sph1 = (spharm_r(l,m,th,ph))*((-1.d0)**dble(abs(m)))/(sqrt(4.d0*pc_pi*dble(2*l+1)))
                     vxcofs(l,m) = vxcofs(l,m) + (modlya%vex(i,j,k)*sph1)*dmas!*cos(ph)
                     vycofs(l,m) = vycofs(l,m) + (modlya%vey(i,j,k)*sph1)*dmas
                     vzcofs(l,m) = vzcofs(l,m) + (modlya%vez(i,j,k)*sph1)*dmas
                     dencofs(l,m) = dencofs(l,m) + (modlya%den(i,j,k)*sph1)*dmas
                  enddo
               enddo
            end do
         end do
         do k=1,nz
            do j=1,ny
               yat%vex(i,j,k) = epx(yat%yzn(j),yat%zzn(k),vxcofs)
               yat%vey(i,j,k) = epx(yat%yzn(j),yat%zzn(k),vycofs)
               yat%vez(i,j,k) = epx(yat%yzn(j),yat%zzn(k),vzcofs)
               yat%den(i,j,k) = epx(yat%yzn(j),yat%zzn(k),dencofs)
            end do
         end do
         call calc_quads_new_loop(qxx,qyy,qzz,qxy,qxz,qyz,vxcofs,vycofs,vzcofs,dencofs,modlya%xzn(i),dphidr(i,1,1))
         axx = axx + qxx*dr(i); ayy = ayy + qyy*dr(i); azz = azz + qzz*dr(i);
         axy = axy + qxy*dr(i); axz = axz + qxz*dr(i); ayz = ayz + qyz*dr(i);
      end do
      !$OMP END PARALLEL DO
   end if


   amp2(1) = 2.0d0*gc2vf * geofac * (axx - 1.0d0/3.0d0*(axx+ayy+azz)) / c2
   amp2(2) = 2.0d0*gc2vf * geofac * (ayy - 1.0d0/3.0d0*(axx+ayy+azz)) / c2
   amp2(3) = 2.0d0*gc2vf * geofac * (azz - 1.0d0/3.0d0*(axx+ayy+azz)) / c2
   amp2(4) = gc2vf * geofac * axy / c2
   amp2(5) = gc2vf * geofac * axz / c2
   amp2(6) = gc2vf * geofac * ayz / c2

   !Call the routine used to calculate the signal on our mock data
   amplya = 0.0d0; amplyi = 0.0d0
   if(PRESENT(modlyi)) then
      CALL gw(yat,amplya,1.0d0)
      CALL gw(yit,amplyi,-1.0d0)
      amp1 = (amplya+amplyi)
   else
      CALL gw(yat,amplya,1.0d0)
      amp1 = amplya
   end if
   deALLOCATE (dth, STAT = ast)
   deALLOCATE (dph, STAT = ast)  
   deALLOCATE (sph, STAT = ast)
   deALLOCATE (cph, STAT = ast)
   deALLOCATE (cth, STAT = ast)
   deALLOCATE (sth, STAT = ast)

 end subroutine integral_test

  
  subroutine gw(modl,amp,yf)
    implicit none
    TYPE(fruit), INTENT(IN) :: modl
    REAL(KIND=8), DIMENSION(6), INTENT(OUT) :: amp 
    INTEGER :: nx,ny,nz,i,j,k,ast,ic
    REAL(KIND=8) :: geofac
    REAL(KIND=8),intent(in) :: yf
    parameter(geofac = 1.d0/3.d0)

    REAL(KIND=8) :: pc_cl
    parameter(pc_cl = 2.99792458d+10)

    REAL(KIND=8) :: pc_gc
    parameter(pc_gc = 6.67259d-08)

    REAL(KIND=8) :: pc_pi
    parameter(pc_pi = dacos(-1d0))

    REAL(KIND=8) ::  gc2vf,c2,A_xx,A_yy,A_zz,A_xy,A_yz,A_xz, &
         x,y,z,v_x,v_y,v_z,rsin, &
         dPdx, dPdy,DPdz,dmas,wint,trace
    REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: dphidr(:,:,:),dphidth(:,:,:), & 
         dphidphi(:,:,:)
    REAL(KIND=8), dimension(:),ALLOCATABLE :: dph(:),dth(:),dr(:),cth,sth,cph,sph,axt,grav
    REAL(KIND=8) :: corr_sph,corr_cth,corr_sth,corr_cph
    gc2vf = pc_gc/pc_cl**2d0
    c2 = pc_cl**2d0

    nx=size(modl%den(:,1,1))
    ny=size(modl%den(1,:,1))
    nz=size(modl%den(1,1,:))
    ALLOCATE (dth(ny), STAT = ast)
    ALLOCATE (dph(nz), STAT = ast)
    dth = dcos(modl%yzl(:))-dcos(modl%yzr(:))
    dph = modl%zzr(:)-modl%zzl(:)

    corr_sph = SUM(sin(modl%zzn(:))**2d0*dph(:))
    corr_cph = SUM(cos(modl%zzn(:))**2d0*dph(:))
    
    corr_cth = SUM(cos(modl%yzn(:))**2d0*dth(:))
    corr_sth = SUM(sin(modl%yzn(:))**2d0*dth(:))

    corr_sph = sqrt(pc_pi/corr_sph)
    corr_cth = sqrt(2d0/3d0/corr_cth)
    corr_sth = sqrt(4d0/3d0/corr_sth) 
    corr_cph = sqrt(pc_pi/corr_cph)
    ! corr_cth = 0.99995140318460529d0! dsqrt(2d0/3d0/corr_cth)
    ! corr_sth = 1.000024301064834d0!dsqrt(4d0/3d0/corr_sth)
    ! corr_sph = 0.99999978633434816d0 !Ewalds factor

    ! corr_cth = 0.99995140318460529
    ! corr_sth = 1.000024301064834 

    ! corr_cth = 0.99995140318460529
    ! corr_sth = 1.000024301064834

    ALLOCATE (sph(nz), STAT = ast)
    ALLOCATE (cph(nz), STAT = ast)

    ALLOCATE (cth(ny), STAT = ast)
    ALLOCATE (sth(ny), STAT = ast)

    do j = 1,ny
       cth(j) = dcos(modl%yzn(j))*corr_cth
       sth(j) = dsin(modl%yzn(j))*corr_sth
    end do
    do k =1,nz
       cph(k) = dcos(modl%zzn(k))*corr_cph
       sph(k) = dsin(modl%zzn(k))*corr_sph
    end do

    ALLOCATE (dphidr(nx,ny,nz), STAT = ast)
    ALLOCATE (dphidth(nx,ny,nz), STAT = ast)
    ALLOCATE (dphidphi(nx,ny,nz), STAT = ast)
    dphidth = 0.0d0
    dphidphi = 0.0d0

    do i=1,nx-1
       dphidr(i,:,:) = (modl%gpo(i+1,1:ny,:)-modl%gpo(i,1:ny,:))/(modl%xzr(i)-modl%xzl(i))
    enddo
    dphidr(nx,:,:) = (modl%gpo(nx,1:ny,:)-modl%gpo(nx-1,1:ny,:))/(modl%xzr(nx)-modl%xzl(nx))

    A_xx = 0d0
    A_yy = 0d0
    A_zz = 0d0
    A_xy = 0d0
    A_xz = 0d0
    A_yz = 0d0
    ALLOCATE (dr(nx), STAT = ast)
    ALLOCATE(axt(nx*ny*nz), STAT = ast)
    dr = modl%xzr(:)**3-modl%xzl(:)**3
    do k=1,nz
       do j=1,ny
          do i=100,100!nx

             if(yf .eq. -1.0d0) then
                rsin = modl%xzn(i) * sth(j)
                y = modl%xzn(i)*cth(j)
                x = -rsin * cph(k)
                z = rsin * sph(k)

                v_x = -((sth(j)*modl%vex(i,j,k)+cth(j) * & 
                     modl%vey(i,j,k))*cph(k)  - sph(k)*modl%vez(i,j,k))
                v_z = (sth(j)*modl%vex(i,j,k)+cth(j)* &
                     modl%vey(i,j,k))*sph(k) + cph(k)*modl%vez(i,j,k)
                v_y = cth(j)*modl%vex(i,j,k)-sth(j)* &
                     modl%vey(i,j,k)


                dPdx = -((sth(j)*dphidr(i,j,k) + cth(j)*dphidth(i,j,k) ) * &
                     cph(k) - sph(k)*dphidphi(i,j,k))
                dPdz = (sth(j)*dphidr(i,j,k) + cth(j)*dphidth(i,j,k) ) * &
                     sph(k) + cph(k)*dphidphi(i,j,k)
                dPdy = cth(j)*dphidr(i,j,k) - sth(j)*dphidth(i,j,k)
             else
                rsin = modl%xzn(i) * sth(j)
                z = modl%xzn(i)*cth(j)
                x = rsin * cph(k)
                y = rsin * sph(k)
                
                v_x = ( (sth(j)*modl%vex(i,j,k)+cth(j) * & 
                     modl%vey(i,j,k) )*cph(k)  - sph(k)*modl%vez(i,j,k))
                v_y = (sth(j)*modl%vex(i,j,k)+cth(j)* &
                     modl%vey(i,j,k))*sph(k) + cph(k)*modl%vez(i,j,k)
                v_z = cth(j)*modl%vex(i,j,k)-sth(j)* &
                     modl%vey(i,j,k)


                dPdx =((sth(j)*dphidr(i,j,k) + cth(j)*dphidth(i,j,k) ) * &
                     cph(k) - sph(k)*dphidphi(i,j,k))
                dPdy =(sth(j)*dphidr(i,j,k) + cth(j)*dphidth(i,j,k) ) * &
                     sph(k) + cph(k)*dphidphi(i,j,k)
                dPdz = cth(j)*dphidr(i,j,k) - sth(j)*dphidth(i,j,k)

             end if

             dmas = modl%den(i,j,k)*dr(i)*dth(j)*dph(k)!*modl%yyw(j,k)
             A_xx = A_xx + dmas * ( (v_x*v_x - x*dPdx))
             A_yy = A_yy + dmas * ( (v_y*v_y - y*dPdy))
             A_zz = A_zz + dmas * ( (v_z*v_z - z*dPdz))
             A_xy = A_xy + dmas * (2d0*v_x*v_y - x*dPdy - y*dPdx)
             A_xz = A_xz + dmas * (2d0*v_x*v_z - x*dPdz - z*dPdx)
             A_yz = A_yz + dmas * (2d0*v_y*v_z - y*dPdz - z*dPdy)
          enddo
       enddo
    enddo
    trace = (A_xx + A_yy + A_zz)
    A_xx = 2d0*geofac * (A_xx - 1.0d0/3.0d0 * trace)
    A_yy = 2d0*geofac * (A_yy - 1.0d0/3.0d0 * trace)
    A_zz = 2d0*geofac * (A_zz - 1.0d0/3.0d0 * trace)

    A_xy = geofac * A_xy
    A_xz = geofac * A_xz
    A_yz = geofac * A_yz

    amp(1) = gc2vf * A_xx / c2
    amp(2) = gc2vf * A_yy / c2
    amp(3) = gc2vf * A_zz / c2
    amp(4) = gc2vf * A_xy / c2
    amp(5) = gc2vf * A_xz / c2
    amp(6) = gc2vf * A_yz / c2
    deALLOCATE (dth, STAT = ast)
    deALLOCATE (dph, STAT = ast)  
    deALLOCATE (sph, STAT = ast)
    deALLOCATE (cph, STAT = ast)
    deALLOCATE (cth, STAT = ast)
    deALLOCATE (sth, STAT = ast)
    deALLOCATE (dphidr, STAT = ast)
    deALLOCATE (dphidth, STAT = ast)
    deALLOCATE (dphidphi, STAT = ast)
    deALLOCATE (dr, STAT = ast)
    deALLOCATE(axt, STAT = ast)
  end subroutine gw


  function epx(th,ph,cfs) result(out)
    USE sph_leg_pol
    implicit none
    integer l,m,ast
    real(kind=8),intent(in) :: cfs(:,:)
    real(kind=8) :: cof(0:lmax,-lmax:lmax)
    real(kind=8) :: spa,ph,th,out
    cof = cfs
    out = 0.0d0
    do l = 0,lmax
       do m = -l,l
          out = out + (spharm_r(l,m,th,ph))*cof(l,m)
       end do
    end do
  end function epx

  subroutine trap_test(modl,amp)
    implicit none
    TYPE(fruit), INTENT(IN) :: modl
    REAL(KIND=8), INTENT(OUT) :: amp 
    INTEGER :: nx,ny,nz,i,j,k,ast

    REAL(KIND=8) :: pc_pi
    parameter(pc_pi = dacos(-1d0))

    REAL(KIND=8) ::  A_xx,A_yy,dmas,wint
    REAL(KIND=8), dimension(:),ALLOCATABLE :: dph(:),dth(:),cth,sth,cph,sph
    REAL(KIND=8) :: corr_sph,corr_cth,corr_sth,corr_cph

    nx=size(modl%den(:,1,1))
    ny=size(modl%den(1,:,1))
    nz=size(modl%den(1,1,:))
    !write(*,*) nx,ny,nz
    ALLOCATE (dth(ny), STAT = ast)
    ALLOCATE (dph(nz), STAT = ast)
    dth = dcos(modl%yzl(:))-dcos(modl%yzr(:))
    dph = modl%zzr(:)-modl%zzl(:)

    corr_sph = SUM(sin(modl%zzn(:))**2d0*dph(:))
    corr_cph = SUM(cos(modl%zzn(:))**2d0*dph(:))

    corr_cth = SUM(cos(modl%yzn(:))**2d0*dth(:))
    corr_sth = SUM(sin(modl%yzn(:))**2d0*dth(:))
    
    corr_sph = sqrt(pc_pi/corr_sph)
    corr_cth = sqrt(2d0/3d0/corr_cth)
    corr_sth = sqrt(4d0/3d0/corr_sth) 
    corr_cph = sqrt(pc_pi/corr_cph)

    ALLOCATE (sph(nz), STAT = ast)
    ALLOCATE (cph(nz), STAT = ast)

    ALLOCATE (cth(ny), STAT = ast)
    ALLOCATE (sth(ny), STAT = ast)

    do j = 1,ny
       cth(j) = dcos(modl%yzn(j))*corr_cth
       sth(j) = dsin(modl%yzn(j))*corr_sth
    end do
    do k =1,nz
       cph(k) = dcos(modl%zzn(k))*corr_cph
       sph(k) = dsin(modl%zzn(k))*corr_sph
    end do

    A_xx = 0d0
    A_yy = 0d0
    do k=1,nz
       do j=1,ny
          do i=100,100

             dmas = dth(j)*dph(k)
             A_xx = A_xx + dmas * (3.0d0/512.0d0)*sqrt(pc_pi)*cph(k)**4 * sth(j)**4
             A_yy = A_yy + dmas * (3.0d0/512.0d0)*sqrt(pc_pi)*(cph(k)*sph(k))**2 * sth(j)**4

          enddo
       enddo
    enddo
    write(*,*) A_xx, A_yy
    deALLOCATE (dth, STAT = ast)
    deALLOCATE (dph, STAT = ast)  
    deALLOCATE (sph, STAT = ast)
    deALLOCATE (cph, STAT = ast)
    deALLOCATE (cth, STAT = ast)
    deALLOCATE (sth, STAT = ast)
    amp = 0.0d0
  end subroutine trap_test


  






end module sphexp

