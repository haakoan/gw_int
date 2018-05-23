module quads
  use data_types
  use tblinf
  implicit none
  real(kind=8),private :: pc_pi, trd,frd,erd,it1,it2
  !real(kind=8), allocatable, dimension(:,:),private :: vr,vph,vth,rho
  !real(kind=8),private :: r,pot
  parameter(pc_pi = acos(-1.d0))
  parameter(trd = 2.d0/3.d0)
  parameter(frd = 4.d0/3.d0)
  parameter(erd = 8.d0/3.d0)
  parameter(it1 = sqrt(4.0d0*pc_pi))
  parameter(it2 = (4.0d0*pc_pi))
  
contains

  ! subroutine set_vars(vel1,vel2,vel3,den,radi,gpo)
  !   real(kind=8),dimension(:,:), intent(in) :: vel1,vel2,vel3,den
  !   real(kind=8), intent(in) :: radi,gpo
  !   integer :: i,ast
  !   ALLOCATE (vr(0:lmax,-lmax:lmax), STAT = ast)
  !   ALLOCATE (vph(0:lmax,-lmax:lmax), STAT = ast)
  !   ALLOCATE (vth(0:lmax,-lmax:lmax), STAT = ast)
  !   ALLOCATE (rho(0:lmax,-lmax:lmax), STAT = ast)

  !   vr = vel1 ; vph = vel3; vth = vel2; rho = den;
  !   r = radi; pot = gpo;
  ! end subroutine set_vars

    

  subroutine calc_quads_new_loop(Q11,Q22,Q33,Q12,Q13,Q23,vel1,vel2,vel3,den,radi,gpo)
    implicit none
    real(kind=8),dimension(:,:), intent(in) :: vel1,vel2,vel3,den
    real(kind=8), intent(in) :: radi,gpo
    integer :: i,j,ast,lim,m_idx,l_idx
   
    real(kind=8), intent(out) :: Q11,Q22,Q33,Q12,Q13,Q23

    real(kind=8) :: res,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
    real(kind=8) :: tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16
    real(kind=8) :: tmp17,tmp18,tmp19,tmp20,tmp21,tmp22,tmp23,tmp24

    real(kind=8) :: qt11,qt22,qt33,qt12,qt13,qt23
    real(kind=8), allocatable, dimension(:,:) :: vr,vph,vth,rho
    real(kind=8) :: r,pot, pre_factor
    REAL(KIND=8) :: pc_cl
    parameter(pc_cl = 2.99792458d+10)
    REAL(KIND=8) :: pc_gc
    parameter(pc_gc = 6.67259d-08)

    integer :: l1,m1,l2,m2,l3,m3
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14
    integer :: i15,i16,i17,i18,i19,i20,i21,i22,i23,i24
    integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14
    integer :: j15,j16,j17,j18,j19,j20,j21,j22,j23,j24
   
    ALLOCATE (vr(0:lmax,-lmax:lmax), STAT = ast)
    ALLOCATE (vph(0:lmax,-lmax:lmax), STAT = ast)
    ALLOCATE (vth(0:lmax,-lmax:lmax), STAT = ast)
    ALLOCATE (rho(0:lmax,-lmax:lmax), STAT = ast)
    vr = vel1 ; vph = vel3; vth = vel2;
    rho = den;
    r = radi; pot = gpo;

    pre_factor =1.0d0! pc_gc/pc_cl**4

    i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0; i6 = 0; i7 = 0; i8 = 0; i9 = 0;
    i10 = 0; i11 = 0; i12 = 0; i13 = 0; i14 = 0; i15 = 0; i16 = 0;
    i17 = 0; i18 = 0; i19 = 0; i20 = 0;
    j1 = 0; j2 = 0; j3 = 0; j4 = 0; j5 = 0; j6 = 0; j7 = 0; j8 = 0; j9 = 0;
    j10 = 0; j11 = 0; j12 = 0; j13 = 0; j14 = 0; j15 = 0; j16 = 0;
    j17 = 0; j18 = 0; j19 = 0; j20 = 0;
    Q11 = 0.0d0; Q22 = 0.0d0; Q33 = 0.0d0; Q12 = 0.0d0; Q13 = 0.0d0; Q23 = 0.0d0;
    qt11 = 0.0d0; qt22 = 0.0d0; qt33 = 0.0d0; qt12 = 0.0d0; qt13 = 0.0d0; qt23 = 0.0d0;
    lim = lmax+1

    do l_idx = 0,lim**3-1
       l1 = l_idx/(lim*lim)
       l2 = MOD(l_idx,lim*lim)/lim
       l3 = MOD(MOD(l_idx,lim*lim),lim)
       ! implement quick return
       do m_idx = 0,(2*l1+1)*(2*l2+1)*(2*l3+1)-1
          m1 = m_idx / ((2*l2+1)*(2*l3+1)) - l1
          m2 = MOD(m_idx,(2*l2+1)*(2*l3+1))/(2*l3+1) - l2
          m3 = MOD(MOD(m_idx,(2*l2+1)*(2*l3+1)),(2*l3+1)) - l3
          ! implement quick return

          tmp1 = 0.0d0; tmp2 = 0.0d0; tmp3 = 0.0d0;
          tmp4 = 0.0d0; tmp5 = 0.0d0; tmp6 = 0.0d0;
          tmp7 = 0.0d0; tmp8 = 0.0d0; tmp9 = 0.0d0;
          tmp10 = 0.0d0; tmp11 = 0.0d0; tmp12 = 0.0d0;
          tmp13 = 0.0d0; tmp14 = 0.0d0; tmp15 = 0.0d0;
          tmp16 = 0.0d0; tmp17 = 0.0d0; tmp18 = 0.0d0;
          tmp19 = 0.0d0; tmp20 = 0.0d0;
          
          j1 = incr(l1, m1, l2, m2, l3, m3, 0, 2, 0, 2,i1)
          if(j1 .gt.i1) then
             i1=j1
             tmp1 = ts1(i1)!ILU(l1, m1, l2, m2, l3, m3, 2, 0, 0, 2)
          endif

          j2 = incr(l1, m1, l2, m2, l3, m3, 0, 0, 2, 0,i2)
          if(j2 .gt.i2) then
             i2=j2
             tmp2 = ts2(i2)!ILU(l1, m1, l2, m2, l3, m3, 0, 0, 2, 0)
          endif

          j3 = incr(l1, m1, l2, m2, l3, m3, 1, 1, 0, 2,i3)
          if(j3 .gt.i3) then
             i3=j3
             tmp3 = ts3(i3)!ILU(l1, m1, l2, m2, l3, m3, 1, 1, 0, 2)
          endif

          j4 = incr(l1, m1, l2, m2, l3, m3, 2, 0, 0, 2,i4)
          if(j4 .gt.i4) then
             i4=j4
             tmp4 = ts4(i4)!ILU(l1, m1, l2, m2, l3, m3, 0, 2, 0, 2)
          endif

          j5 = incr(l1, m1, l2, m2, l3, m3, 0, 2, 2, 0,i5)
          if(j5 .gt.i5) then
             i5=j5
             tmp5 = ts5(i5)!ILU(l1, m1, l2, m2, l3, m3, 2, 0, 2, 0)
          endif

          j6 = incr(l1, m1, l2, m2, l3, m3, 0, 0, 0, 2,i6)
          if(j6 .gt.i6) then
             i6=j6
             tmp6 = ts6(i6)!ILU(l1, m1, l2, m2, l3, m3, 0, 0, 0, 2)
          endif

          j7 = incr(l1, m1, l2, m2, l3, m3, 1, 1, 2, 0,i7)
          if(j7 .gt.i7) then
             i7=j7
             tmp7 = ts7(i7)!ILU(l1, m1, l2, m2, l3, m3, 1, 1, 2, 0)
          endif

          j8 = incr(l1, m1, l2, m2, l3, m3, 2, 0, 2, 0,i8)
          if(j8 .gt.i8) then
             i8=j8
             tmp8 = ts8(i8)!ILU(l1, m1, l2, m2, l3, m3, 0, 2, 2, 0)
          endif

          j9 = incr(l1, m1, l2, m2, l3, m3, 2, 0, 0, 0,i9)
          if(j9 .gt.i9) then
             i9=j9
             tmp9 = ts9(i9)!ILU(l1, m1, l2, m2, l3, m3, 0, 2, 0, 0)
          endif

          j10 = incr(l1, m1, l2, m2, l3, m3, 1, 1, 0, 0,i10)
          if(j10 .gt.i10) then
             i10=j10
             tmp10 = ts10(i10)!ILU(l1, m1, l2, m2, l3, m3, 1, 1, 0, 0)
          endif

          j11 = incr(l1, m1, l2, m2, l3, m3, 0, 2, 0, 0,i11)
          if(j11 .gt.i11) then
             i11=j11
             tmp11 = ts11(i11)!ILU(l1, m1, l2, m2, l3, m3, 2, 0, 0, 0)
          endif

          j12 = incr(l1, m1, l2, m2, l3, m3, 0, 1, 0, 2,i12)
          if(j12 .gt.i12) then
             i12=j12
             tmp12 = ts12(i12)!ILU(l1, m1, l2, m2, l3, m3, 1, 0, 0, 2)
          endif

          j13 = incr(l1, m1, l2, m2, l3, m3, 0, 1, 2, 0,i13)
          if(j13 .gt.i13) then
             i13=j13
             tmp13 = ts13(i13)!ILU(l1, m1, l2, m2, l3, m3, 1, 0, 2, 0)
             if(tmp13 .eq. 0.0d0) then
                write(*,*) i13,l1,m1,l2,m2,l3,m3
             end if
          endif

          j14 = incr(l1, m1, l2, m2, l3, m3, 1, 0, 0, 2,i14)
          if(j14 .gt.i14) then
             i14=j14
             tmp14 = ts14(i14)!ILU(l1, m1, l2, m2, l3, m3, 0, 1, 0, 2)
          endif

          j15 = incr(l1, m1, l2, m2, l3, m3, 1, 0, 2, 0,i15)
          if(j15 .gt.i15) then
             i15=j15
             tmp15 = ts15(i15)!ILU(l1, m1, l2, m2, l3, m3, 0, 1, 2, 0)
          endif

          j16 = incr(l1, m1, l2, m2, l3, m3, 1, 1, 0, 1,i16)
          if(j16 .gt.i16) then
             i16=j16
             tmp16 = ts16(i16)!ILU(l1, m1, l2, m2, l3, m3, 1, 1, 0, 1)
          endif

          j17 = incr(l1, m1, l2, m2, l3, m3, 2, 0, 0, 1,i17)
          if(j17 .gt.i17) then
             i17=j17
             tmp17 = ts17(i17)!ILU(l1, m1, l2, m2, l3, m3, 0, 2, 0, 1)
             if(tmp17 .eq. 0.0d0) then
                write(*,*) i17,l1,m1,l2,m2,l3,m3
             end if
         
          endif

          j18 = incr(l1, m1, l2, m2, l3, m3, 0, 2, 0, 1,i18)
          if(j18 .gt.i18) then
             i18=j18
             tmp18 = ts18(i18)!ILU(l1, m1, l2, m2, l3, m3, 2, 0, 0, 1)
          endif

          j19 = incr(l1, m1, l2, m2, l3, m3, 1, 0, 0, 1,i19)
          if(j19 .gt.i19) then
             i19=j19
             tmp19 = ts19(i19)!ILU(l1, m1, l2, m2, l3, m3, 0, 1, 0, 1)
          endif
          
          j20 = incr(l1, m1, l2, m2, l3, m3, 0, 1, 0, 1,i20)
          if(j20 .gt.i20) then
             i20=j20
             tmp20 = ts20(i20)!ILU(l1, m1, l2, m2, l3, m3, 1, 0, 0, 1)
          endif

          qt11 = qt11 + &   
               !pot*rho*tmp1*r +&
               (vph(l2,m2)*vph(l3,m3)*tmp2 +&
               vr(l2,m2)*vr(l3,m3)*tmp1 +&
               2.0d0*vr(l2,m2)*vth(l3,m3)*tmp3 +&
               vth(l2,m2)*vth(l3,m3)*tmp4)*rho(l1,m1)
          qt22 = qt22 +&
               !pot*rho(l1,m1)*tmp5*r +&
               (vph(l2,m2)*vph(l3,m3)*tmp6 +&
               vr(l2,m2)*vr(l3,m3)*tmp5 +&
               2.0d0*vr(l2,m2)*vth(l3,m3)*tmp7 +&
               vth(l2,m2)*vth(l3,m3)*tmp8)*rho(l1,m1)

          qt33 = qt33 +&
               !pot*rho(l1,m1)*tmp9*r +&
               (vr(l2,m2)*vr(l3,m3)*tmp9 -&
               2.0d0*vr(l2,m2)*vth(l3,m3)*tmp10 +&
               vth(l2,m2)*vth(l3,m3)*tmp11)*rho(l1,m1)

          qt12 = qt12 + &               
               (vph(l2,m2)*vr(l3,m3)*tmp12 -&
               vph(l2,m2)*vr(l3,m3)*tmp13 +&
               vph(l2,m2)*vth(l3,m3)*tmp14 -&
               vph(l2,m2)*vth(l3,m3)*tmp15)*2.0d0*rho(l1,m1)
          !write(*,*) l1,m1,l2,m2,l3,m3
          !write(*,*) tmp12,tmp13,tmp14,tmp15
          !write(*,*) "----------------------------------"

               
          qt13 = qt13 +&
               !2*pot*rho(l1,m1)*tmp16*r +&
               (vr(l2,m2)*vr(l3,m3)*tmp16 +&
               vr(l2,m2)*vth(l3,m3)*tmp17 -&
               vr(l2,m2)*vth(l3,m3)*tmp18 -&
               vth(l2,m2)*vth(l3,m3)*tmp16)*2.0d0*rho(l1,m1)


          qt23 = qt23 + &
               (vph(l2,m2)*vr(l3,m3)*tmp19 -&
               vph(l2,m2)*vth(l3,m3)*tmp20)*2.0d0*rho(l1,m1)
       end do
    end do
    i21 = 0; i22 = 0; i23 = 0; i24 = 0;  
    j21 = 0; j22 = 0; j23 = 0; j24 = 0;

    do l1 = 0,lmax       
       do m1 = -l1,l1
          tmp21 = 0.0d0; tmp22 = 0.0d0; tmp23 = 0.0d0;
          tmp24 = 0.0d0;
          !pot*rho*tmp1*r +&
          !pot*rho(l1,m1)*tmp5*r +&
          !pot*rho(l1,m1)*tmp9*r +&
          !2*pot*rho(l1,m1)*tmp16*r +&

          j21 = incr(l1, m1, 0, 0, 0, 0, 0, 2, 0, 2,i21)
          if(j21 .gt.i21) then
             i21=j21
             tmp21 = ts21(i21)!ILU(l1, m1, 0, 0, 0, 0, 2, 0, 0, 2)
          endif

          j22 = incr(l1, m1, 0, 0, 0, 0, 0, 2, 2, 0,i22)
          if(j22 .gt.i22) then
             i22=j22
             tmp22 = ts22(i22)!ILU(l1, m1, 0, 0, 0, 0, 2, 0, 2, 0)
          endif

          j23 = incr(l1, m1, 0, 0, 0, 0, 2, 0, 0, 0,i23)
          if(j23 .gt.i23) then
             i23=j23
             tmp23 = ts23(i23)!ILU(l1, m1, 0, 0, 0, 0, 0, 2, 0, 0)
          endif

          j24 = incr(l1, m1, 0, 0, 0, 0, 1, 1, 0, 1,i24)
          if(j24 .gt.i24) then
             i24=j24
             tmp24 = ts24(i24)!ILU(l1, m1, 0, 0, 0, 0, 1, 1, 0, 1)
          endif

          qt11 = qt11 -pot*tmp21*r*it2*rho(l1,m1)
          qt22 = qt22 -pot*tmp22*r*it2*rho(l1,m1)
          qt33 = qt33 -pot*tmp23*r*it2*rho(l1,m1)
          qt13 = qt13 -2.0d0*pot*tmp24*r*it2*rho(l1,m1)


       end do
    end do

    Q11 = Q11+ pre_factor*qt11 ; Q22 = Q22+ pre_factor*qt22 ; Q33 = Q33+ pre_factor*qt33
    Q12 = Q12+ pre_factor*qt12 ; Q13 = Q13+ pre_factor*qt13 ; Q23 = Q23+ pre_factor*qt23

  end subroutine calc_quads_new_loop

end module quads

