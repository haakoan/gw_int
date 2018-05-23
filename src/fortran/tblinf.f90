module tblinf
  use data_types
  real(kind=8), allocatable, public :: ts1(:),ts2(:),ts3(:),ts4(:),ts5(:),ts6(:),ts7(:),ts8(:),ts9(:)
  real(kind=8), allocatable, public :: ts10(:),ts11(:),ts12(:),ts13(:),ts14(:),ts15(:),ts16(:)
  real(kind=8), allocatable, public :: ts17(:),ts18(:),ts19(:),ts20(:),ts21(:),ts22(:),ts23(:),ts24(:)

  interface
     function lookforkey_c(l1,m1,l2,m2,l3,m3,n1,n2,n3,n4) bind(C, name="lookforkey_")
       use iso_c_binding, only : c_double, c_int
       implicit none
       real(kind=c_double) :: lookforkey_c
       integer :: l1,m1,l2,m2,l3,m3,n1,n2,n3,n4
     end function lookforkey_c

     subroutine set_tbl_c() bind(C, name="set_tbl_")
       !Nothing her but an interface
     end subroutine set_tbl_c
  end interface

contains
  function ILU(l1,m1,l2,m2,l3,m3,n1,n2,n3,n4)
    use iso_c_binding, only : c_double, c_int
    implicit none
    real(kind=c_double) :: ILU
    integer, intent(in) :: l1,m1,l2,m2,l3,m3,n1,n2,n3,n4
    integer:: s1,q1,s2,q2,s3,q3   
    logical :: lps
    integer:: lm1(2,3),lm2(2),i,j
    ILU = 0.0d0
           
    if( MODULO((m1+m2+m3+l1+l2+l3+n1),2) .eq. 0) then
       lm1(1,1) = l1 ; lm1(2,1) = m1 ; lm1(1,2) = l2 ; lm1(2,2) = m2 ; lm1(1,3) = l3 ; lm1(2,3) = m3
       lps = .TRUE.

       do while(lps)
          if(lm1(1,1) .gt. lm1(1,2)) then
             lm2 = lm1(:,1)
             lm1(:,1) = lm1(:,2)
             lm1(:,2) = lm2
             lps = .FALSE.
          end if

          if(lm1(1,2) .gt. lm1(1,3)) then
             lm2 = lm1(:,2)
             lm1(:,2) = lm1(:,3)
             lm1(:,3) = lm2
             lps = .FALSE.
          end if

          if((lm1(2,1) .gt. lm1(2,2)) .and. ( lm1(1,1) .EQ. lm1(1,2))) then
             lm2 = lm1(:,1)
             lm1(:,1) = lm1(:,2)
             lm1(:,2) = lm2
             lps = .FALSE.
          end if

          if((lm1(2,2) .gt. lm1(2,3))  .and. ( lm1(1,2) .EQ. lm1(1,3))) then
             lm2 = lm1(:,2)
             lm1(:,2) = lm1(:,3)
             lm1(:,3) = lm2
             lps = .FALSE.
          end if
          lps = (.NOT. lps)
       end do
       s1 =  lm1(1,1) ; q1= lm1(2,1); s2 = lm1(1,2); q2= lm1(2,2); s3= lm1(1,3); q3= lm1(2,3)

       if((n3 .eq. 0 .and. n4 .eq. 2) .and. ((q1+q2+q3 .eq. -2) .or. (q1+q2+q3 .eq. 2) .or. (q1+q2+q3 .eq. 0) )) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       elseif((n3 .eq. 0 .and. n4 .eq. 0) .and. (q1+q2+q3 .eq. 0)) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       elseif((n3 .eq. 0 .and. n4 .eq. 1) .and. ((q1+q2+q3 .eq. 1) .or. (q1+q2+q3 .eq. -1))) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       elseif((n3 .eq. 2 .and. n4 .eq. 0) .and. ((q1+q2+q3 .eq. -2) .or. (q1+q2+q3 .eq. 2) .or. (q1+q2+q3 .eq. 0) )) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       elseif((n3 .eq. 2 .and. n4 .eq. 1) .and. &
         ((q1+q2+q3 .eq. 1) .or. (q1+q2+q3 .eq. -1) .or. (q1+q2+q3 .eq. 3) .or. (q1+q2+q3 .eq. -3))) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       elseif((n3 .eq. 2 .and. n4 .eq. 2) .and. ((q1+q2+q3 .eq. -4) .or. (q1+q2+q3 .eq. 4) .or. (q1+q2+q3 .eq. 0) )) then
          ILU= lookforkey_c(s1,q1,s2,q2,s3,q3,n1,n2,n3,n4)
       end if
    end if
    
 
    
  end function ILU

  subroutine set_tbl()
    call set_tbl_c()
  end subroutine set_tbl

  function incr(l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,i)
    integer, intent(in) :: l1,m1,l2,m2,l3,m3,n1,n2,n3,n4
    integer:: i,j,incr
    j = i
    if( MODULO((m1+m2+m3+l1+l2+l3+n1),2) .eq. 0) then
       if((n3 .eq. 0 .and. n4 .eq. 2) .and. ((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2) .or. (m1+m2+m3 .eq. 0) )) then
          j=i+1
       elseif((n3 .eq. 0 .and. n4 .eq. 1) .and. ((m1+m2+m3 .eq. 1) .or. (m1+m2+m3 .eq. -1))) then
          j=i+1
       elseif((n3 .eq. 0 .and. n4 .eq. 0) .and. (m1+m2+m3 .eq. 0)) then
          j=i+1	
       elseif((n3 .eq. 2 .and. n4 .eq. 0) .and. ((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2) .or. (m1+m2+m3 .eq. 0) )) then
          j=i+1
       elseif((n3 .eq. 2 .and. n4 .eq. 1) .and. &
         ((m1+m2+m3 .eq. 1) .or. (m1+m2+m3 .eq. -1) .or. (m1+m2+m3 .eq. 3) .or. (m1+m2+m3 .eq. -3))) then
          j=i+1
       elseif((n3 .eq. 2 .and. n4 .eq. 2) .and. ((m1+m2+m3 .eq. -4) .or. (m1+m2+m3 .eq. 4)  .or. (m1+m2+m3 .eq. 0))) then
          j=i+1
       end if
    end if
    incr = j
  end function incr

  subroutine set_arrys()
    implicit none
    real(kind=8), allocatable :: ys1(:)
    integer :: i4
    integer :: ast,arl
    arl = ((lmax+1)*(2*lmax+1))**3

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 2, 0, 2)
    ALLOCATE(ts1(i4),STAT=ast)
    ts1=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 0, 2, 0)
    ALLOCATE(ts2(i4),STAT=ast)
    ts2=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 1, 0, 2)
    ALLOCATE(ts3(i4),STAT=ast)
    ts3=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 2, 0, 0, 2)
    ALLOCATE(ts4(i4),STAT=ast)
    ts4=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 2, 2, 0)
    ALLOCATE(ts5(i4),STAT=ast)
    ts5=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 0, 0, 2)
    ALLOCATE(ts6(i4),STAT=ast)
    ts6=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 1, 2, 0)
    ALLOCATE(ts7(i4),STAT=ast)
    ts7=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 2, 0, 2, 0)
    ALLOCATE(ts8(i4),STAT=ast)
    ts8=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 2, 0, 0, 0)
    ALLOCATE(ts9(i4),STAT=ast)
    ts9=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 1, 0, 0)
    ALLOCATE(ts10(i4),STAT=ast)
    ts10=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 2, 0, 0)
    ALLOCATE(ts11(i4),STAT=ast)
    ts11=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 1, 0, 2)
    ALLOCATE(ts12(i4),STAT=ast)
    ts12=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 1, 2, 0)
    ALLOCATE(ts13(i4),STAT=ast)
    ts13=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 0, 0, 2)
    ALLOCATE(ts14(i4),STAT=ast)
    ts14=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 0, 2, 0)
    ALLOCATE(ts15(i4),STAT=ast)
    ts15=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 1, 0, 1)
    ALLOCATE(ts16(i4),STAT=ast)
    ts16=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 2, 0, 0, 1)
    ALLOCATE(ts17(i4),STAT=ast)
    ts17=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 2, 0, 1)
    ALLOCATE(ts18(i4),STAT=ast)
    ts18=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 1, 0, 0, 1)
    ALLOCATE(ts19(i4),STAT=ast)
    ts19=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array(ys1,i4, 0, 1, 0, 1)
    ALLOCATE(ts20(i4),STAT=ast)
    ts20=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)


    !Outside of the first loop!
    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array1(ys1,i4, 0, 2, 0, 2)
    ALLOCATE(ts21(i4),STAT=ast)
    ts21=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array1(ys1,i4, 0, 2, 2, 0)
    ALLOCATE(ts22(i4),STAT=ast)
    ts22=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array1(ys1,i4, 2, 0, 0, 0)
    ALLOCATE(ts23(i4),STAT=ast)
    ts23=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)

    i4=0
    ALLOCATE(ys1(arl),STAT=ast)
    call set_array1(ys1,i4, 1, 1, 0, 1)
    ALLOCATE(ts24(i4),STAT=ast)
    ts24=ys1(1:i4)
    DEALLOCATE(ys1,STAT=ast)
    

  end subroutine set_arrys


  
  subroutine set_array(arry,i3,n1,n2,n3,n4)
    implicit none
    integer,intent(in) :: n1,n2,n3,n4
    real(kind=8),intent(inout) :: arry(:)
    integer :: arl,ast,l_idx,m_idx,lim
    integer :: l1,m1,l2,m2,l3,m3,i1,i2
    integer, intent(out) :: i3
    i1 = 0 ; i2 = 0
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
          i2 = incr(l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,i1)
          if(i2 .gt. i1) then
             i1 = i2
             arry(i1) = ILU(l1, m1, l2, m2, l3, m3, n1, n2, n3, n4)
          end if
       end do
    end do
    i3 = i1
  end subroutine set_array

  subroutine set_array1(arry,i3,n1,n2,n3,n4)
    implicit none
    integer,intent(in) :: n1,n2,n3,n4
    real(kind=8),intent(inout) :: arry(:)
    integer :: arl,ast
    integer :: l1,m1,i1,i2
    integer, intent(out) :: i3
    i1 = 0 ; i2 = 0
    do l1 = 0,lmax
       do m1 = -l1,l1
          i2 = incr(l1,m1,0,0,0,0,n1,n2,n3,n4,i1)
          if(i2 .gt. i1) then
             i1 = i2
             arry(i1) = ILU(l1, m1, 0, 0, 0, 0, n1, n2, n3, n4)
          end if
         end do
    end do
    i3 = i1
  end subroutine set_array1
 
  
  
end module tblinf

