program makecode
  use fint
  implicit none
  real(kind=8) :: agr,ft,res
  integer :: l1,l2,l3,m1,m2,m3,n1,n2,n3,n4
  integer :: er,r(6),i
  integer :: nlines,ctr
  integer, allocatable :: cbs(:,:)
  real(kind=8) :: pc_pi
  parameter(pc_pi = acos(-1.d0))
  !call cint(1,2,0,1,2,3,4,-1,0,-3,agr)
  !write(*,*) agr
  !call cint(1,2,0,1,2,3,4,-1,1,-3,agr)
  !write(*,*) agr
2000 format(6(I2))
2001 format((I2,I2,I2,I2,I2,I2,I2,I2,I2,I2,1E20.10))
  open(352, file="../unik_combs.dat",action='read')
  nlines = 0
  do
     read(352,'(a)',iostat=er) res
     if (er/=0) EXIT
     nlines = nlines+1
  end do
  rewind(352)
  write(*,*) nlines
  ALLOCATE (cbs(nlines,6), STAT = er)
  do i = 1,nlines
     read(352,*) cbs(i,:)
  enddo
  ctr = 0
  write(*,*) "We start with a total of ", nlines
  open(353, file="../results.dat",action='write')
  do i = 1, nlines
!     read(352,*) l1,m1,l2,m2,l3,m3
     l1 = cbs(i,1) ; m1 = cbs(i,2)
     l2 = cbs(i,3) ; m2 = cbs(i,4)
     l3 = cbs(i,5) ; m3 = cbs(i,6)
     do n1 = 0,2
        if(mod(l1+l2+l3+m1+m2+m3+n1,2) .eq. 0) then
           do n2 = 0,2
              do n3 = 0,2
                 do n4 = 0,2
                    select case(n3)
                    case(0)
                       select case(n4)
                       case(0)
                          res = 0d0
                       case(1)
                          if((m1+m2+m3 .eq. -1) .or. (m1+m2+m3 .eq. 1)) then
                             ft = pc_pi
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                             !write(*,*) "meee",res, m1+m2+m3
                             !write(*,*) l1,m1,l2,m2,l3,m3,n1,n2
                          end if
                       case(2)
                          if((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2)) then
                             ft = pc_pi*0.5d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                          if((m1+m2+m3 .eq. 0)) then
                             ft = pc_pi
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                       end select
                    case(1)
                       res = 0d0
                    case(2)
                       select case(n4)
                       case(0)
                          if((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2)) then
                             ft = pc_pi*0.5d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                          if((m1+m2+m3 .eq. 0)) then
                             ft = pc_pi
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                       case(1)
                          if((m1+m2+m3 .eq. -3) .or. (m1+m2+m3 .eq. 3)) then
                             ft = -pc_pi*0.25d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                          if((m1+m2+m3 .eq. -1) .or. (m1+m2+m3 .eq. 1)) then
                             ft = pc_pi*0.25d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                       case(2)
                          if((m1+m2+m3 .eq. -4) .or. (m1+m2+m3 .eq. 4)) then
                             ft = -pc_pi*0.125d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                          if((m1+m2+m3 .eq. 0)) then
                             ft = pc_pi*0.125d0
                             call cint(l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,agr)
                             res = ft*agr
                             write(353,2001) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                             ctr = ctr+1
                          end if
                       end select
                    end select
                 end do
              end do
           end do
        end if
     end do
  end do
  write(*,*) "At the end we have", ctr
  write(*,*) res
  close(352)
  close(353)

end program makecode
