program gen_table
  use ffastwigxj
  use fwigxjpf
  use al_int
 

  !call make_table()
  call test_sphint()


contains

  subroutine test_sphint()
    implicit none
    real(kind=8) :: agr,ft,res
    integer :: l1,l2,l3,m1,m2,m3,n1,n2,n3,n4,i
    real(kind=8) :: pc_pi
    parameter(pc_pi = acos(-1.d0))
    call ffastwigxj_load("table_50.3j",3);

    
    ! res = ISPH(0,0,0,0,1,-1,0,0,0,1)
    ! ft = 0.135675235229675186130d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,0,0,1,-1,0,0,0,1
    ! write(*,'(1E20.10)') res,ft,res/ft

    ! res = ISPH(0,0,0,0,0,0,0,2,0,0)
    ! ft = 0.188063194515919d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,0,0,0,0,0,2,0,0
    ! write(*,'(1E20.10)') res,ft,res/ft



    ! res = ISPH(1,1,0,0,1,1,0,0,0,2)
    ! ft = 0.07052369794346953d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 1,1,0,0,1,1,0,0,0,2
    ! write(*,'(1E20.10)') res,ft,res/ft
   
    ! res = ISPH(1,1,0,0,1,1,0,0,2,0)
    ! ft = -0.07052369794346953d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 1,1,0,0,1,1,0,0,2,0
    ! write(*,'(1E20.10)') res,ft,res/ft
   
    ! res = ISPH(1,1,1,0,1,0,0,0,2,1)
    ! ft = -0.02543910661d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 1,1,1,0,1,0,0,0,2,1
    ! write(*,'(1E20.10)') res,ft,res/ft
   
    ! res = ISPH(2,2,2,2,2,0,0,0,2,2)
    ! ft = 0.01126398447d0
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 1,1,1,0,1,0,0,0,2,1
    ! write(*,'(1E20.10)') res,ft,res/ft



    ! res = ISPH(0,0,1,1,1,1,0,2,0,2)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,1,0,2,0,2
    ! write(*,'(1E20.10)') res

    ! res = ISPH(0,0,1,1,1,-1,0,2,0,2)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,-1,0,2,0,2
    ! write(*,'(1E20.10)') res

    ! res = ISPH(0,0,1,1,1,1,2,2,0,0)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,1,2,2,0,0
    ! write(*,'(1E20.10)') res

    ! res = ISPH(0,0,1,1,1,1,2,0,0,2)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,1,2,0,0,2
    ! write(*,'(1E20.10)') res


    ! res = ISPH(0,0,1,1,1,1,0,2,2,0)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,1,0,2,2,0
    ! write(*,'(1E20.10)') res

    ! res = ISPH(1,1,1,-1,1,1,1,0,0,1)
    ! write(*,'(A)') "l1,m1,l2,m2,l3,m3,n1,n2,n3,n4"
    ! write(*,'(10(I1))') 0,0,1,1,1,1,0,0,2,2
    ! write(*,'(1E20.10)') res


    res = I3(1,-1,1,-1,0,0,0,2)
    write(*,'(1E20.10)') res
    res = ISPH(1,-1,1,-1,0,0,2,0,0,2)
    write(*,'(1E20.10)') res

    !Speed test


    ! do i = 1,10000000
    !    res = ISPH(1,1,1,0,1,0,0,0,2,1)
    !    res = ISPH(5,2,10,-9,2,2,0,0,2,2)
    ! end do
    call ffastwigxj_unload(3);
    call fwig_temp_free();
    call fwig_table_free();
  end subroutine test_sphint


  subroutine make_table()
    implicit none
    real(kind=8) :: agr,ft,res
    integer :: l1,l2,l3,m1,m2,m3,n1,n2,n3,n4
    integer :: er,r(6),i
    integer :: nlines,ctr
    integer, allocatable :: cbs(:,:)
    real(kind=8) :: pc_pi
    parameter(pc_pi = acos(-1.d0))

2000 format(6(I2))
2001 format((I4,I4,I4,I4,I4,I4,I4,I4,I4,I4,E20.10))
    open(352, file="../../unik_combs.dat",action='read')
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
    open(353, file="../../table1v2.dat",action='write')
    call ffastwigxj_load("table_100.3j",3);
    call fwig_table_init(2*200,9)
    call fwig_temp_init(2*200)
    do i = 1, nlines
       l1 = cbs(i,1) ; m1 = cbs(i,2)
       l2 = cbs(i,3) ; m2 = cbs(i,4)
       l3 = cbs(i,5) ; m3 = cbs(i,6)  
       do n1 = 0,2
          do n2 = 0,2
             do n3 = 0,2
                do n4 = 0,2
                   res = ISPH(l1,m1,l2,m2,l3,m3,n1,n2,n3,n4)
                   if(abs(res) .ne. 0.0d0) then
                      write(*,*) l1,m1,l2,m2,l3,m3,n1,n2,n3,n4
                      write(353,2001),l1,m1,l2,m2,l3,m3,n1,n2,n3,n4,res
                   end if
                end do
             end do
          end do
       end do
    end do
    call ffastwigxj_unload(3);
    call ffastwigxj_unload(3);
    call fwig_temp_free();
    call fwig_table_free();
    close(353)
    close(352)
  end subroutine make_table

end program gen_table

  ! write(*,*) I3(2,-1,2,1,2,2,0,2)
  ! write(*,*) I3(2,-1,2,1,2,2,0,1)
  ! write(*,*) I3(2,-2,2,1,2,2,1,1)
  ! write(*,*) I3(2,-2,2,1,2,2,1,0)
  ! write(*,*) I3(2,-1,2,1,2,2,2,0)
  
  ! write(*,*) I3(2,1,4,3,6,4,0,0)
  ! write(*,*) I1(6,4)
  ! !write(*,*) G(2,1,2,2,3,3),G(2,1,2,2,4,3)
  
