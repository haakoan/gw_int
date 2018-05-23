program main
  use sphexp
  use tblinf
  CHARACTER(LEN=90) :: filename,filein
  integer :: clck_counts_beg, clck_counts_end, clck_rate
  real ::  beg_cpu_time, end_cpu_time

  CALL getarg(1, filein)
  filename = 'sph_an'//trim(filein)//'.dat' 	
  open(332,FILE=trim(filename),action="write")
  filename = 'gw_rot'//trim(filein)//'.dat'
  open(333,FILE=trim(filename),action="write")
  close(333)
  close(332)


  call set_tbl()
  call set_arrys()
  write(*,*) "Read table and arrays are set"

  call system_clock ( clck_counts_beg, clck_rate )
  call cpu_time (beg_cpu_time)
  call run_int_test()
  call cpu_time (end_cpu_time)
  call system_clock ( clck_counts_end, clck_rate )
  write (*, *) "V2 Wall:",  (clck_counts_end - clck_counts_beg) / real (clck_rate)
  write (*, *) "V2 CPU", end_cpu_time - beg_cpu_time

  


end program main
