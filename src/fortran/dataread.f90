module dataread
 USE HDF5
 USE data_types
contains
  subroutine readdata(filename,hdl,s,su,data,yy)
    !use iso_c_binding
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: yy
    CHARACTER(LEN=90),INTENT(IN) :: filename
    CHARACTER(len=:),allocatable,INTENT(IN) :: hdl
    CHARACTER(LEN=20) :: gname
    INTEGER(HID_T) :: file,mgr,did,ugr,tp,mgryy ! Handle
    INTEGER(HID_T) :: sid,redo
    INTEGER(hsize_t),INTENT(IN) :: s
    INTEGER(hsize_t),INTENT(OUT) :: su
    INTEGER :: status,ndims,alstatus
    INTEGER(hsize_t) :: idx,sts
    INTEGER :: ret_value
    INTEGER :: stype,nlinks,maxc
    LOGICAL :: isset
    REAL(KIND=8),INTENT(OUT), dimension(:,:,:,:),ALLOCATABLE :: data(:,:,:,:)
    INTEGER(hsize_t), dimension(:),ALLOCATABLE :: dms(:),dmsm(:)
    CALL H5open_f(status)
    !Open file and find grps
    !Here we first need to open the main file. Then we need to 
    !search through the file and find all time-step groups,
    !the name of these groups and the number of them is not the
    !same in every file.
    CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file, status)!Open file
    CALL H5Gopen_f(file,"/", mgr, status)!Open outer grp
    CALL H5Gget_info_f(mgr, stype, nlinks, maxc, status, isset)!Read info about outer grp
    !This line is key. It finds the name of all the grps in the file.
    CALL H5Lget_name_by_idx_f(file, "/", H5_INDEX_NAME_F,H5_ITER_INC_F, &
         s,gname,status,sts)!,H5F_ACC_RDONLY_F)  
    !Open steps and group. Several steps per file 
    su = nlinks
    !If the YinYang grid is used we need to enter subgrps in the main grp.
    !The file will look something like this, in the case of yinyang, 
    !/Step#001
    !|-->/Yang
    !    |-->den
    !    |-->gpo
    !    |-->...
    !|-->/Yin
    !    |-->den
    !    |-->gpo
    !    |-->...
    !|-->time
    !|-->yinyang_weight
    !|-->...
    !/Step#002
    !/...
    ! and something like this in the case of a normal grid
    !/Step#001
    !|-->den
    !|-->gpo
    !|-->time
    !/Step#002
    !/...
    !The way the code is set up we first open the yin part, run our code on that
    !then open the yang part and run the same code on that.
    !This might suboptimal with regards to speed. 
    if(yy .EQ. 1) then
       CALL H5Gopen_f(file, gname, mgryy, status)
       CALL H5Gopen_f(mgryy,"Yin", ugr, status) 
    elseif(yy .EQ. 2) then
       CALL H5Gopen_f(file, gname, mgryy, status)
       CALL H5Gopen_f(mgryy,"Yang", ugr, status)
       !If we use no yinyang grid just open in the "normal way"
    else
       CALL H5Gopen_f(file, gname, ugr, status)
    end if

    !if(hdl .eq. "ish") then
    !   tp = H5T_NATIVE_INTEGER
    !else
    tp = H5T_NATIVE_DOUBLE
    !tp = H5T_IEEE_F64LE
    redo = H5F_ACC_RDONLY_F
    !end if
    !yingyang_weight is used as a weight when integrating over the grid, if we use the normal
    !grid then this is equal to one for every cell.
    !Since the yinyang_weight dataset lives above the yin and yang grps we need to take care
    !when opening it. This is what I do here. The same is true for the time array.
    if(((hdl .eq. "yinyang_weight") .or. (hdl .eq. "time") &
         .or. (hdl .eq. "xzr") .or. (hdl .eq. "xzn") .or. (hdl .eq. "xzl") &
         .or. (hdl .eq. "yzr") .or. (hdl .eq. "yzn") .or. (hdl .eq. "yzl") &
         .or. (hdl .eq. "zzr") .or. (hdl .eq. "zzn") .or. (hdl .eq. "zzl"))  &
         .AND. ((yy .eq. 1) .or. (yy .eq.2))) then
       CALL H5Dopen_f(mgryy, hdl, did, status, redo) 
       !    CALL H5Dget_type_f(did,tp,status)
    else
       CALL H5Dopen_f(ugr, hdl, did, status, redo) 
       !     CALL H5Dget_type_f(did,tp,status)
    endif
    if(hdl .eq. "time") then
       ALLOCATE(data(1,1,1,1), STAT = alstatus)
       ALLOCATE (dms(1), STAT = alstatus)
       dms = 1
       CALL H5Dread_f(did,tp,data,dms,status)

       !-----------CLOSE ALL HANDLERS-----------------------!
       CALL H5Dclose_f(did, status)
       CALL H5Gclose_f(mgr, status)
       CALL H5Gclose_f(ugr, status)
       CALL H5Fclose_f(file, status)
       CALL H5close_f(status)
       return
    endif

    !open datasets and find the dimensions. Then allocate the 
    !array to put the data in. After that it is sent back to the
    !set data routine and read in to the modl strct.
    CALL H5Dget_space_f(did, sid, status)
    CALL H5Sget_simple_extent_ndims_f(sid,ndims,status)
    ALLOCATE (dms(ndims), STAT = alstatus)
    ALLOCATE (dmsm(ndims), STAT = alstatus)
    CALL H5Sget_simple_extent_dims_f(sid,dms,dmsm,status)
    !Get allocate data to return
    if(ndims .EQ. 3) then
       ALLOCATE ( data(dms(1),dms(2),dms(3),1), STAT = alstatus)
    elseif(ndims .EQ. 4) then
       ALLOCATE ( data(dms(1),dms(2),dms(3),dms(4)), STAT = alstatus)
    elseif(ndims .EQ. 1) then
       ALLOCATE ( data(dms(1),1,1,1), STAT = alstatus)
    elseif(ndims .EQ. 0) then
       ALLOCATE ( data(1,1,1,1), STAT = alstatus)
    elseif(ndims .EQ. 2) then
       ALLOCATE ( data(dms(1),dms(2),1,1), STAT = alstatus)
    else
       write(*,*) 'DIMENSION PROBLEMS'
       CALL EXIT(STATUS)
    endif
    CALL H5Dread_f(did,tp,data,dms,status) !DATA READ

    !-----------CLOSE ALL HANDLERS-----------------------!
    CALL H5Sclose_f(sid, status)
    CALL H5Dclose_f(did, status)
    if(yy .EQ. 1 .OR. yy .EQ. 2) then
       CALL H5Gclose_f(mgryy,status)
    endif
    CALL H5Gclose_f(mgr, status)
    CALL H5Gclose_f(ugr, status)

    CALL H5Fclose_f(file, status)

    CALL H5close_f(status)
  end subroutine readdata

  subroutine cpmodl(modl,cpy,i1,i2)
    TYPE(fruit), INTENT(IN) :: modl
    TYPE(fruit), INTENT(OUT) :: cpy
    INTEGER, INTENT(IN) :: i1,i2

    cpy%time = modl%time
    cpy%xzn = modl%xzn(i1:i2)
    cpy%xzr = modl%xzr(i1:i2)
    cpy%xzl = modl%xzl(i1:i2)
    cpy%yzn = modl%yzn
    cpy%yzr = modl%yzr
    cpy%yzl = modl%yzl
    cpy%zzn = modl%zzn
    cpy%zzr = modl%zzr
    cpy%zzl = modl%zzl
    cpy%den = modl%den(i1:i2,:,:)
    cpy%gpo = modl%gpo(i1:i2,:,:)
    cpy%vex = modl%vex(i1:i2,:,:)
    cpy%vey = modl%vey(i1:i2,:,:)
    cpy%vez = modl%vez(i1:i2,:,:)
    cpy%yyw = modl%yyw
  end subroutine cpmodl


  subroutine setdata(filename,j,modle,sp,yy)
    IMPLICIT NONE
    TYPE(fruit),INTENT(OUT) :: modle
    INTEGER, INTENT(IN) :: yy
    CHARACTER(LEN=90),INTENT(IN) :: filename
    INTEGER(hsize_t):: s,l
    CHARACTER (len=:), allocatable :: dum
    INTEGER(hsize_t),INTENT(IN):: j
    INTEGER(hsize_t),INTENT(OUT):: sp
    INTEGER :: i,ck
    REAL(KIND=8), dimension(:,:,:,:),ALLOCATABLE :: data(:,:,:,:)
    INTEGER(hsize_t), dimension(:),ALLOCATABLE :: dms(:),dmsm(:)
    dum = "time"
    CALL readdata(filename,dum,j,sp,data,0)
    modle%time = data(1,1,1,1)
    DEALLOCATE(data, STAT = ck)

    dum = "xzn"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%xzn(size(data(:,1,1,1))), STAT = ck)
    modle%xzn = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "xzr"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%xzr(size(data(:,1,1,1))), STAT = ck)
    modle%xzr = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "xzl"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%xzl(size(data(:,1,1,1))), STAT = ck)
    modle%xzl = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "yzn"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%yzn(size(data(:,1,1,1))), STAT = ck)
    modle%yzn = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "yzr"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%yzr(size(data(:,1,1,1))), STAT = ck)
    modle%yzr = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "yzl"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%yzl(size(data(:,1,1,1))), STAT = ck)
    modle%yzl = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "zzn"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%zzn(size(data(:,1,1,1))), STAT = ck)
    modle%zzn = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "zzr"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%zzr(size(data(:,1,1,1))), STAT = ck)
    modle%zzr = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "zzl"
    CALL readdata(filename,dum,j,sp,data,0)
    ALLOCATE(modle%zzl(size(data(:,1,1,1))), STAT = ck)
    modle%zzl = data(:,1,1,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "den"
    CALL readdata(filename,dum,j,sp,data,yy)
    ALLOCATE(modle%den(size(data(:,1,1,1)),size(data(1,:,1,1)),size(data(1,1,:,1))), STAT = ck)
    modle%den = data(:,:,:,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "gpo"
    CALL readdata(filename,dum,j,sp,data,yy)
    ALLOCATE(modle%gpo(size(data(:,1,1,1)),size(data(1,:,1,1)),size(data(1,1,:,1))), STAT = ck)
    modle%gpo = data(:,:,:,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "vex"
    CALL readdata(filename,dum,j,sp,data,yy)
    ALLOCATE(modle%vex(size(data(:,1,1,1)),size(data(1,:,1,1)),size(data(1,1,:,1))), STAT = ck)
    modle%vex = data(:,:,:,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "vey"
    CALL readdata(filename,dum,j,sp,data,yy)
    ALLOCATE(modle%vey(size(data(:,1,1,1)),size(data(1,:,1,1)),size(data(1,1,:,1))), STAT = ck)
    modle%vey = data(:,:,:,1) 
    DEALLOCATE(data, STAT = ck)

    dum = "vez"
    CALL readdata(filename,dum,j,sp,data,yy)
    ALLOCATE(modle%vez(size(data(:,1,1,1)),size(data(1,:,1,1)),size(data(1,1,:,1))), STAT = ck)
    modle%vez = data(:,:,:,1) 
    DEALLOCATE(data, STAT = ck)
    if(yy .EQ. 1) then
       dum = "yinyang_weight"
       CALL readdata(filename,dum,j,sp,data,0)
       ALLOCATE(modle%yyw(size(modle%yzn),size(modle%zzn)),STAT=ck)
       modle%yyw = data(:,:,1,1) 
       DEALLOCATE(data, STAT = ck)

    elseif(yy .EQ. 2) then
       dum = "yinyang_weight"
       CALL readdata(filename,dum,j,sp,data,0)
       ALLOCATE(modle%yyw(size(modle%yzn),size(modle%zzn)),STAT=ck)
       modle%yyw = data(:,:,1,1) 
       DEALLOCATE(data, STAT = ck)
    else
       ALLOCATE(modle%yyw(size(modle%yzn),size(modle%zzn)),STAT=ck)
       modle%yyw = 1.0d0
    end if

  end subroutine setdata
end module dataread
