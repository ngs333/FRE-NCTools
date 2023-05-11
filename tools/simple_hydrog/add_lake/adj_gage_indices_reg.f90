program adj_gage_indices_reg

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: ngagemax= 500
integer, parameter :: len_bshort= 4

real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, ktr, rcode, rcode2, varid, ncid, attnum, id, jd
integer :: latid, lonid, ndm, j1, i1, n1, i2, j2, ngdim, ngid, igid
integer :: jgid, bgid, chdim, ngage, ngage_rg
real :: pi, dtr, sum, mval_in, reglt1, reglt2, regln1, regln2
character(len=4)   :: b1
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: river_input_file, gage_file, fname

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ngagemax)          :: igage, jgage, igage_rg, jgage_rg

real, dimension (ngagemax)             :: lat_rg, lon_rg

character(len=4), dimension (ngagemax) :: bshort, bshort_rg

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:)     :: lat, lon, land_frac, cell_area


pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.adj_gage_indices_reg', form= 'formatted')

read (5,'(a)') river_input_file
read (5,'(a)') gage_file
read (5,*) reglt1, reglt2
read (5,*) regln1, regln2
close (5)

write (6,*) 'fort.5 data read'


! ----------------------------------------------------------------------
! read gage data -- basin name; i,j,n of grid cell containing gage
! ----------------------------------------------------------------------

open (23, form= 'formatted')
ktr= 0 ; ngage= 0
21 ktr= ktr + 1
   read (23,'(a4,42x,3i6)',end=22) b1, i1, j1, n1
   if (n1 /= 1) then
       write (6,*) "ERROR: program cannot handle multiple tiles, n= ", n1
       write (6,*) "...exiting..."
       stop 21
   endif
   ngage= ngage + 1
   bshort(ngage)= b1
   igage (ngage)= i1
   jgage (ngage)= j1
   go to 21
22 continue
write (6,*)
write (6,*) 'number of gages = ', ktr - 1, ngage
close (23)



! ----------------------------------------------------------------------
!  get lon and lat dims from river file 
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(river_input_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open river netcdf file"  
    write (6,*) trim(river_input_file)
    stop 1
endif

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 2
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jd)
write (6,*) "jd= ", jd

allocate (lat_idx(jd))
start= 1 ;  count= 1 ;  count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat_idx)

       
rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 3
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
write (6,*) "id= ", id

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)
  
rcode= nf_close (ncid)


! ----------------------------------------------------------------------
! now open river file -- read lat,lon grids, land_frac, cellarea
! ----------------------------------------------------------------------

allocate (lat(id,jd), lon(id,jd), land_frac(id,jd), cell_area(id,jd))

!write (6,'(a)')  trim(river_input_file)
write (10,'(/a)') trim(river_input_file)

rcode= NF_OPEN (trim(river_input_file), NF_NOWRITE, ncid)
if  (rcode /= 0) then
     write (6,*) "ERROR: cannot open netcdf file"  ; stop 12
endif

start= 1 ;  count= 1

! ----------------------------------------------------------------------
! get latitudes and longitudes of input river file
! ----------------------------------------------------------------------

rcode= nf_inq_varid (ncid, 'y', latid)         ! lat field
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find lat variable (y)" ; stop 20
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= id) then
    write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 25
endif
rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
if (ndm /= jd) then
    write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 25
endif

start= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat)
       
       
rcode= nf_inq_varid (ncid, 'x', lonid)         ! lon field
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find lon variable (x)" ; stop 30
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= id) then
    write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 35
endif
rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
if (ndm /= jd) then
    write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 35
endif

start= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_double (ncid, lonid, start, count, lon)


!write (6,*) 'read land_frac'
rcode= nf_inq_varid (ncid, 'land_frac', varid)     ! land fraction field
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find land_frac" ; stop 50
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= id) then
    write (6,*) "ERROR: inconsistent land_frac dimension, ", ndm, id ;  stop 55
endif
rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
if (ndm /= jd) then
    write (6,*) "ERROR: inconsistent land_frac dimension, ", ndm, jd ;  stop 55
endif

start= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_double (ncid, varid, start, count, land_frac)

var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
!write (6,*) 'units= ', var_units

mval_in= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
endif
!write (6,*) 'mval= ', mval_in

where (land_frac == mval_in) land_frac= mval_mdl


!write (6,*) 'read cellarea'
rcode= nf_inq_varid (ncid, 'cellarea', varid)     ! cell area field
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find cellarea" ; stop 60
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= id) then
    write (6,*) "ERROR: inconsistent cellarea dimension, ", ndm, id ;  stop 65
endif
rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
if (ndm /= jd) then
    write (6,*) "ERROR: inconsistent cellarea dimension, ", ndm, jd ;  stop 65
endif

start= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_double (ncid, varid, start, count, cell_area)

var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
!write (6,*) 'units= ', var_units

mval_in= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
endif
!write (6,*) 'mval= ', mval_in

where (cell_area == mval_in) cell_area= 0.


rcode= nf_close (ncid)
   
!write (10,'(/"river lats")')
!do j= 1,jd
!   write (10,*) 'j= ', j
!   write (10,'(10f10.4)') (lat(i,j), i= 1,id)
!enddo

!write (10,'(/"river lons")')
!do j= 1,jd
!   write (10,*) 'j= ', j
!   write (10,'(10f10.4)') (lon(i,j), i= 1,id)
!enddo

!write (10,'(/"cell_area")')
!do j= 1,jd
!   write (10,*) 'j= ', j
!   write (10,'(10f10.1)') (cell_area(i,j)/1.e6, i= 1,id)
!enddo


write (10,*) 'gage data:'
do l= 1,ngage
   write (10,'(a,2i6,2f10.2)') bshort(l), igage(l), jgage(l), lat(igage(l),jgage(l)), lon(igage(l),jgage(l))
enddo


sum= 0.
do j= 1,jd
   do i= 1,id
      sum= sum + cell_area(i,j)
   enddo
enddo
write (6,'(/"sum of cell_area= ",f15.1)') sum/1.e+6

do j= 1,jd
   if (abs(lat(1,j)-reglt1) < 0.001) j1= j
   if (abs(lat(1,j)-reglt2) < 0.001) j2= j
enddo

do i= 1,id
   if (abs(lon(i,1)-regln1) < 0.001) i1= i
   if (abs(lon(i,1)-regln2) < 0.001) i2= i
enddo

ktr= 0
do l= 1,ngage
   if ( igage(l) >= i1 .and. igage(l) <= i2 .and. &
        jgage(l) >= j1 .and. jgage(l) <= j2 ) then
         ktr= ktr + 1
         igage_rg(ktr)= igage(l)-i1+1
         jgage_rg(ktr)= jgage(l)-j1+1
         bshort_rg(ktr)= bshort(l)
         lat_rg(ktr)= lat(1,jgage(l))
         lon_rg(ktr)= lon(igage(l),1)
   endif
enddo
ngage_rg= ktr
   
write (6,*) "ngage_rg= ", ngage_rg
do l= 1,ngage_rg
   write (6,'(a,2i6,2f10.2)') bshort_rg(l), igage_rg(l), jgage_rg(l), lat_rg(l), lon_rg(l)
enddo
   
write (10,*) "ngage_rg= ", ngage_rg
do l= 1,ngage_rg
   write (10,'(a,2i6,2f10.2)') bshort_rg(l), igage_rg(l), jgage_rg(l), lat_rg(l), lon_rg(l)
enddo


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

write (fname, '(a)') 'gage_data.nc'
rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))
   
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

 ! define character-position dimension
rcode= NF_DEF_DIM (ncid, 'ngage',  ngage_rg, ngdim)
rcode= NF_DEF_DIM (ncid, 'nchar',  len_bshort, chdim)

!    gages:
ndims(1)= ngdim
rcode= NF_DEF_VAR (ncid, 'igageloc',  NF_DOUBLE, 1, ngdim,  igid)
rcode= NF_PUT_ATT_TEXT (ncid, igid, 'long_name', 20, 'i index of gage cell')
rcode= NF_PUT_ATT_TEXT (ncid, igid, 'units', 4, 'none')

rcode= NF_DEF_VAR (ncid, 'jgageloc',  NF_DOUBLE, 1, ngdim,  jgid)
rcode= NF_PUT_ATT_TEXT (ncid, jgid, 'long_name', 20, 'j index of gage cell')
rcode= NF_PUT_ATT_TEXT (ncid, jgid, 'units', 4, 'none')

ndims(1)= chdim ;  ndims(2)= ngdim
rcode= NF_DEF_VAR (ncid, 'bname_gage',  NF_CHAR, 2, ndims,  bgid)
rcode= NF_PUT_ATT_TEXT (ncid, bgid, 'long_name', 23, 'basin name at gage cell')
rcode= NF_PUT_ATT_TEXT (ncid, bgid, 'units', 4, 'none')

!  leave define mode
rcode= NF_ENDDEF (ncid)

start= 1 ;  count= 1 ;  count(1)= ngage_rg
rcode= NF_PUT_VARA_DOUBLE (ncid, igid, start, count, real(igage_rg))

start= 1 ;  count= 1 ;  count(1)= ngage_rg
rcode= NF_PUT_VARA_DOUBLE (ncid, jgid, start, count, real(jgage_rg))

do l= 1,ngage_rg
   start(1)= 1 ;  count(1)= len_bshort
   start(2)= l ;  count(2)= 1
   rcode= NF_PUT_VARA_TEXT (ncid, bgid, start, count, bshort_rg(l))
enddo

!  close netcdf file
rcode= NF_CLOSE (ncid)

close (10)

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, land_frac, cell_area)

   
stop

end



