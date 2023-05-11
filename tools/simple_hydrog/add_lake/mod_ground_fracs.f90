program mod_ground_fracs

! program inserts lake-fractions into m45 3D cover type, ground type files 

implicit none

integer, parameter :: maxdims= 3

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, nlake, rcode2, ntype, zdim, zid
real :: pi, dtr, sum, mval_frac
real*4 :: mval_out
character(len=8)   :: var
character(len=40)  :: var_units, zax_name
character(len=100) :: long_name, input_cover_file

integer, dimension (maxdims)           :: start, count, dimids, ndims

real, allocatable, dimension (:)       :: lat, latb, lon, lonb, zdat
real, allocatable, dimension (:,:,:)   :: frac, frac_new

real*4, allocatable, dimension (:,:,:)     :: adat3

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.mod_ground_fracs', form= 'formatted')

read (5,'(a)') input_cover_file
close (5)

! ----------------------------------------------------------------------
! now open cover file
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(input_cover_file), NF_NOWRITE, ncid)
if  (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"  ; stop 12
endif

start= 1
count= 1

! ----------------------------------------------------------------------
! get latitudes and longitudes of input cover file
! ----------------------------------------------------------------------

rcode= nf_inq_varid (ncid, 'LAT', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'lat', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 22
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jd)
write (6,*) 'jd= ', jd

allocate (lat(jd))

count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat)

if (lat(1) > 0.) then
    write (6,*) "ERROR: lats in wrong direction"
    stop 10
endif

rcode= nf_inq_varid (ncid, 'LATB', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'latb', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find latb variable" ; stop 22
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jdp1)
    
allocate (latb(jdp1))
count(1)= jdp1
rcode= nf_get_vara_double (ncid, latid, start, count, latb)
write (6,*) 'jdp1= ', jdp1

rcode= nf_inq_varid (ncid, 'LON', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'lon', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 23
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
write (6,*) 'id= ', id

allocate (lon(id))

count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon)

if (lon(1) < 0.) then
    write (6,*) "ERROR: lons in wrong direction"
    stop 11
endif

rcode= nf_inq_varid (ncid, 'LONB', lonid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'lonb', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lonb variable" ; stop 23
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), idp1)
    
allocate (lonb(idp1))
count(1)= idp1
rcode= nf_get_vara_double (ncid, lonid, start, count, lonb)
write (6,*) 'idp1= ', idp1


write (10,'(/"input file lats")')
do j= 1,jd
   write (10,'(i5,6f12.4)') j, lat(j), latb(j), latb(j + 1)
enddo

write (10,'(/"input file lons")')
do i= 1,id
   write (10,'(i5,4f12.4)') i, lon(i), lonb(i), lonb(i + 1)
enddo
write (10,'(/)')


! ----------------------------------------------------------------------
! get z values (cover types)
! ----------------------------------------------------------------------

write (zax_name,'(a,i2)') 'zax1_10'
rcode= nf_inq_varid (ncid, trim(zax_name), zid)         ! number of lats
if (rcode /= 0) then
    write (zax_name,'(a,i2)') 'ZAX1_10'
    rcode2 = nf_inq_varid (ncid, trim(zax_name), zid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find z variable" ; stop 26
    endif
endif
rcode= nf_inq_vardimid (ncid, zid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ntype)

allocate (zdat(ntype))

count(1)= ntype
rcode= nf_get_vara_double (ncid, zid, start, count, zdat)

write (6,*) 'zdat= ', zdat


! ----------------------------------------------------------------------
!  now read cover field
! ----------------------------------------------------------------------

write (6,*) 'read frac'
rcode= nf_inq_varid (ncid, 'frac', varid)
if (rcode /= 0) then
    print *, "rcode in ncvid is ", rcode ;  stop 30
endif

rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= id) then
    write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', id= ', id
    stop 31
endif
rcode= nf_inq_dimlen (ncid, dimids(2), n)
if (n /= jd) then
    write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', jd= ', jd
    stop 32
endif
rcode= nf_inq_dimlen (ncid, dimids(3), n)
if (n /= ntype) then
    write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', ntype= ', ntype
    stop 33
endif


var_units= ' '
if (rcode == 0) then
    rcode= nf_get_att_text (ncid, varid, "units", var_units)
endif
write (6,*) 'units= ', trim(var_units)

long_name= ' '
rcode= nf_get_att_text (ncid, varid, "long_name", long_name)
write (6,*) 'long_name= ', trim(long_name)

mval_frac= -99999.
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_frac)
endif
write (6,*) 'mval= ', mval_frac

allocate (frac(id,jd,ntype), frac_new(id,jd,ntype))

frac= mval_frac
start= 1; count= 1; count(1)= id ; count(2)= jd ; count(3)= ntype
rcode= nf_get_vara_double (ncid, varid, start, count, frac)

rcode= nf_close (ncid)


! now modify frac such that z=9 (ice) data is added to z=2, z=9 set to zero
frac_new= frac
frac_new(:,:,2)=frac(:,:,2)+frac(:,:,9)
frac_new(:,:,9)=0.

! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

rcode= NF_CREATE ('cover_file_lake.nc', NF_CLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', 17, 'cover_file_lake.nc')
   
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lonb, lat, latb, cover type)
rcode= NF_DEF_DIM (ncid, 'lon',  id,   londim)
rcode= NF_DEF_DIM (ncid, 'lat',  jd,   latdim)
rcode= NF_DEF_DIM (ncid, 'lonb', idp1, lonbdim)
rcode= NF_DEF_DIM (ncid, 'latb', jdp1, latbdim)

write (zax_name,'(a,i2)') 'zax1_10'
rcode= NF_DEF_DIM (ncid, trim(zax_name), 10, zdim)

!  create coordinate variables
rcode= NF_DEF_VAR (ncid, 'lon',  NF_DOUBLE, 1, londim,  lonid)
rcode= NF_DEF_VAR (ncid, 'lat',  NF_DOUBLE, 1, latdim,  latid)
rcode= NF_DEF_VAR (ncid, 'lonb', NF_DOUBLE, 1, lonbdim, lonbid)
rcode= NF_DEF_VAR (ncid, 'latb', NF_DOUBLE, 1, latbdim, latbid)
rcode= NF_DEF_VAR (ncid, trim(zax_name), NF_DOUBLE, 1, zdim, zid)

!  create attributes for coordinate variables
!    longitude:
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 9, 'longitude')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'edges', 4, 'lonb')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'point_spacing', 4, 'even')

!    latitude:
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'edges', 4, 'latb')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'point_spacing', 6, 'uneven')

!    longitude edges:
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'long_name', 15, 'longitude edges')
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'units', 9, 'degrees_E')

!    latitude edges:
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'long_name', 14, 'latitude edges')
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'units', 9, 'degrees_N')

!    z axis
rcode= NF_PUT_ATT_TEXT (ncid, zid, 'axis', 1, 'Z')
rcode= NF_PUT_ATT_TEXT (ncid, zid, 'point_spacing', 4, 'even')

!    create data variable and attributes
mval_out= mval_frac
ndims(1)= londim ;  ndims(2)= latdim ;  ndims(3)= zdim
rcode= NF_DEF_VAR (ncid, 'frac', NF_FLOAT, 3, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', len_trim(long_name), trim(long_name))
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
rcode= NF_PUT_ATT_REAL (ncid, varid, 'missing_value', NF_FLOAT, 1, mval_out)
 
!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data
start= 1 ;  count= 1
      
count(1)= id
rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start(1), count(1), lon)

count(1)= jd
rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start(1), count(1), lat)

count(1)= idp1
rcode= NF_PUT_VARA_DOUBLE (ncid, lonbid, start(1), count(1), lonb)

count(1)= jdp1
rcode= NF_PUT_VARA_DOUBLE (ncid, latbid, start(1), count(1), latb)

count(1)= ntype
rcode= NF_PUT_VARA_DOUBLE (ncid, zid, start(1), count(1), zdat)

!    new cover data with inserted lakes
allocate (adat3(id,jd,ntype))
adat3= frac_new
count(1)= id ;  count(2)= jd ;  count(3)= ntype
rcode= NF_PUT_VARA_REAL (ncid, varid, start, count, adat3)

!  close netcdf file
rcode= NF_CLOSE (ncid)

deallocate (lat, latb, lon, lonb, zdat)
deallocate (frac, frac_new)
deallocate (adat3)
   
   
stop

end



