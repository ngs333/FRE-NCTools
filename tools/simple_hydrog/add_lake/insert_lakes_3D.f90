program insert_lakes_3D

! program inserts lake-fractions into m45 3D cover type, ground type files 

implicit none

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx= 1000

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20

include '/usr/local/include/netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, idx_lake, nlake, rcode2, ntype, zdim, zid
integer :: jdg, idg, jdx, idx
real :: pi, dtr, sum, mval_frac, casp_ratio
real*4 :: mval_out
character(len=8)   :: var
character(len=40)  :: var_units, zax_name
character(len=100) :: long_name, input_cover_file, lake_file, gridspec_file

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, ncell_lake, lake_type

real, dimension (nlkmx)                :: lat_lake, lon_lake, lake_area

real, allocatable, dimension (:)       :: lat, latb, lon, lonb, zdat
real, allocatable, dimension (:,:)     :: cell_area, latg, long
real, allocatable, dimension (:,:,:)   :: frac, frac_new

real*4, allocatable, dimension (:,:,:)     :: adat3

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.insert_lakes_3D', form= 'formatted')

casp_ratio= 1.
read (5,'(a)') input_cover_file
read (5,*) idx_lake
read (5,'(a)') gridspec_file
read (5,*) casp_ratio
close (5)

! ----------------------------------------------------------------------
! read lat and lon of lake to be inserted in cover field
!   lake_type is 1 for internal drain points -- where we'll assign
!     a lake fraction = 0.2; associated lake_area is meaningless
!   lake_type is 2 for lakes draining to the ocean; lake_area is
!     the area of lake in the grid cell
! ----------------------------------------------------------------------

write (10,*) 'input lake info:'
open (20, form= 'formatted')
ktr= 0
1 ktr= ktr + 1
    read (20,*,end=2) lat_lake(ktr), lon_lake(ktr), lake_type(ktr), lake_area(ktr)
    if (lon_lake(ktr) < 0.) lon_lake(ktr)= lon_lake(ktr) + 360.
    write (10,'(2f10.2,i6,f12.0)') lat_lake(ktr), lon_lake(ktr), &
        lake_type(ktr), lake_area(ktr)
    go to 1
2 continue
nlake= ktr - 1
write (6,*) 'number of lakes to be inserted = ', nlake
close (20)


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
    rcode2 = nf_inq_varid (ncid, 'lon', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 23
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
write (6,*) 'id= ', id

allocate (lon(id), cell_area(id,jd))

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

! compute area of latitude
cell_area= 0.
do j= 1,jd
   do i= 1,id
      cell_area(i,j)= erad*erad*(lonb(i+1)-lonb(i))*dtr* &
          (sin(latb(j+1)*dtr) - sin(latb(j)*dtr))
   enddo
enddo

write (6,*) 'erad= ', erad, ', dtr= ', dtr
sum= 0.
do j= 1,jd
   do i= 1,id
      sum= sum + cell_area(i,j)
   enddo
enddo
write (6,*) 'sum= ', sum

write (10,'(/"input file lats")')
do j= 1,jd
   write (10,'(i5,6f12.4)') j, lat(j), latb(j), latb(j + 1), cell_area(1,j)
enddo

write (10,'(/"input file lons")')
do i= 1,id
   write (10,'(i5,4f12.4)') i, lon(i), lonb(i), lonb(i + 1)
enddo
write (10,'(/)')


! ----------------------------------------------------------------------
! get z values (cover types)
! ----------------------------------------------------------------------

write (zax_name,'(a,i2)') 'zax1_', idx_lake
rcode= nf_inq_varid (ncid, trim(zax_name), zid)         ! number of lats
if (rcode /= 0) then
    write (zax_name,'(a,i2)') 'ZAX1_', idx_lake
    rcode2 = nf_inq_varid (ncid, trim(zax_name), zid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find z variable" ; stop 26
    endif
endif
rcode= nf_inq_vardimid (ncid, zid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ntype)
if (ntype /= idx_lake) then
    write (6,*) 'ERROR: inconsistent dimensions, ntype= ', ntype, ', idx_lake= ', idx_lake
    stop 27
endif

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

allocate (frac(id,jd,ntype))

frac= mval_frac
start= 1; count= 1; count(1)= id ; count(2)= jd ; count(3)= ntype
rcode= nf_get_vara_double (ncid, varid, start, count, frac)

rcode= nf_close (ncid)


! ----------------------------------------------------------------------
! okay now overwrite lats and lons with those from input gridspec file
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(gridspec_file), NF_NOWRITE, ncid)
if  (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"  ; stop 12
endif

start= 1 ;  count= 1

rcode= nf_inq_varid (ncid, 'Y', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'y', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable, gridspec file" ; stop 122
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), idg)
rcode= nf_inq_dimlen (ncid, dimids(2), jdg)
write (6,*) 'gridspec, x, jdg= ', jdg, ', idg= ', idg
if (jdg /= (jd*2+1)) then
    write (6,*) 'ERROR: gridspec has inconsistent dimensions, y, jdg= ', jdg, ', jd= ', jd
    stop 131
endif
if (idg /= (id*2+1)) then
    write (6,*) 'ERROR: gridspec has inconsistent dimensions, x, idg= ', idg, ', id= ', id
    stop 132
endif

allocate (latg(idg,jdg))
count(1)= idg ;  count(2)= jdg
rcode= nf_get_vara_double (ncid, latid, start, count, latg)

if (latg(1,1) > 0.) then
    write (6,*) "ERROR: lats in wrong direction"
    stop 110
endif

do j= 1,jdg-1,2
   jdx= (j+1)/2
   lat(jdx)= latg(1,j+1)
   latb(jdx)= latg(1,j)
enddo
latb(jdp1)= latg(1,jdg)

rcode= nf_inq_varid (ncid, 'X', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable, gridspec file" ; stop 123
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), idg)
rcode= nf_inq_dimlen (ncid, dimids(2), jdg)
write (6,*) 'gridspec, y, jdg= ', jdg, ', idg= ', idg
if (jdg /= (jd*2+1)) then
    write (6,*) 'ERROR: gridspec has inconsistent dimensions, y, jdg= ', jdg, ', jd= ', jd
    stop 141
endif
if (idg /= (id*2+1)) then
    write (6,*) 'ERROR: gridspec has inconsistent dimensions, x, idg= ', idg, ', id= ', id
    stop 142
endif

allocate (long(idg,jdg))
count(1)= idg ;  count(2)= jdg
rcode= nf_get_vara_double (ncid, lonid, start, count, long)

if (long(1,1) < 0.) then
    write (6,*) "ERROR: lons in wrong direction"
    stop 111
endif

do i= 1,idg-1,2
   idx= (i+1)/2
   lon(idx)= long(i+1,1)
   lonb(idx)= long(i,1)
enddo
lonb(idp1)= long(idg,1)

deallocate (latg, long)

! compute area of latitude
cell_area= 0.
do j= 1,jd
   do i= 1,id
      cell_area(i,j)= erad*erad*(lonb(i+1)-lonb(i))*dtr* &
          (sin(latb(j+1)*dtr) - sin(latb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jd
   do i= 1,id
      sum= sum + cell_area(i,j)
   enddo
enddo
write (6,*) 'sum= ', sum

write (10,'(/"gridspec file lats")')
do j= 1,jd
   write (10,'(i5,6f12.4)') j, lat(j), latb(j), latb(j + 1), cell_area(1,j)
enddo

write (10,'(/"gridspec file lons")')
do i= 1,id
   write (10,'(i5,4f12.4)') i, lon(i), lonb(i), lonb(i + 1)
enddo
write (10,'(/)')

rcode= nf_close (ncid)


! ----------------------------------------------------------------------
!  find i,j of cover grid cell in which input lake point is located
! ----------------------------------------------------------------------
do l= 1,nlake
   do j= 1,jd
      if (lat_lake(l) >= latb(j) .and. lat_lake(l) < latb(j+1)) then
          jlake(l)= j
          go to 10
      endif
   enddo
10 continue

   do i= 1,id
      if (lon_lake(l) >= lonb(i) .and. lon_lake(l) < lonb(i+1)) then
          ilake(l)= i
          go to 20
      endif
   enddo
20 continue
enddo

   
! ----------------------------------------------------------------------
!  now insert lake fraction in frac_new field
!    do not alter Caspian fractions -- already present in frac field
! ----------------------------------------------------------------------

allocate (frac_new(id,jd,ntype))

write (10,'(/"insert lakes:")')
frac_new= frac

write (6,*) 'Caspian drain point, frac= ', frac(21,67,idx_lake)

! multiply all Caspian lake-area by pre-determined factor; do this
!  only for Caspian, the only lake present in the frac field
do j= 1,jd
   do i= 1,id
      if (frac(i,j,idx_lake) > 0.) then
          frac_new(i,j,idx_lake)= frac(i,j,idx_lake)*casp_ratio
          if (frac_new(i,j,idx_lake) > 1.) then
              write (6,'(a,f15.6,2i6,2f12.2)') 'ERROR: Casp, frac > 1, ', &
                 frac_new(i,j,idx_lake), j, i, lat(j), lon(i)
          endif
      endif
   enddo
enddo

write (6,*) 'Caspian drain point, frac_new= ', frac_new(21,67,idx_lake)

!do l= 1,nlake
!   if (frac(ilake(l),jlake(l),idx_lake) > 0.) then
!       write (6,*) "lake already present at ", lon_lake(l), lat_lake(l)
!   else
!       if (lake_type(l) == 1) then
!           frac_new(ilake(l),jlake(l),idx_lake)= lfrac_const
!       else if (lake_type(l) == 2) then
!           frac_new(ilake(l),jlake(l),idx_lake)= &
!             lake_area(l)/cell_area(ilake(l),jlake(l))
!           write (10,'(i4,2f10.2,2f12.0,f10.2)') l, lat_lake(l), lon_lake(l), lake_area(l), &
!             cell_area(ilake(l),jlake(l)), lake_area(l)/cell_area(ilake(l),jlake(l))
!           if (lake_area(l) > cell_area(ilake(l),jlake(l))) then
!               write (6,'(a,2f10.2,2f12.0)') "ERROR: lake area exceeds grid cell area, ", &
!                 lon_lake(l), lat_lake(l), lake_area(l), cell_area(ilake(l),jlake(l))
!!               stop 40
!           endif
!       else
!           write (6,*) "ERROR: invalid lake_type, ", l, lake_type(l)
!           stop 41
!       endif
!   endif
!enddo


do l= 1,nlake
   if (lake_type(l) == 1) then
       if (frac(ilake(l),jlake(l),idx_lake) > 0.) then
           write (6,*) "lake already present at ", lon_lake(l), lat_lake(l)
       else
           frac_new(ilake(l),jlake(l),idx_lake)= lfrac_const
       endif
   else if (lake_type(l) == 2) then
       if (frac(ilake(l),jlake(l),idx_lake) > 0.) write (6,*) "lake already present at ", &
           lon_lake(l), lat_lake(l)
       frac_new(ilake(l),jlake(l),idx_lake)= lake_area(l)/cell_area(ilake(l),jlake(l))
       
       write (10,'(i4,2f10.2,2f12.0,f15.6)') l, lat_lake(l), lon_lake(l), lake_area(l), &
         cell_area(ilake(l),jlake(l)), lake_area(l)/cell_area(ilake(l),jlake(l))
       if (lake_area(l) > cell_area(ilake(l),jlake(l))) then
           write (6,'(a,2f10.2,2f12.0)') "ERROR: lake area exceeds grid cell area, ", &
             lon_lake(l), lat_lake(l), lake_area(l), cell_area(ilake(l),jlake(l))
!           stop 40
       endif
   else
       write (6,*) "ERROR: invalid lake_type, ", l, lake_type(l)
       stop 41
   endif
enddo

write (6,*) 'Caspian drain point, frac_new= ', frac_new(21,67,idx_lake)


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

write (zax_name,'(a,i2)') 'zax1_', idx_lake
rcode= NF_DEF_DIM (ncid, trim(zax_name), idx_lake, zdim)

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

deallocate (lat, latb, lon, lonb, zdat, cell_area)
deallocate (frac, frac_new)
deallocate (adat3)
   
   
stop

end



