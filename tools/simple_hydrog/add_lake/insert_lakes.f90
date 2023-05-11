program insert_lakes

! program inserts 1-cell lakes in 0.5x0.5 or 1x1 cover type, ground type files

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: nlkmx= 1000
integer, parameter :: ncellmx= 12

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.
real, parameter :: res_m45= 2.5

!!TODO:
!!include '/usr/local/include/netcdf.inc'
include 'netcdf.inc'


integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, idx_lake, nlake, ncl_fill, nc, mval_cover
integer :: bmin, bmax, rcode2, id_lake, jd_lake, i1, j1, nip, njp
integer :: input_type, ii, jj
real :: pi, dtr, sum, reslon, res
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: long_name, input_cover_file, lake_file

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, ncell_lake

real, dimension (nlkmx)                :: lat_lake, lon_lake, chk_lake

integer, dimension (ncellmx)           :: jinc= &
(/ 0,  0, -1, -1,  0, -1,  1,  1,  1,  1,  0, -1 /)

integer, dimension (ncellmx)           :: iinc= &
(/ 0,  1,  0,  1, -1, -1, -1,  0,  1,  2,  2,  2 /)

integer, allocatable, dimension (:,:)  :: cover, cover_lake
real, allocatable, dimension (:)       :: lat, latb, lon, lonb, arlat, &
                                          lats_lk, lons_lk
real, allocatable, dimension (:,:)     :: lake_dat

real*4, allocatable, dimension (:)     :: adat1

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.insert_lakes', form= 'formatted')

read (5,'(a)') input_cover_file
read (5,*) idx_lake
read (5,*) input_type
if (input_type == 2) then
   read (5,'(a)') lake_file
endif
close (5)

! ----------------------------------------------------------------------
! read lat land lon of lake to be inserted in cover field
! ----------------------------------------------------------------------

if (input_type == 1) then
   open (20, form= 'formatted')
   ktr= 0
   1 ktr= ktr + 1
     read (20,*,end=2) lat_lake(ktr), lon_lake(ktr)
     if (lon_lake(ktr) < 0.) lon_lake(ktr)= lon_lake(ktr) + 360.
     go to 1
2    continue
   nlake= ktr - 1
   write (6,*) 'number of lakes to be inserted = ', nlake
   close (20)
else if (input_type == 2) then
!  get lats and lons of input data
   rcode= NF_OPEN (trim(lake_file), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file"  ; stop 1
   endif

   start= 1 ;  count= 1

   rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find lat variable" ; stop 2
       endif
   endif
   rcode= nf_inq_vardimid (ncid, latid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), jd_lake)
   write (6,*) 'jd_lake= ', jd_lake

   allocate (lats_lk(jd_lake))

   count(1)= jd_lake
   rcode= nf_get_vara_double (ncid, latid, start, count, lats_lk)

   if (lats_lk(1) > 0.) then
       write (6,*) "ERROR: lats in wrong direction"
       stop 2
   endif

   rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find lon variable" ; stop 3
       endif
   endif
   rcode= nf_inq_vardimid (ncid, lonid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), id_lake)
   write (6,*) 'id_lake= ', id_lake

   allocate (lons_lk(id_lake))

   count(1)= id_lake
   rcode= nf_get_vara_double (ncid, lonid, start, count, lons_lk)

   if (lons_lk(1) < 0.) then
       write (6,*) "ERROR: lons in wrong direction"
       stop 3
   endif

   write (10,'(/"lake file lats")')
   do j= 1,jd_lake
      write (10,'(i5,6f12.4)') j, lats_lk(j)
   enddo

   write (10,'(/"lake file lons")')
   do i= 1,id_lake
      write (10,'(i5,4f12.4)') i, lons_lk(i)
   enddo
   write (10,'(/)')

   rcode= nf_close (ncid)

! now read grid of data indicating which grid cells need more lake cells
   allocate (lake_dat(id_lake,jd_lake))

   open (20, form= 'formatted')
   read (20,*) i1, nip
   read (20,*) j1, njp
   do i= 1,4
      read (20,'(10x)')
   enddo
   lake_dat= -1.
   do j= j1,j1+njp-1
      read (20,*) (lake_dat(i,j), i= i1,i1+nip-1)
   enddo
   close (20)

! for this first time at least, all but 2 grid cells have "lake_dat" < 5
   ktr= 0
   do j= 1,jd_lake
      do i= 1,id_lake
         if (lake_dat(i,j) < 0.0) go to 15
         if (lake_dat(i,j) > real(ncellmx-1)) go to 15   ! don't use the really big ones
         do n= 1,50
            if (lake_dat(i,j) >= real(n-1) .and. lake_dat(i,j) < real(n)) then
                ktr= ktr + 1
                lat_lake(ktr)= lats_lk(j)
                lon_lake(ktr)= lons_lk(i)
                chk_lake(ktr)= lake_dat(i,j)
                ncell_lake(ktr)= int(lake_dat(i,j)+1.)
                go to 15
            endif
         enddo
         write (6,*) "ERROR: lake_dat not classed, ", j, i, lake_dat(i,j)
         stop 5
15       continue
      enddo
   enddo
   nlake= ktr

   write (6,*) 'number of lakes in grid= ', nlake
   do n= 1,ktr
      write (6,'(i6,2f9.2,f10.4,i6)') n, lat_lake(n), lon_lake(n), chk_lake(n), ncell_lake(n)
   enddo

   deallocate (lats_lk, lons_lk, lake_dat)
else
   write (6,*) "ERROR: invalid input_type, ", input_type
   stop 6
endif


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
jdp1= jd + 1
write (6,*) 'jd= ', jd, ', jdp1= ', jdp1

allocate (lat(jd), latb(jdp1), arlat(jd))

count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat)

if (lat(1) > 0.) then
    write (6,*) "ERROR: lats in wrong direction"
    stop 10
endif

latb(1)= lat1
latb(jdp1)= lat2
do j= 2,jd
   latb(j)= 0.5*(lat(j-1) + lat(j))
enddo

rcode= nf_inq_varid (ncid, 'LON', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'lon', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 23
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
idp1= id + 1
write (6,*) 'id= ', id, ', idp1= ', idp1

allocate (lon(id), lonb(idp1))

count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon)

if (lon(1) < 0.) then
    write (6,*) "ERROR: lons in wrong direction"
    stop 11
endif

lonb(1)= lon1
lonb(idp1)= lon2
do i= 2,id
   lonb(i)= 0.5*(lon(i-1) + lon(i))
enddo

reslon= 360./real(id)
write (6,*) "reslon= ", reslon
! compute area of latitude
arlat= 0.
do j= 1,jd
   arlat(j)= erad*erad*reslon*dtr*(sin(latb(j+1)*dtr) - sin(latb(j)*dtr))
enddo

sum= 0.
do j= 1,jd
   sum= sum + arlat(j)*real(id)
enddo
write (6,*) 'sum= ', sum

write (10,'(/"input file lats")')
do j= 1,jd
   write (10,'(i5,6f12.4)') j, lat(j), latb(j), latb(j + 1), arlat(j)
enddo

write (10,'(/"input file lons")')
do i= 1,id
   write (10,'(i5,4f12.4)') i, lon(i), lonb(i), lonb(i + 1)
enddo
write (10,'(/)')


! ----------------------------------------------------------------------
!  now read cover field
! ----------------------------------------------------------------------

allocate (cover(id,jd))

write (6,*) 'read cover'
rcode= nf_inq_varid (ncid, 'cover', varid)
if (rcode /= 0) then
    print *, "rcode in ncvid is ", rcode ;  stop 10
endif

var_units= ' '
if (rcode == 0) then
    rcode= nf_get_att_text (ncid, varid, "units", var_units)
endif
write (6,*) 'units= ', trim(var_units)

long_name= ' '
rcode= nf_get_att_text (ncid, varid, "long_name", long_name)
write (6,*) 'long_name= ', trim(long_name)

mval_cover= -99999
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_int (ncid, varid, 'missing_value', mval_cover)
endif
write (6,*) 'mval= ', mval_cover

cover= mval_cover
start= 1; count= 1; count(1)= id ; count(2)= jd
rcode= nf_get_vara_int (ncid, varid, start, count, cover)

rcode= nf_close (ncid)

bmax= -99999. ;  bmin= 99999.
do j= 1,jd
   do i= 1,id
      if (cover(i,j) /= mval_cover) then
          bmax= max(cover(i,j),bmax)
          bmin= min(cover(i,j),bmin)
      endif
   enddo
enddo
write (6,*) 'bmin= ', bmin, ', bmax= ', bmax


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
!  now insert lake in cover field
! ----------------------------------------------------------------------

allocate (cover_lake(id,jd))

! this code is for creating 3x3 or 5x5 clusters of lake cells
!res= 360./real(id)

!if (amod(res_m45,res) == 0.) then
!    ncl_fill= int(res_m45/res)
!else
!    ncl_fill= int(res_m45/res) + 1
!endif

!if (mod(ncl_fill,2) == 0) ncl_fill= ncl_fill + 1

!nc= ncl_fill/2

!cover_lake= cover
!do l= 1,nlake
!   do j= jlake(l)-nc,jlake(l)+nc
!      do i= ilake(l)-nc,ilake(l)+nc
!         if (cover(i,j) /= mval_cover) then
!             cover_lake(i,j)= real(idx_lake)
!         endif
!      enddo
!   enddo
!enddo

if (input_type == 1) then
    cover_lake= cover
    do l= 1,nlake
       cover_lake(ilake(l),jlake(l))= real(idx_lake)
    enddo
else
    cover_lake= cover
    do l= 1,nlake
       do n= 1,ncell_lake(l)
          ii= ilake(l)+iinc(n)
          jj= jlake(l)+jinc(n)
          if (ii < 1) ii= id + ii
          if (ii > id) ii= ii - id
          cover_lake(ii,jj)= real(idx_lake)
       enddo
    enddo
endif


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

rcode= NF_CREATE ('cover_file_lake.nc', NF_CLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', 17, 'cover_file_lake.nc')

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lonb, lat, latb)
rcode= NF_DEF_DIM (ncid, 'lon', id, londim)
rcode= NF_DEF_DIM (ncid, 'lat', jd, latdim)

!  create coordinate variables
rcode= NF_DEF_VAR (ncid, 'lon', NF_DOUBLE, 1, londim, lonid)
rcode= NF_DEF_VAR (ncid, 'lat', NF_DOUBLE, 1, latdim, latid)

!  create attributes for coordinate variables
!    longitude:
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 9, 'longitude')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'cartesian_axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')

!    latitude:
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'cartesian_axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')

!    create data variable and attributes
ndims(1)= londim
ndims(2)= latdim
rcode= NF_DEF_VAR (ncid, 'cover', NF_INT, 2, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', len_trim(long_name), trim(long_name))
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
rcode= NF_PUT_ATT_INT (ncid, varid, 'missing_value', NF_INT, 1, mval_cover)

!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data
start= 1 ;  count= 1

count(1)= id
rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start(1), count(1), lon)

count(1)= jd
rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start(1), count(1), lat)

!    new cover data with inserted lakes
count(1)= id ;  count(2)= jd
rcode= NF_PUT_VARA_INT (ncid, varid, start, count, cover_lake)

!  close netcdf file
rcode= NF_CLOSE (ncid)

deallocate (lat, latb, lon, lonb, arlat)
deallocate (cover, cover_lake)


stop

end



