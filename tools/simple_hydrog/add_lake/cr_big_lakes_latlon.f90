program cr_big_lakes_latlon

! =========================================================================
!   program reads lat, lon, cellarea, land_frac, and fields, and
!     computes the fields: celllength, internal
! =========================================================================

implicit none

include "param.h"

integer, parameter :: maxdims= 2
integer, parameter :: ntilmx= 6
!integer, parameter :: idr= 5760, jdr= 2880
!!TODO: where shoudl idr, ijr come from? param.h?


integer, parameter :: idrp1= idr+1, jdrp1= jdr+1

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, k, rcode, varid, ncid, attnum, id, jd, ntiles
integer :: latid, lonid, latdim, londim, ii, jj, varid2, rcode2, ndm
integer :: latgid, longid, nxp, nyp, idp1, jdp1
real :: pi, dtr, sum, mval1, resr
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname

integer, dimension (maxdims)            :: start, count, dimids, ndims

real, dimension (idr)                   :: lonr
real, dimension (jdr)                   :: latr
real, dimension (idrp1)                 :: lonrb
real, dimension (jdrp1)                 :: latrb

character(len=100), dimension (ntilmx)  :: river_input_file

real, allocatable, dimension (:,:)      :: lonrg, latrg, arlatr, lndfr_r, lktau_r
real, allocatable, dimension (:,:,:)    :: cell_a, land_fr, lake_tau, lat, lon

pi= 4.*atan(1.)
dtr= pi/180.
resr= 360./real(idr)

open (10, file= 'out.cr_big_lakes_latlon', form= 'formatted')

read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
close (5)

! ---------------------------------------------------------------------------------------
!  get lon and lat dims from first hydrography file -- should be identical for all files
! ---------------------------------------------------------------------------------------
rcode= NF_OPEN (trim(river_input_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"
    write (6,*) trim(river_input_file(1))
    stop 1
endif

start= 1 ; count= 1
rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 2
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jd)
jdp1= jd+1

write (6,*) 'jd= ', jd

rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 3
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
idp1= id+1

write (6,*) 'id= ', id

rcode= nf_close (ncid)

allocate (lat(id,jd,ntiles), lon(id,jd,ntiles))
allocate (cell_a(id,jd,ntiles), land_fr(id,jd,ntiles), lake_tau(id,jd,ntiles))

! ----------------------------------------------------------------------
!  now get lons and lats from all input hydrography files
! ----------------------------------------------------------------------
do n= 1,ntiles

   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
   write (10,'(i6,3x,a)') n, trim(river_input_file(n))

   rcode= NF_OPEN (trim(river_input_file(n)), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file" ; stop 10
   endif

   start= 1 ; count= 1

   if (ntiles == 1) then
!     regular grid
       rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
       if (rcode /= 0) then
           rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
           if (rcode2 /= 0) then
               write (6,*) "ERROR: cannot find lat variable" ; stop 20
           endif
       endif
       rcode= nf_inq_vardimid (ncid, latid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= jd) then
           write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 25
       endif

       count(1)= jd
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(:,:,n))

       rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
       if (rcode /= 0) then
           rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
           if (rcode2 /= 0) then
               write (6,*) "ERROR: cannot find lon variable" ; stop 30
           endif
       endif
       rcode= nf_inq_vardimid (ncid, lonid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= id) then
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 35
       endif

       count(1)= id
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,:,n))

   else

!     cubic sphere -- assume no edge data
       rcode= nf_inq_varid (ncid, 'y', latid)         ! number of lats
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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(:,:,n))

       rcode= nf_inq_varid (ncid, 'x', lonid)         ! number of lons
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
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, jd ;  stop 35
       endif

       start= 1 ;  count(1)= id ;  count(2)= jd
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,:,n))

   endif

! ----------------------------------------------------------------------
!  now read cell area, land frac, and lake_tau fields
! ----------------------------------------------------------------------

!   write (6,*) 'read cellarea'
   rcode= nf_inq_varid (ncid, 'land_area', varid)
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'cellarea', varid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find land_area/cellarea variable" ; stop 22
       endif
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval1= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval1)
   endif

   cell_a(:,:,n)= mval1
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_a(:,:,n))

!   where (cell_a(:,:,n) == mval1) cell_a(:,:,n)= mval_mdl
   where (cell_a(:,:,n) == mval1) cell_a(:,:,n)= 0.


!   write (6,*) 'read land_frac'
   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval1= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval1)
   endif

   land_fr(:,:,n)= mval1
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_fr(:,:,n))

   where (land_fr(:,:,n) == mval1) land_fr(:,:,n)= mval_mdl


!   write (6,*) 'read lake_tau'
   rcode= nf_inq_varid (ncid, 'lake_tau', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval1= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval1)
   endif

   lake_tau(:,:,n)= mval1
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, lake_tau(:,:,n))

   where (lake_tau(:,:,n) == mval1) lake_tau(:,:,n)= mval_mdl

   rcode= nf_close (ncid)
enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (cell_a(i,j,n) /= mval_mdl) sum= sum + cell_a(i,j,n)
      enddo
   enddo
enddo
write (6,*) "cell_area sum= ", sum/1.e6

! ----------------------------------------------------------------------
! set up regular lats and lons
! ----------------------------------------------------------------------

allocate (latrg(idr,jdr), lonrg(idr,jdr), arlatr(idr,jdr))

do j= 1,jdr
   latr(j)= lat1 + resr*0.5 + real(j-1)*resr
   latrg(:,j)= latr(j)
enddo

latrb(1)= lat1
latrb(jdrp1)= -(lat1)
do j= 2,jdr
   latrb(j)= 0.5*(latr(j)+latr(j-1))
enddo

do i= 1,idr
   lonr(i)= lon1 + resr*0.5 + real(i-1)*resr
   lonrg(i,:)= lonr(i)
enddo

lonrb(1)= lon1
lonrb(idrp1)= lon2
do i= 2,idr
   lonrb(i)= 0.5*(lonr(i)+lonr(i-1))
enddo

! compute area of latitude
arlatr= 0.
do j= 1,jdr
   do i= 1,idr
      arlatr(i,j)= erad*erad*(lonrb(i+1)-lonrb(i))*dtr* &
           (sin(latrb(j+1)*dtr) - sin(latrb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdr
   do i= 1,idr
      sum= sum + arlatr(i,j)
   enddo
enddo
write (6,*) 'sum of regular-grid area= ', sum

write (10,'(/a)') 'regular-grid lats'
do j= 1,jdr
   write (10,'(i6,5f10.3)') j, latr(j), latrb(j), latrb(j+1), arlatr(1,j)
enddo

write (10,'(/a)') 'regular-grid lons'
do i= 1,idr
   write (10,'(i6,3f10.3)') i, lonr(i), lonrb(i), lonrb(i+1)
enddo

deallocate (arlatr)


allocate (lndfr_r(idr,jdr), lktau_r(idr,jdr))

lndfr_r= 0. ;  lktau_r= mval_mdl
do n= 1,ntiles
   do j= 1,jd
      if (mod(j,200) == 0) write (6,*) 'j= ', j, ', tile= ', n
      do i= 1,id
         if (land_fr(i,j,n) < 0.5) go to 170

             do jj= 1,jdr
                if (lat(i,j,n) >= latrb(jj) .and. lat(i,j,n) < latrb(jj+1)) then
                    do ii= 1,idr
                       if (lon(i,j,n) >= lonrb(ii) .and. lon(i,j,n) < lonrb(ii+1)) then
                           lndfr_r(ii,jj)= 1.
                           if (lake_tau(i,j,n) > 0.) then
                               lktau_r(ii,jj)= lake_tau(i,j,n)
                           endif
                           go to 170
                       endif
                    enddo
                endif
             enddo

170      continue
      enddo
   enddo
enddo


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

write (fname, '(a)') 'river_network.nc'
rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lat)
rcode= NF_DEF_DIM (ncid, 'grid_x',  idr,   londim)
rcode= NF_DEF_DIM (ncid, 'grid_y',  jdr,   latdim)

!  create coordinate variables
rcode= NF_DEF_VAR (ncid, 'grid_x',  NF_DOUBLE, 1, londim,  lonid)
rcode= NF_DEF_VAR (ncid, 'grid_y',  NF_DOUBLE, 1, latdim,  latid)

!  create attributes for coordinate variables
!    longitude:
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 16, 'T-cell longitude')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'cartesian_axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 12, 'degrees_east')

!    latitude:
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 13, 'degrees_north')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'cartesian_axis', 1, 'Y')

!    create data variable and attributes
ndims(1)= londim ;  ndims(2)= latdim
rcode= NF_DEF_VAR (ncid, 'land_frac', NF_DOUBLE, 2, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 13, 'land fraction')
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

rcode= NF_DEF_VAR (ncid, 'lake_tau', NF_DOUBLE, 2, ndims, varid2)
rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 8, 'lake_tau')
rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 1, 's')
rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)

rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 12, 'degrees_east')

rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 13, 'degrees_north')

!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data
start= 1 ;  count= 1

count(1)= idr
rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lonr)

count(1)= jdr
rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, latr)

start= 1 ;  count(1)= idr ;  count(2)= jdr
rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lonrg)

start= 1 ;  count(1)= idr ;  count(2)= jdr
rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, latrg)

!    land_frac data
count(1)= idr ;  count(2)= jdr
rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, lndfr_r)

!    lake_tau data
count(1)= idr ;  count(2)= jdr
rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, lktau_r)

!  close netcdf file
rcode= NF_CLOSE (ncid)


deallocate (lat, lon)
deallocate (cell_a, land_fr, lake_tau)
deallocate (latrg, lonrg, lndfr_r, lktau_r)

stop

end



