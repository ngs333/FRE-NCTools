program cp_slope_fields

! =========================================================================
!   program reads lat, lon, cellarea, land_frac, tocell, travel,
!     celllength, min elevation, and river-length fields, and
!     computes the fields: slope_to_next, slope_to_ocean,
!     max_slope_to_next, and max_slope_to_ocean
! =========================================================================

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: ntilmx= 6

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latdim, londim, ii, jj, ip, jp, varid2, varid3
integer :: varid4, varid5, varid6, varid7, i1, j1, k, ktr, rcode2, n1
integer :: np, varid10, ntiles, ndm, latgid, longid, idp2, jdp2, varid8
integer :: varid9, npts_emin
real :: pi, dtr, mval_in, dlat, dlon, csum, clen
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname
logical :: mod_elev_min

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ntilmx)            :: itw, ite, its, itn

real, dimension (ni_cells,nj_cells)    :: out_flow
!!TODO: Initialization moved below to compile.
!!  (/  8.,   4.,   2., &
!!       16.,   0.,   1., &
!!       32.,  64., 128. /)

character(len=100), dimension (ntilmx)  :: river_input_file

integer, allocatable, dimension(:)      :: nmod_emin, imod_emin, jmod_emin
real, allocatable, dimension (:)        :: lat_idx, lon_idx, val_emin
real, allocatable, dimension (:,:,:)    :: cell_a, tocell, land_fr, travel, &
                                           cell_l, rivlen, lat, lon, sin_lat, &
                                           cos_lat, elev_min, slope_n, slope_o, &
                                           slpn_max, slpo_max, elev


out_flow = reshape ((/  8.,   4.,   2., 16.,   0.,   1.,  32.,  64., 128. /), shape (out_flow))

pi= 4.*atan(1.)
dtr= pi/180.

mod_elev_min= .false.

do j= 1,nj_cells
   write (6,'(3f8.0)') (out_flow(i,j), i= 1,ni_cells)
enddo

open (10, file= 'out.cp_slope_fields', form= 'formatted')

read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') mod_elev_min
if (mod_elev_min) then
    read (5,*) npts_emin
    allocate (nmod_emin(npts_emin), imod_emin(npts_emin), jmod_emin(npts_emin))
    allocate (val_emin(npts_emin))
    read (5,*) (nmod_emin(l), l= 1,npts_emin)
    read (5,*) (imod_emin(l), l= 1,npts_emin)
    read (5,*) (jmod_emin(l), l= 1,npts_emin)
    read (5,*) (val_emin(l), l= 1,npts_emin)
endif
close (5)

write (6,*) 'ntiles= ', ntiles
do n= 1,ntiles
   write (6,'(a)') 'input file= ', trim(river_input_file(n))
enddo

! ---------------------------------------------------------------------------
!  get lon and lat dims from first file -- should be identical for all files
! ---------------------------------------------------------------------------
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

allocate (lat_idx(jd))
start= 1 ;  count= 1 ;  count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat_idx)

jdp1= jd + 1
jdp2= jd + 2
write (6,*) 'jd= ', jd, ', jdp1= ', jdp1, ', jdp2= ', jdp2

rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 3
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)

idp1= id + 1
idp2= id + 2
write (6,*) 'id= ', id, ', idp1= ', idp1, ', idp2= ', idp2

rcode= nf_close (ncid)

allocate (lat(idp2,jdp2,ntiles), lon(idp2,jdp2,ntiles))

allocate (cell_a(idp2,jdp2,ntiles), land_fr(idp2,jdp2,ntiles), tocell(idp2,jdp2,ntiles))
allocate (cell_l(idp2,jdp2,ntiles), elev_min(idp2,jdp2,ntiles), travel(idp2,jdp2,ntiles))
allocate (rivlen(idp2,jdp2,ntiles), elev(idp2,jdp2,ntiles))

! ----------------------------------------------------------------------
!  now get lons and lats from all input river files
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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(2,2:jdp1,n))
       do i= 3,idp1
          lat(i,:,:)= lat(2,:,:)
       enddo

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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2,n))
       do j= 3,jdp1
          lon(:,j,:)= lon(:,2,:)
       enddo

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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(2:idp1,2:jdp1,n))

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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2:jdp1,n))

   endif

!   write (10,'(/"river lats, tile", i4)') n
!   do j= 1,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lat(i,j,n), i= 1,idp2)
!   enddo

!   write (10,'(/"river lons, tile", i4)') n
!   do j= 1,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lon(i,j,n), i= 1,idp2)
!   enddo

! ----------------------------------------------------------------------
!  now read cell area, land frac, and tocell fields
! ----------------------------------------------------------------------

   write (6,*) 'read cellarea'
   rcode= nf_inq_varid (ncid, 'land_area', varid)
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'cellarea', varid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find land_area/cellarea variable" ; stop 22
       endif
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   cell_a(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_a(2:idp1,2:jdp1,n))

   where (cell_a(:,:,n) == mval_in) cell_a(:,:,n)= 0.


   write (6,*) 'read land_frac'
   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   land_fr(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_fr(2:idp1,2:jdp1,n))

   where (land_fr(:,:,n) == mval_in) land_fr(:,:,n)= mval_mdl


   write (6,*) 'read tocell'
   rcode= nf_inq_varid (ncid, 'tocell', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   tocell(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(2:idp1,2:jdp1,n))

   where (tocell(:,:,n) == mval_in) tocell(:,:,n)= mval_mdl


   write (6,*) 'read elev_min'
   rcode= nf_inq_varid (ncid, 'elev_min', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   elev_min(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, elev_min(2:idp1,2:jdp1,n))

   where (elev_min(:,:,n) == mval_in) elev_min(:,:,n)= mval_mdl

   write (6,*) 'read elev'
   rcode= nf_inq_varid (ncid, 'elev', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   elev(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, elev(2:idp1,2:jdp1,n))

   where (elev(:,:,n) == mval_in) elev(:,:,n)= mval_mdl

   write (6,*) 'read celllength'
   rcode= nf_inq_varid (ncid, 'celllength', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   cell_l(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_l(2:idp1,2:jdp1,n))

   where (cell_l(:,:,n) == mval_in) cell_l(:,:,n)= mval_mdl


   write (6,*) 'read travel'
   rcode= nf_inq_varid (ncid, 'travel', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   travel(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, travel(2:idp1,2:jdp1,n))

   where (travel(:,:,n) == mval_in) travel(:,:,n)= mval_mdl


   write (6,*) 'read river length'
   rcode= nf_inq_varid (ncid, 'rivlen', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif
   write (6,*) 'mval= ', mval_in

   rivlen(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, rivlen(2:idp1,2:jdp1,n))

   where (rivlen(:,:,n) == mval_in) rivlen(:,:,n)= mval_mdl

   rcode= nf_close (ncid)
enddo

if (mod_elev_min) then
    do l= 1,npts_emin
       write (6,'(3(a,i4))') "MODIFYING ELEV_MIN: tile= ", nmod_emin(l), ", &
                              i= ", imod_emin(l), ", j= ", jmod_emin(l)
       write (6,*) "original elev_min= ", elev_min(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))
       elev_min(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))= val_emin(l)
       write (6,*) "new elev_min= ", elev_min(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))
       write (6,*) 'tocell value= ', tocell(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))
       if (val_emin(l) > elev(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))) then
           write (6,*)
           write (6,*) "*** WARNING: elev_min > elev, ", elev(imod_emin(l)+1,jmod_emin(l)+1,nmod_emin(l))
           write (6,*)
       endif
    enddo
endif


! okay, now set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

if (ntiles == 1) then
    itw(n)= 1 ;  ite(n)= 1 ;  its(n)= 1 ;  itn(n)= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_l)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, elev_min)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, rivlen)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
else
! define tiles to the west, east, south, and north
    do n= 1,ntiles
       if (mod(n,2) == 0) then
           itw(n)= mod(n+5,ntiles) ; ite(n)= mod(n+2,ntiles)
           its(n)= mod(n+4,ntiles) ; itn(n)= mod(n+1,ntiles)
       else
           itw(n)= mod(n+4,ntiles) ; ite(n)= mod(n+1,ntiles)
           its(n)= mod(n+5,ntiles) ; itn(n)= mod(n+2,ntiles)
       endif
       if (itw(n) == 0) itw(n)= ntiles
       if (ite(n) == 0) ite(n)= ntiles
       if (its(n) == 0) its(n)= ntiles
       if (itn(n) == 0) itn(n)= ntiles
       write (6,'("tile ",i4, ", itw= ",i4, ", ite= ",i4, ", its= ",i4, ", itn= ",i4)') &
           n, itw(n), ite(n), its(n), itn(n)
    enddo

    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_l)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, elev_min)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, rivlen)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
endif


!do n= 1,ntiles
!   write (10,'(/"filled tocell, tile", i4)') n
!   do j= 1,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f8.1)') (tocell(i,j,n), i= 1,idp2)
!   enddo
!
!   write (10,'(/"filled lat, tile", i4)') n
!   do j= 1,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lat(i,j,n), i= 1,idp2)
!   enddo
!
!   write (10,'(/"filled lon, tile", i4)') n
!   do j= 1,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lon(i,j,n), i= 1,idp2)
!   enddo
!enddo


allocate (cos_lat(idp2,jdp2,ntiles), sin_lat(idp2,jdp2,ntiles))

! ----------------------------------------------------------------------
! compute sin and cos of latitudes
! ----------------------------------------------------------------------
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         sin_lat(i,j,n)= sin(lat(i,j,n)*dtr)
         cos_lat(i,j,n)= cos(lat(i,j,n)*dtr)
      enddo
   enddo
enddo

allocate (slope_n(idp2,jdp2,ntiles),  slope_o(idp2,jdp2,ntiles))
allocate (slpn_max(idp2,jdp2,ntiles), slpo_max(idp2,jdp2,ntiles))

slope_n= mval_mdl ;  slope_o= mval_mdl
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 165

             do jj= 1,nj_cells
                jp= j+jj-2

                do ii= 1,ni_cells
                   ip=i+ii-2

                    if (tocell(i,j,n) == out_flow(ii,jj)) then

                        if (ip /= i .or. jp /= j) then
                            dlat= (lat(i,j,n)-lat(ip,jp,n))*dtr
                            dlon= (lon(i,j,n)-lon(ip,jp,n))*dtr
                            clen= 2.*asin( sqrt( (sin(dlat*0.5))**2. + &
                               cos_lat(i,j,n)*cos_lat(ip,jp,n)*(sin(dlon*0.5))**2. ) )* &
                               erad*1000.
                        else
                            clen= 0.
                        endif
                        if (abs(clen-cell_l(i,j,n)) > 1.e-8) then
                            write (6,*) "ERROR: inconsistent celllength"
                            write (6,*) j, i, n, clen, cell_l(i,j,n)
                            stop 165
                        endif

                        if (elev_min(i,j,n) /= mval_mdl .and. cell_l(i,j,n) /= mval_mdl) then
                            if (travel(i,j,n) > 1. .and. cell_l(i,j,n) /= 0. .and. &
                                elev_min(ip,jp,n) /= mval_mdl) then
                                slope_n(i,j,n)= (elev_min(i,j,n)-elev_min(ip,jp,n))/cell_l(i,j,n)
                            else if (travel(i,j,n) == 1. .and. cell_l(i,j,n) /= 0.) then
                                slope_n(i,j,n)= elev_min(i,j,n)/cell_l(i,j,n)
!                            else if (travel(i,j,n) == 0.) then
!                                slope_n(i,j,n)= 0.
                            endif
                        endif

                        if (elev_min(i,j,n) /= mval_mdl .and. rivlen(i,j,n) /= mval_mdl) then
                            if (travel(i,j,n) > 0.  .and. rivlen(i,j,n) /= 0.) then
                                slope_o(i,j,n)= elev_min(i,j,n)/rivlen(i,j,n)
 !                           else if (travel(i,j,n) == 0.) then
 !                               slope_o(i,j,n)= 0.
                            endif
                        endif

                    endif
                enddo
             enddo

165      continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop
write (6,*) 'slope variables computed'


! ----------------------------------------------------------------------
! follow river downstream; max_slope_to_next and max_slope_to_ocean
! ----------------------------------------------------------------------

! first find all the grid cells where tocell=0
ktr= 0
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 60
         if (tocell(i,j,n) == 0.) then
             ktr= ktr + 1
         endif
60       continue
      enddo
   enddo
enddo
write (6,*) 'number of grid cells where tocell is zero= ', ktr

!write (6,*) 'tocell= ', tocell(109,10), cell_a(109,10)

slpn_max= -999999. ;  slpo_max= -999999.
do n= 1,ntiles
   do j= 2,jdp1
      if (mod(j,200) == 0) write (6,*) 'j= ', j, ', tile= ', n
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 170
!             write (6,*) n, j, i, tocell(i,j,n)
             i1= i ; j1= j ; n1= n
             ktr= 0 ;  csum= 0.

             if (cell_a(i,j,n) == mval_mdl) then
                 write (6,'(a,2i5,2f10.3,2f10.0)') 'cell_a is missing value, ', j, i, &
                     lat(j,j,n), lon(i,j,n), tocell(i,j,n), cell_a(i,j,n)
                 stop 120
             endif

120          continue
             do jj= 1,nj_cells
                jp= j1+jj-2

                do ii= 1,ni_cells
                   ip=i1+ii-2

                if (tocell(i1,j1,n1) == out_flow(ii,jj)) then
                    ktr= ktr + 1

                    if (ip /= i1 .or. jp /= j1) then
                        dlat= (lat(i1,j1,n)-lat(ip,jp,n))*dtr
                        dlon= (lon(i1,j1,n)-lon(ip,jp,n))*dtr
                        clen=      2.*asin( sqrt( (sin(dlat*0.5))**2. + &
                           cos_lat(i1,j1,n)*cos_lat(ip,jp,n)*(sin(dlon*0.5))**2. ) )* &
                           erad*1000.
                    else
                        clen= 0.
                    endif
                    if (ktr == 1) then
                        if (abs(clen-cell_l(i,j,n)) > 1.e-8) then
                            write (6,*) "ERROR: inconsistent celllength"
                            write (6,*) j, i, n, clen, cell_l(i,j,n)
                            stop 265
                        endif
                    endif
                    csum= csum + clen

                    slpn_max(i,j,n)= max(slpn_max(i,j,n),slope_n(i1,j1,n1))
                    slpo_max(i,j,n)= max(slpo_max(i,j,n),slope_o(i1,j1,n1))

                    if (tocell(i1,j1,n1) == 0.) then
                        if (real(ktr-1) /= travel(i,j,n)) then
                            write (6,*) "ERROR: inconsistent travel"
                            write (6,*) j, i, n, ktr-1, travel(i,j,n)
                            stop 275
                        endif
!                        if (abs((csum-clen)-rivlen(i,j,n)) > 1.e-18) then  ! 28 feb 2022, change this due to ifort change
                        if (abs((csum-clen)-rivlen(i,j,n)) > 1.e-8) then
                            write (6,*) "ERROR: inconsistent rivlen"
                            write (6,'(3i8,2f30.18)') j, i, n, csum-clen, rivlen(i,j,n)
                            stop 285
                        endif
                        go to 170
                    else
                        i1= ip ; j1= jp ;  np= n1
                        if ( (ip == idp2 .and. (jp == 1 .or. jp == jdp2)) .or. &
                             (ip == 1    .and. (jp == 1 .or. jp == jdp2)) .or. &
                             (jp == 1    .and. (ip == 1 .or. ip == idp2)) .or. &
                             (jp == jdp2 .and. (ip == 1 .or. ip == idp2))) then
                               write (6,*) "WARNING: i and j edge, ", i, j, n, ip, jp
                        endif
                        if (ntiles == 1) then
                            if (ip == idp2) then
                                i1= 2
                            else if (ip == 1) then
                                i1= idp1
                            endif
                        else
                            if (mod(n1,2) == 0) then
                                if (ip == idp2) then
                                    i1= idp1-jp+2 ; j1= 2          ; n1= ite(np)
                                else if (ip == 1) then
                                    i1= idp1      ; j1= jp         ; n1= itw(np)
                                endif
                                if (jp == jdp2) then
                                    i1= ip        ; j1= 2          ; n1= itn(np)
                                else if (jp == 1) then
                                    i1= idp1      ; j1= jdp1-ip+2  ; n1= its(np)
                                endif
                            else
                                if (ip == idp2) then
                                    i1= 2         ; j1= jp         ; n1= ite(np)
                                else if (ip == 1) then
                                    i1= idp1-jp+2 ; j1= jdp1       ; n1= itw(np)
                                endif
                                if (jp == jdp2) then
                                    i1= 2         ; j1= jdp1-ip+2  ; n1= itn(np)
                                else if (jp == 1) then
                                    i1= ip        ; j1= jdp1       ; n1= its(np)
                                endif
                            endif
                        endif
                        go to 120
                    endif

                else if (tocell(i1,j1,n1) == mval_mdl) then

          ! tocell undefined, should not occur
                    write (6,'(a,3i5,2f10.3,3i5,2f10.3,f15.1,f7.0)') "WARNING: tocell has missing value, ", &
                       j, i, n, lat(i,j,n), lon(i,j,n), j1, i1, n1, lat(i1,j1,n1), lon(i1,j1,n1), &
                       cell_a(i1,j1,n1), land_fr(i1,j1,n1)
                    stop 170
                endif

                enddo
             enddo

170      continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop
write (6,*) 'network variables computed'

where (tocell == mval_mdl)
   slope_n= mval_mdl
   slope_o= mval_mdl
   slpn_max= mval_mdl
   slpo_max= mval_mdl
endwhere


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'river_network.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lat)
   rcode= NF_DEF_DIM (ncid, 'grid_x',  id,   londim)
   rcode= NF_DEF_DIM (ncid, 'grid_y',  jd,   latdim)

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
   rcode= NF_DEF_VAR (ncid, 'tocell', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 28, 'direction to downstream cell')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'travel', NF_DOUBLE, 2, ndims, varid4)
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 42, &
             'cells left to travel before reaching ocean')
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'celllength', NF_DOUBLE, 2, ndims, varid5)
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', 11, 'cell length')
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'land_frac', NF_DOUBLE, 2, ndims, varid6)
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', 13, 'land fraction')
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 12, 'degrees_east')

   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 13, 'degrees_north')

   rcode= NF_DEF_VAR (ncid, 'rivlen', NF_DOUBLE, 2, ndims, varid2)
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 12, 'river length')
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'elev_min', NF_DOUBLE, 2, ndims, varid8)
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'long_name', 17, 'minimum elevation')
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid8, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'slope_to_next', NF_DOUBLE, 2, ndims, varid7)
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'long_name', 37, 'change in grid cell minimum elevation')
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid7, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'slope_to_ocean', NF_DOUBLE, 2, ndims, varid9)
   rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'long_name', 33, 'change in river minimum elevation')
   rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid9, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'max_slope_to_next', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 26, 'max value of slope_to_next')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'max_slope_to_ocean', NF_DOUBLE, 2, ndims, varid10)
   rcode= NF_PUT_ATT_TEXT (ncid, varid10, 'long_name', 27, 'max value of slope_to_ocean')
   rcode= NF_PUT_ATT_TEXT (ncid, varid10, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid10, 'missing_value', NF_DOUBLE, 1, mval_mdl)

!  leave define mode
   rcode= NF_ENDDEF (ncid)

!  write coordinate data
   start= 1 ;  count= 1

   count(1)= id
   rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

   count(1)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lon(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, lat(2:idp1,2:jdp1,n))

!    tocell data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, tocell(2:idp1,2:jdp1,n))

!    travel data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, travel(2:idp1,2:jdp1,n))

!    cell length data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, cell_l(2:idp1,2:jdp1,n))

!    land fraction data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, land_fr(2:idp1,2:jdp1,n))

!    min elevation data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid8, start, count, elev_min(2:idp1,2:jdp1,n))

!    river length data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, rivlen(2:idp1,2:jdp1,n))

!    slope_to_next data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, slope_n(2:idp1,2:jdp1,n))

!    slope_to_ocean data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid9, start, count, slope_o(2:idp1,2:jdp1,n))

!    max_slope_to_next data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, slpn_max(2:idp1,2:jdp1,n))

!    max_slope_to_ocean data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid10, start, count, slpo_max(2:idp1,2:jdp1,n))

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, sin_lat, cos_lat)
deallocate (cell_a, tocell, land_fr, rivlen)
deallocate (travel, cell_l, elev_min, elev)
deallocate (slope_n, slope_o, slpn_max, slpo_max)

contains


! ----------------------------------------------------------------------
subroutine create_halo (ntl, id, jd, itw, ite, its, itn, field)
! ----------------------------------------------------------------------

implicit none

integer :: n, i, j, ip1, ip2, jp1, jp2

integer, intent(in)     :: ntl, id, jd
integer, intent(in)     :: itw(:), ite(:), its(:), itn(:)

real, intent(inout)     :: field(:,:,:)

ip1= id + 1  ;  jp1= jd + 1
ip2= id + 2  ;  jp2= jd + 2


if (ntl == 1) then
    field(1,:,ntl)= field(ip1,:,ntl)
    field(ip2,:,ntl)= field(2,:,ntl)
else
    do n= 1,ntl
! define tiles to the west, east, south, and north
       if (mod(n,2) == 0) then
           do j= 2,jp1
              field(1,j,n)=     field(ip1,j,itw(n))  ! western edge
           enddo

           do j= 2,jp1
              i= ip1-j+2
              field(ip2,j,n)=  field(i,2,ite(n))     ! eastern edge
           enddo

           do i= 2,ip1
              j= jp1-i+2
              field(i,1,n)=     field(ip1,j,its(n))  ! southern edge
           enddo

           do i= 2,ip1
              field(i,jp2,n)=  field(i,2,itn(n))     ! northern edge
           enddo
       else
           do j= 2,jp1
              i= ip1-j+2
              field(1,j,n)=     field(i,jp1,itw(n))  ! western edge
           enddo

           do j= 2,jp1
              field(ip2,j,n)=  field(2,j,ite(n))     ! eastern edge
           enddo

           do i= 2,ip1
              field(i,1,n)=     field(i,jp1,its(n))  ! southern edge
           enddo

           do i= 2,ip1
              j= jp1-i+2
              field(i,jp2,n)=  field(2,j,itn(n))     ! northern edge
           enddo
       endif
    enddo

endif

end subroutine create_halo

end



