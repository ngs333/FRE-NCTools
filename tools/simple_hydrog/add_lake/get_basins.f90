program get_basins

! program gets basins upstream of gages (on native grid only)


implicit none

!!TODo:
include 'param_basin.h'

integer, parameter :: maxdims= 3
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: nclmx= 10000
integer, parameter :: len_bname= 4

real, parameter :: mval_mdl= 1.e+20
real, parameter :: erad= 6.371e+6

include 'netcdf.inc'


integer :: k, n, rcode, varid, i, j, ncid, ktr, id, jd, idp1, jdp1
integer :: latid, lonid, idp2, jdp2, rcode2, ndm, l, i1, j1, n1, ip
integer :: jp, np, ngage, nchar, ii, jj, ntiles, m, iwrt, jwrt
integer :: latdim, londim, latgid, longid, varid2
real :: mval1, sum1, pi, dtr, glat, glon, elev
real*4 :: mval_out
character(len=8)   :: huc, ibc
character(len=120) :: fname

integer, dimension (nbasin)           :: igage, jgage, tgage, ncell, &
                                         idx_gage
integer, dimension (maxdims)          :: start, count, dimids, ndims
integer, dimension (nclmx,nbasin)     :: ib_cell, jb_cell, nb_cell

real, dimension (nbasin)              :: obarea, arbas_ntv, arbas_ntv1

character(len=4),  dimension (nbasin) :: bshort
character(len=8),  dimension (nbasin) :: ibas
character(len=40), dimension (nbasin) :: bname_in

real, dimension (ni_cells,nj_cells)    :: out_flow
!!TODO: Moved down to compile
!!    (/  8.,   4.,   2., &
!!      16.,   0.,   1., &
!!       32.,  64., 128. /)

character(len=9) :: area_var = 'land_area'

integer, allocatable, dimension (:)       :: itw, ite, its, itn, itile
integer, allocatable, dimension (:,:)     :: bmask
integer, allocatable, dimension (:,:,:)   :: ipt, jpt
integer, allocatable, dimension (:,:,:,:) :: gage_grd

real, allocatable, dimension (:)          :: ig1, jg1, lat_idx, lon_idx
real, allocatable, dimension (:,:,:)      :: latm, lonm, area, cellarea, &
                                             land_frac, tocell, suba
real*4, allocatable, dimension (:,:)      :: dat_out2

character(len=4),   allocatable, dimension (:)     :: bname
character(len=120), allocatable, dimension (:)     :: hydrography_file

out_flow = reshape ((/  8.,   4.,   2., 16.,   0.,   1.,  32.,  64., 128. /), shape (out_flow))


iwrt= 337 ; jwrt= 229


pi= 4.*atan(1.)
dtr= pi/180.

read (5,*) ntiles
allocate (hydrography_file(ntiles))
do n= 1,ntiles
   read (5,'(a)') hydrography_file(n)
enddo
close (5)

!  read list of basins to be analyzed
open (15, form= 'formatted')
do k= 1,nbasin
   read (15,*) ibas(k), bshort(k), bname_in(k), glat, glon, elev, &
       huc, obarea(k)
enddo
close (15)


! ----------------------------------------------------------------------
!  READ HYDROGRAPHY FILES
! ----------------------------------------------------------------------
!  get lon and lat index dims from first model file
rcode= NF_OPEN (trim(hydrography_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"
    write (6,*) trim(hydrography_file(1))
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

jdp1= jd + 1 ;  jdp2= jd + 2
write (6,*) 'jd= ', jd, ', jdp1= ', jdp1, ', jdp2= ', jdp2

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

idp1= id + 1 ;  idp2= id + 2
write (6,*) 'id= ', id, ', idp1= ', idp1, ', idp2= ', idp2

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)

rcode= nf_close (ncid)

! ----------------------------------------------------------------------
!  now get lons and lats from all input HYDROGRAPHY files
! ----------------------------------------------------------------------
allocate (latm(idp2,jdp2,ntiles), lonm(idp2,jdp2,ntiles))
allocate (cellarea(idp2,jdp2,ntiles), land_frac(idp2,jdp2,ntiles))
allocate (tocell(idp2,jdp2,ntiles), suba(idp2,jdp2,ntiles))

igage= 0 ;  jgage= 0 ;  tgage= 0  ; idx_gage= 0
latm= 0. ; lonm= 0.
do n= 1,ntiles

!   write (6,'(/i6,3x,a)')  n, trim(hydrography_file(n))
   write (10,'(/i6,3x,a)') n, trim(hydrography_file(n))

   rcode= NF_OPEN (trim(hydrography_file(n)), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file" ; stop 10
   endif

   start= 1 ; count= 1

   if (ntiles == 1) then
! regular grid
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
       rcode= nf_get_vara_double (ncid, latid, start, count, latm(2,2:jdp1,n))
       do i= 3,idp1
          latm(i,:,:)= latm(2,:,:)
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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lonm(2:idp1,2,n))
       do j= 3,jdp1
          lonm(:,j,:)= lonm(:,2,:)
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
       start= 1 ;  count= 1 ;  count(1)= id ;  count(2)= jd
       rcode= nf_get_vara_double (ncid, latid, start, count, latm(2:idp1,2:jdp1,n))

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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lonm(2:idp1,2:jdp1,n))

   endif

! ----------------------------------------------------------------------
!  now read cell area, land frac, and tocell fields
! ----------------------------------------------------------------------

! read cell area
   rcode= nf_inq_varid (ncid, 'land_area', varid)
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'cellarea', varid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find land_area/cellarea variable" ; stop 42
       endif
   endif

   call read_2dgrid (ncid, n, varid, mval1, 2, idp1, 2, jdp1, cellarea)
   where (cellarea(:,:,n) == mval1) cellarea(:,:,n)= 0.

! read land_frac
   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find land_frac variable" ; stop 43
   endif

   call read_2dgrid (ncid, n, varid, mval1, 2, idp1, 2, jdp1, land_frac)
   where (land_frac(:,:,n) == mval1) land_frac(:,:,n)= mval_mdl

! read tocell
   rcode= nf_inq_varid (ncid, 'tocell', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find tocell variable" ; stop 44
   endif

   call read_2dgrid (ncid, n, varid, mval1, 2, idp1, 2, jdp1, tocell)
   where (tocell(:,:,n) == mval1) tocell(:,:,n)= mval_mdl

! read suba
   rcode= nf_inq_varid (ncid, 'subA', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find subA variable" ; stop 48
   endif

   call read_2dgrid (ncid, n, varid, mval1, 2, idp1, 2, jdp1, suba)
   where (suba(:,:,n) == mval1) suba(:,:,n)= mval_mdl

!  read gage data
   rcode= nf_inq_varid (ncid, 'bname_gage', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find bname_gage variable" ; stop 45
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), nchar)
   rcode= nf_inq_dimlen (ncid, dimids(2), ngage)
   if (nchar /= len_bname) then
       write (6,*) "ERROR: wrong length of bname text string, ", nchar ; stop 46
   endif

   allocate (bname(ngage), ig1(ngage), jg1(ngage))

   start= 1 ;  count= 1
   do l= 1,ngage
      start(1)= 1 ;  count(1)= nchar
      start(2)= l ;  count(2)= 1
      rcode= nf_get_vara_text (ncid, varid, start, count, bname(l))
   enddo

   rcode= nf_inq_varid (ncid, 'igageloc', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find i_gage variable" ; stop 47
   endif
   start= 1 ;  count= 1 ;  count(1)= ngage
   rcode= nf_get_vara_double (ncid, varid, start, count, ig1(1:ngage))

   rcode= nf_inq_varid (ncid, 'jgageloc', varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find j_gage variable" ; stop 48
   endif
   start= 1 ;  count= 1 ;  count(1)= ngage
   rcode= nf_get_vara_double (ncid, varid, start, count, jg1(1:ngage))

   do l= 1,ngage
      if (trim(bname(l)) == 'dumm') go to 50
      do k= 1,nbasin
         if (trim(bname(l)) == bshort(k)) then
             igage(k)= int(ig1(l))
             jgage(k)= int(jg1(l))
             tgage(k)= n
             idx_gage(k)= k
             go to 50
         endif
      enddo
!      write (6,*) "ERROR: invalid basin name, ", trim(bname(l)) ;  stop 50
50    continue
   enddo

   deallocate (bname, ig1, jg1)

   rcode= nf_close (ncid)
enddo

write (10,*) "GAGE DATA"
do k= 1,nbasin
   write (10,*) bshort(k), igage(k), jgage(k), tgage(k)
enddo

write (6,*)

! okay, now set up the 'halo' for each tile
!   for tocell, land_frac, cellarea, lat, and lon

allocate (itw(ntiles), ite(ntiles), its(ntiles), itn(ntiles))

if (ntiles == 1) then
    itw= 1 ;  ite= 1 ;  its= 1 ;  itn= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_frac)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cellarea)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, latm)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lonm)
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
    enddo

    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_frac)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cellarea)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, latm)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lonm)
endif


! ----------------------------------------------------------------------
!  GET BASIN GRID CELLS USING TOCELL FIELD
! ----------------------------------------------------------------------

allocate (gage_grd(idp2,jdp2,ntiles,nbasin))

! first find all the grid cells containing a gage
gage_grd= -1
do k= 1,nbasin
   if (igage(k) > 0 .and. jgage(k) > 0 .and. tgage(k) > 0) then
       i1= igage(k)+1 ;  j1= jgage(k)+1 ;  n1= tgage(k)
       gage_grd(i1,j1,n1,k)= k
   endif
enddo

ncell= 0 ;  ib_cell= -1 ;  jb_cell= -1 ;  nb_cell= -1
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 66
         i1= i ; j1= j ; n1= n
         ktr= 0

65       continue
         do jj= 1,nj_cells
            jp= j1+jj-2

            do ii= 1,ni_cells
               ip=i1+ii-2

               if (tocell(i1,j1,n1) == out_flow(ii,jj)) then
                   do l= 1,nbasin
                      if (gage_grd(i1,j1,n1,l) > 0) then
                          k= gage_grd(i1,j1,n1,l)
                          ncell(k)= ncell(k) + 1
                          ib_cell(ncell(k),k)= i
                          jb_cell(ncell(k),k)= j
                          nb_cell(ncell(k),k)= n
                      endif
                   enddo
                   if (tocell(i1,j1,n1) == 0.) then
                       go to 66
                   else
                       i1= ip ; j1= jp ;  np= n1
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
                       go to 65
                   endif
               endif

            enddo
         enddo

66       continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop


deallocate (gage_grd)

ktr= 0
do k= 1,nbasin
   write(10,'(a,3i6)') bshort(k), ncell(k)
   do l= 1,ncell(k)
      i= ib_cell(l,k) ; j= jb_cell(l,k) ; n= nb_cell(l,k)
      write (10,'(4i6,2f12.4)') l, i-1, j-1, n, &
         latm(i,j,n), lonm(i,j,n)-360.
      ktr= ktr+1
   enddo
enddo

allocate (area(idp2,jdp2,ntiles))
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         area(i,j,n)= cellarea(i,j,n)*land_frac(i,j,n)
      enddo
   enddo
enddo

arbas_ntv1= 0.
do k= 1,nbasin
   if (ncell(k) == 0) go to 77
       do l= 1,ncell(k)
          i= ib_cell(l,k) ;  j= jb_cell(l,k) ;  n= nb_cell(l,k)
          if (area(i,j,n) /= mval_mdl) arbas_ntv1(k)= arbas_ntv1(k) + area(i,j,n)
       enddo
77 continue
enddo

arbas_ntv= 0.
do k= 1,nbasin
   if (ncell(k) == 0) go to 79
       do l= 1,ncell(k)
          i= ib_cell(l,k) ;  j= jb_cell(l,k) ;  n= nb_cell(l,k)
          if (area(i,j,n) /= mval_mdl) arbas_ntv(k)= arbas_ntv(k) + area(i,j,n)
       enddo
79 continue
enddo


write (10,'(/a)') "AREA DATA:"
do k= 1,nbasin
   write (10,'(a,2f15.2,i6)') bshort(k), obarea(k), arbas_ntv(k)/1.e6, ncell(k)
enddo


! WRITE BASIN-MASK TILE DATA

allocate (bmask(idp2,jdp2))
allocate (dat_out2(idp2,jdp2))

do n= 1,ntiles

   bmask= -1.
   where (land_frac(:,:,n) > 0.) bmask= 0.

   do k= 1,nbasin
      if (ncell(k) > 0) then
          do l= 1,ncell(k)
             if (nb_cell(l,k) == n) go to 105
          enddo
      endif
   enddo
   go to 103
105 continue

   do k= 1,nbasin
      if (ncell(k) > 0) then
          do l= 1,ncell(k)
             if (nb_cell(l,k) == n) then
                 i= ib_cell(l,k) ;  j= jb_cell(l,k)
                 bmask(i,j)= idx_gage(k)
             endif
          enddo
      endif
   enddo

   write (fname, '(a,i1,a)') 'basins_NE.tile', n, '.nc'
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
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 12, 'degrees_east')

!    latitude:
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 15, 'T-cell latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 13, 'degrees_north')

!    create data variable and attributes
   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'basg', NF_INT, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 13, 'basin at gage')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
   rcode= NF_PUT_ATT_INT (ncid, varid, 'missing_value', NF_INT, 1, -1)

   rcode= NF_DEF_VAR (ncid, 'cellarea', NF_FLOAT, 2, ndims, varid2)
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 9, 'cell area')
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 2, 'm2')

   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 12, 'degrees_east')

   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 13, 'degrees_north')

!  leave define mode
   rcode= NF_ENDDEF (ncid)

   count(1)= id
   rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

   count(1)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lonm(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, latm(2:idp1,2:jdp1,n))

!    basin data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_INT (ncid, varid, start, count, bmask(2:idp1,2:jdp1))

!    cell area data
   count(1)= id ;  count(2)= jd
   dat_out2(:,:)= cellarea(:,:,n)
   rcode= NF_PUT_VARA_REAL (ncid, varid2, start, count, dat_out2(2:idp1,2:jdp1))

!  close netcdf file
   rcode= NF_CLOSE (ncid)

103 continue

enddo

deallocate (dat_out2)
deallocate (bmask)



deallocate (lat_idx, lon_idx)
deallocate (itw, ite, its, itn)
deallocate (latm, lonm, cellarea, land_frac, tocell)
deallocate (area)

contains

! **********************************************************************

subroutine create_halo (ntl, id, jd, itw, ite, its, itn, field)

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


! **********************************************************************

subroutine read_2dgrid (ncd, n1, vid, mval_in, i1, i2, j1, j2, gdat)

implicit none

integer, parameter :: maxdims= 3

integer :: rcode, attnum
character(len=40)  :: vunits

integer, intent(in)    :: ncd, n1, vid, i1, i2, j1, j2
real, intent(out)      :: mval_in

real, intent(out)      :: gdat(:,:,:)

integer, dimension (maxdims)   :: start, count

vunits= ' '
rcode= nf_get_att_text (ncd, vid, "units", vunits)
!write (6,*) 'units= ', vunits

mval_in= 1.e+20
rcode= nf_inq_attid (ncd, vid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncd, vid, 'missing_value', mval_in)
endif
!write (6,*) 'mval= ', mval_in

gdat(:,:,n1)= mval_in
start= 1; count= 1; count(1)= i2-i1+1 ; count(2)= j2-j1+1
rcode= nf_get_vara_double (ncd, vid, start, count, gdat(i1:i2,j1:j2,n1))

end subroutine read_2dgrid


! **********************************************************************

end
