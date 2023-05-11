program cr_lake_files_gage

! code is the same as cr_lake_files_new except that gage data is
!   read for each tile, and included in the output netcdf file

implicit none

include "param.h"  !  get ntiles here

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx_all= 1000000
integer, parameter :: nlkmx= 100000
integer, parameter :: nlake_def= 16, nlake_defp1= nlake_def+1
integer, parameter :: iaral= 15, ibalk= 13, ichad= 16, icaspian= nlake_def+1
integer, parameter :: imich= 1, ihuron= 2, ierie= 10, iont= 12
integer, parameter :: id_fine= 2880, jd_fine= 1440
integer, parameter :: idp1_fine= id_fine+1, jdp1_fine= jd_fine+1

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20
real, parameter :: mval_mdl= -9999.
real, parameter :: casp_area= 3.71e11  ! m2
real, parameter :: lake_depth_small= 2.
real, parameter :: ldepth_chad= 2.   ! depth of lake chad in m
real, parameter :: ldepth_max= 50.
real, parameter :: res_fine= 0.0625

include '/usr/local/include/netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, k, rcode2, ndm,  ii, jj, idp2, jdp2, ktr2
integer :: londim, latdim
real :: pi, dtr, sum, mval_frac, mval_tocell, mval_landf
real :: mval_cella, mval_ctn
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname, ctn_latlon_field

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlake_defp1)       :: ktr
integer, dimension (nlkmx_all)         :: il1, jl1, nl1, idx1
integer, dimension (nlkmx,nlake_defp1) :: ilake, jlake, nlake
integer, dimension (id_fine,jd_fine)   :: big_lake

real, dimension (id_fine)              :: lonf
real, dimension (jd_fine)              :: latf
real, dimension (idp1_fine)            :: lonfb
real, dimension (jdp1_fine)            :: latfb
real, dimension (nlkmx_all)            :: fr1
real, dimension (id_fine,jd_fine)      :: arlatf, ctn
real, dimension (nlkmx,nlake_defp1)    :: frlake

!  Lakes (in order) are Michigan, Huron, Superior, Victoria, Tanganyika, Baikal,
!    Great Bear, Malawi, Great Slave, Erie, Winnipeg, Ontario, Balkhash, Ladoga,
!    Aral Sea, Chad, Caspian  ; areas in km2

real, dimension (nlake_defp1)            :: boundw = &
(/ 271.00, 275.56, 267.0, 31.0, 28.0, 102.0, 234.0,  33.0, 240.0, 276.0, &
   260.0,  279.5,  72.0,  28.0, 57.0,  12.0,  45.0 /)

real, dimension (nlake_defp1)            :: bounde = &
(/ 275.56, 280.50, 275.9, 35.0, 31.5, 111.0, 244.0,  36.0, 252.0, 281.5, &
   264.,   285.0,  80.0,  34.0, 63.0,  16.0,  55.0 /)

real, dimension (nlake_defp1)            :: bounds = &
(/   41.0,  42.90,  46.3, -3.0, -9.0,  50.0,  64.0, -16.0,  60.0,  41.0, &
     50.0,  43.0,   44.0, 59.0, 43.0,  12.0,  35.0 /)

real, dimension (nlake_defp1)            :: boundn = &
(/   46.3,  46.65,  49.5,  1.0, -3.0,  58.0,  68.0,  -9.0,  64.0,  43.0, &
     54.0,  44.4,   48.0, 63.0, 48.0,  16.0,  50. /)


real, dimension (nlake_defp1)            :: lake_def_area = &
(/ 58020., 60700., 83270., 68800., 32900., 30500., 31790., 28900., 28440., &
   25680., 24510., 19230., 17400., 18390., 68000., 26000., 3.71e5 /)

!  Lake area from The Water Encyclopedia 2nd ed. (km2)
real, dimension (nlake_defp1)            :: lake_def_area_wenc = &
(/ 58100., 59800., 82680., 69000., 32900., 31500., 30200., 30900., 27200., &
   25700., 24600., 19000., 18200., 17700., 64100., 16600., 374000. /)

!  Lake volume from The Water Encyclopedia 2nd ed. (km3)
real, dimension (nlake_defp1)            :: lake_def_vol_wenc = &
(/  4680.,  3580., 11600.,  2700., 18900., 23000.,  1010.,  7725.,  1070.,  &
     545.,   127.,  1710.,   112.,   908.,  1020.,   44.4,  78200. /)

   character(len=11), dimension (nlake_defp1) :: lake_name = &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ', &
   'Caspian    ' /)
   

character(len=100), dimension (ntiles) :: river_input_file

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:,:)   :: lat, lon, tocell, land_frac, cell_area

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.get_lakes_large', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(a)') ctn_latlon_field
close (5)

write (6,*) 'fort.5 data read'


!  read connected_to_next field that has been fregridded from
!   c360 to 1/8-deg resolution

rcode= NF_OPEN (trim(ctn_latlon_field), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open glcc netcdf file"  
    write (6,*) trim(ctn_latlon_field)
    stop 100
endif

rcode= nf_inq_varid (ncid, 'grid_y', latid)         ! number of lats
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find connected_to_next lat variable" ; stop 102
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jd_fine) then
    write (6,*) "ERROR: inconsistent lat dimension, connected_to_next, jd= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jd_fine
rcode= nf_get_vara_double (ncid, latid, start, count, latf)

latfb(1)= lat1
latfb(jdp1_fine)= lat2
do j= 2,jd_fine
   latfb(j)= 0.5*(latf(j)+latf(j-1))
enddo

rcode= nf_inq_varid (ncid, 'grid_x', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find connected_to_next lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= id_fine) then
    write (6,*) "ERROR: inconsistent lon dimension, connected_to_next, id= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= id_fine
rcode= nf_get_vara_double (ncid, lonid, start, count, lonf)

lonfb(1)= lon1
lonfb(idp1_fine)= lon2
do i= 2,id_fine
   lonfb(i)= 0.5*(lonf(i)+lonf(i-1))
enddo

! compute area of latitude
arlatf= 0.
do j= 1,jd_fine
   do i= 1,id_fine
      arlatf(i,j)= erad*erad*(lonfb(i+1)-lonfb(i))*dtr* &
           (sin(latfb(j+1)*dtr) - sin(latfb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jd_fine
   do i= 1,id_fine
      sum= sum + arlatf(i,j)
   enddo
enddo
write (6,*) 'sum of connected_to_next area= ', sum

write (10,'(/a)') 'lats, connected_to_next'
do j= 1,jd_fine
   write (10,'(i6,5f10.3)') j, latf(j), latfb(j), latfb(j+1), arlatf(1,j)
enddo

write (10,'(/a)') 'lons, connected_to_next'
do i= 1,id_fine
   write (10,'(i6,3f10.3)') i, lonf(i), lonfb(i), lonfb(i+1)
enddo

ctn= mval_mdl
write (6,*) 'read connected_to_next field'
rcode= nf_inq_varid (ncid, 'connected_to_next', varid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find connected_to_next"
    stop 110
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), i)
rcode= nf_inq_dimlen (ncid, dimids(2), l)
if (i /= id_fine .or. l /= jd_fine) then
    write (6,*) "ERROR: inconsistent connected_to_next dimensions, id= ", i, ', jd= ', l
    stop 111
endif
start= 1 ;  count= 1 ;  count(1)= id_fine ;  count(2)= jd_fine
rcode= nf_get_vara_double (ncid, varid, start, count, ctn)

mval_ctn= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_ctn)
write (6,*) 'mval= ', mval_ctn
   
var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
write (6,*) 'units= ', var_units


rcode= nf_close (ncid)



! ----------------------------------------------------------------------
!  get lon and lat dims from first river file -- should be identical
!    for all files
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(river_input_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open river netcdf file"  
    write (6,*) trim(river_input_file(1))
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
write (6,*) 'jd= ', jd
jdp1= jd + 1
jdp2= jd + 2

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
write (6,*) 'id= ', id
idp1= id + 1
idp2= id + 2

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)
  
rcode= nf_close (ncid)



! ----------------------------------------------------------------------
! now open river files -- read lat,lon grids, tocell, land_frac, 
!   cellarea
! ----------------------------------------------------------------------

allocate (lat(id,jd,ntiles), lon(id,jd,ntiles))
allocate (tocell(id,jd,ntiles), land_frac(id,jd,ntiles), cell_area(id,jd,ntiles))

do n= 1,ntiles
   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
   write (10,'(/i6,3x,a)') n, trim(river_input_file(n))

   rcode= NF_OPEN (trim(river_input_file(n)), NF_NOWRITE, ncid)
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
   rcode= nf_get_vara_double (ncid, latid, start, count, lat(:,:,n))
       
       
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
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,:,n))


   write (6,*) 'read tocell'
   rcode= nf_inq_varid (ncid, 'tocell', varid)     ! tocell field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find tocell" ; stop 40
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent tocell dimension, ", ndm, id ;  stop 45
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent tocell dimension, ", ndm, jd ;  stop 45
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(:,:,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_tocell= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_tocell)
   endif
   write (6,*) 'mval= ', mval_tocell
   
   where (tocell(:,:,n) == mval_tocell) tocell(:,:,n)= mval_mdl


   write (6,*) 'read land_frac'
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
   rcode= nf_get_vara_double (ncid, varid, start, count, land_frac(:,:,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_landf= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
   endif
   write (6,*) 'mval= ', mval_landf

   where (land_frac(:,:,n) == mval_landf) land_frac(:,:,n)= mval_mdl


   write (6,*) 'read cellarea'
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
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_area(:,:,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_cella= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_cella)
   endif
   write (6,*) 'mval= ', mval_cella

   where (cell_area(:,:,n) == mval_cella) cell_area(:,:,n)= 0.
   
   rcode= nf_close (ncid)
   
!   write (10,'(/"river lats, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lat(i,j,n), i= 1,id)
!   enddo

!   write (10,'(/"river lons, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lon(i,j,n), i= 1,id)
!   enddo

!   write (10,'(/"cell_area, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f10.1)') (cell_area(i,j,n)/1.e6, i= 1,id)
!   enddo

enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         sum= sum + cell_area(i,j,n)
      enddo
   enddo
enddo
write (6,'(/"sum of cell_area= ",f15.1)') sum/1.e+6

! assign big-lake indices
big_lake= 0
do j= 1,jd_fine
   do i= 1,id_fine
      if (ctn(i,j) > 0.3) then
          do n= 1,nlake_defp1
             if (latf(j) > bounds(n) .and. latf(j) <= boundn(n) .and. &
                 lonf(i) > boundw(n) .and. lonf(i) <= bounde(n)) then
                   big_lake(i,j)= n
             endif
          enddo
      endif
   enddo
enddo

       
! now if ctn lat/lon field has value > 0, assign as lake in c720 field/tile
ktr= 0 ;  ilake= 0 ;  jlake= 0 ;  nlake= 0 ;  ktr2= 0
do j= 1,jd_fine
   if (mod(j,100) == 0) write (6,*) 'j= ', j
   do i= 1,id_fine
      if (ctn(i,j) > 0.3) then
          do n= 1,ntiles
             do jj= 1,jd
                do ii= 1,id
                   if (lat(ii,jj,n) > latfb(j) .and. lat(ii,jj,n) <= latfb(j+1) .and. &
                       lon(ii,jj,n) > lonfb(i) .and. lon(ii,jj,n) <= lonfb(i+1)) then
                           ktr2= ktr2 + 1
                           il1(ktr2)= ii
                           jl1(ktr2)= jj
                           nl1(ktr2)= n
                           idx1(ktr2)= big_lake(i,j)
                           fr1(ktr2)= ctn(i,j)
                           go to 55
                   endif
                enddo
             enddo
          enddo
55        continue
      endif
   enddo
enddo
write (6,*) 'ktr2= ', ktr2

do l= 1,ktr2
   do n= 1,nlake_defp1
      if (idx1(l) == n) then
          ktr(n)= ktr(n) + 1
          ilake(ktr(n),n)= il1(l)
          jlake(ktr(n),n)= jl1(l)
          nlake(ktr(n),n)= nl1(l)
          frlake(ktr(n),n)= fr1(l)
      endif
   enddo
enddo

open (11, file= 'lake_list_16', form= 'formatted')
do l= 1,nlake_def
   write (6,*) l, ', ktr= ', ktr(l), lake_name(l)
   do n= 1,ktr(l)
      if (n > 1) then
          write (11,'(3i5,f6.1,i5)') nlake(n,l), jlake(n,l), ilake(n,l), &
                 frlake(n,l), l
      else
          write (11,'(3i5,f6.1,i5,3x,a)') nlake(n,l), jlake(n,l), ilake(n,l), &
                 frlake(n,l), l, lake_name(l)
      endif
   enddo
enddo
close (11)

open (11, file= 'lake_list_caspian', form= 'formatted')
do l= nlake_defp1,nlake_defp1
   write (6,*) l, ', ktr= ', ktr(l), lake_name(l)
   do n= 1,ktr(l)
      if (n > 1) then
          write (11,'(3i5,f6.1,i5)') nlake(n,l), jlake(n,l), ilake(n,l), &
                 frlake(n,l), l
      else
          write (11,'(3i5,f6.1,i5,3x,a)') nlake(n,l), jlake(n,l), ilake(n,l), &
                 frlake(n,l), l, lake_name(l)
      endif
   enddo
enddo
close (11)


close (10)

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area)

! create netcdf file to check big-lake indices

write (fname, '(a)') 'lake_indices.nc'
rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))
   
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lat)
rcode= NF_DEF_DIM (ncid, 'lon',  id_fine,   londim)
rcode= NF_DEF_DIM (ncid, 'lat',  jd_fine,   latdim)

!  create coordinate variables
rcode= NF_DEF_VAR (ncid, 'lon',  NF_DOUBLE, 1, londim,  lonid)
rcode= NF_DEF_VAR (ncid, 'lat',  NF_DOUBLE, 1, latdim,  latid)

!  create attributes for coordinate variables
!    longitude:
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 9, 'longitude')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')

!    latitude:
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')

!    create data variable and attributes
ndims(1)= londim ;  ndims(2)= latdim
rcode= NF_DEF_VAR (ncid, 'lake_id', NF_INT, 2, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 7, 'lake_id')
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
 
!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data
start= 1 ;  count= 1
      
count(1)= id_fine
rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lonf)

count(1)= jd_fine
rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, latf)
   
start= 1 ;  count(1)= id_fine ;  count(2)= jd_fine
rcode= NF_PUT_VARA_INT (ncid, varid, start, count, big_lake)


!  close netcdf file
rcode= NF_CLOSE (ncid)

   
stop


end



