program cr_ground_type_files

! create ground type file using reynolds analysis of fao map
! use glcc snow_ice field for ice

implicit none

integer, parameter :: id= 4320, jd= 2160
integer, parameter :: idp1= id+1, jdp1= jd+1
integer, parameter :: id_glcc= 2880, jd_glcc= 1440
integer, parameter :: idp1_glcc= id_glcc+1, jdp1_glcc= jd_glcc+1
integer, parameter :: maxdims= 3
integer, parameter :: ntex= 4
integer, parameter :: igc= 1, igm= 2, igf= 3, igcm= 4, igcf= 5
integer, parameter :: igmf= 6, igcmf= 7, igo= 8, iice= 9, ilake= 10
integer, parameter :: npt= 240
integer, parameter :: nexpnd= 3
integer, parameter :: ngtype= 10, ngsoil= 7

real, parameter :: lat1= -90., lon1= 0., lon2= 360.
real, parameter :: mval= -9999.
real, parameter :: erad= 6.371e+3   ! earth's radius in km
real, parameter :: prox_ratio= 5.

include '/usr/local/include/netcdf.inc'

integer :: i, j, k, n, l, ns, nz, j1, rcode, ncid, londim, latdim, lonbdim
integer :: latbdim, lonid, latid, lonbid, latbid, vecdim, fill_val, attnum
integer :: iout, jout, id_ls, idp1_ls, jd_ls, jdp1_ls, i2, j2, varid
integer :: ii1, ii2, jj1, jj2, ktr, ip, jp, zid, zdim, ktr2
integer :: kmin, kmax

real :: scale, pi, dtr, sum, sum2, mval_ls, r, r1
real :: dlon2, dlat2, gmin, gmax

real*4 :: mval_out

logical :: interp_missing_data, exclude_mountain_glacier, use_antarc_greenland_ice_only

character(len=8)   :: ls_var
character(len=60)  :: fname_nc, long_name, var_units
character(len=100) :: fname, fao_file, lsea_mask_file, glcc_file

integer, dimension (maxdims)         :: start, count, ndims, dimids
integer, dimension (id,jd)           :: itex, gmask, igt, jgt, ls_gmask
integer, dimension (id_glcc,jd_glcc) :: icov
integer, dimension (ngtype)          :: kt

real, dimension (id)                 :: lon, lonin, rlon
real, dimension (idp1)               :: lonb
real, dimension (id_glcc)            :: lonl, loninl
real, dimension (idp1_glcc)          :: lonlb
real, dimension (jd)                 :: lat, latin, rlat
real, dimension (jdp1)               :: latb
real, dimension (jd_glcc)            :: latl, latinl
real, dimension (jdp1_glcc)          :: latlb
real, dimension (ngtype)             :: area_st
real, dimension (id,jd)              :: arlat
real, dimension (id,jd)              :: tex_fr, ice_fao
real, dimension (id,jd,ngtype)       :: gtype
real, dimension (id_glcc,jd_glcc)    :: ice_glcc, arlatl

real, allocatable, dimension (:)       :: lat_ls, lon_ls, latb_ls, lonb_ls
real, allocatable, dimension (:,:)     :: lsea, lsea1

real*4, allocatable, dimension (:)     :: adat1
real*4, allocatable, dimension (:,:,:) :: adat3

real*4, dimension (ngtype)     :: zdat= &
(/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)

character(len=15), dimension (ngtype)  :: covertypename= &
(/ 'Coarse         ', 'Medium         ', 'Fine           ', 'Coarse/Medium  ', &
   'Coarse/Fine    ', 'Medium/Fine    ', 'Coarse/Med/Fine', 'Organic        ', &
   'Ice            ', 'Lake           ' /)

character(len=8)  :: vname_glcc=  'Snow_Ice'
   

pi= 4.*atan(1.)
dtr= pi/180.

open (10, form= 'formatted')

interp_missing_data= .false.
exclude_mountain_glacier= .true.
use_antarc_greenland_ice_only= .false.

read (5,'(a)') fao_file
read (5,'(a)') glcc_file
read (5,'(l1)') exclude_mountain_glacier
read (5,'(l1)') use_antarc_greenland_ice_only
read (5,'(l1)') interp_missing_data
if ( interp_missing_data ) then
    read (5,'(a)') lsea_mask_file
    read (5,'(a)') ls_var
endif

write (6,*) 'exclude_mountain_glacier= ', exclude_mountain_glacier
write (6,*) 'use_antarc_greenland_ice_only= ', use_antarc_greenland_ice_only
write (6,*) 'interp_missing_data= ', interp_missing_data
write (6,*) 'ls_var= ', trim(ls_var)

! ----------------------------------------------------------------------
!  open fao netcdf file; get slope data
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(fao_file), NF_NOWRITE, ncid)
if  (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"  ; stop
endif

start= 1 ;  count= 1

! get latitudes and longitudes of soil grid

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jd) then
    write (6,*) "ERROR: wrong number of lats, fao file, ", n
    stop 1
endif
count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat)

if (lat(1) > -89.) then
    write (6,*) "WARNING: fao lats start at ", lat(1)
    if (lat(1) > 89.) then
        write (6,*) "  -- flip lats to start at 90S -- "
        latin= lat
        do j= 1,jd
           j1= jdp1-j
           lat(j)= latin(j1)
        enddo
    else
        stop 2
    endif
    
endif

latb(1)= lat1
latb(jdp1)= -(lat1)
do j= 2,jd
   latb(j)= 0.5*(lat(j)+lat(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= id) then
    write (6,*) "ERROR: wrong number of lons, fao file, ", n
    stop 3
endif
count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon)

if (lon(1) < -179.) then
    write (6,*) "WARNING: fao lons start at ", lon(1)
    write (6,*) "  -- shift lons to start at zero -- "
    lonin= lon
    do i= 1,id/2
       lon(i)= lonin(i+id/2)
    enddo
    do i= id/2+1,id
       lon(i)= lon2 + lonin(i-id/2)
    enddo
endif

lonb(1)= lon1
lonb(idp1)= lon2
do i= 2,id
   lonb(i)= 0.5*(lon(i)+lon(i-1))
enddo


! compute area of latitude
arlat= 0.
do j= 1,jd
   do i= 1,id
      arlat(i,j)= erad*erad*(lonb(i+1)-lonb(i))*dtr* &
           (sin(latb(j+1)*dtr) - sin(latb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jd
   do i= 1,id
      sum= sum + arlat(i,j)
   enddo
enddo
write (6,*) 'area sum= ', sum

write (10,*) 'lats'
do j= 1,jd
   write (10,'(i6,5f9.3)') j, lat(j), latb(j), latb(j+1), arlat(1,j)
enddo

write (10,*) 'lons'
do i= 1,id
   write (10,'(i6,3f9.3)') i, lon(i), lonb(i), lonb(i+1)
enddo


tex_fr= 0.
write (6,*) 'read fao_txt1'
rcode= nf_inq_varid (ncid, 'fao_txt1', varid)       
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
rcode= nf_inq_dimlen (ncid, dimids(2), l)
if (n /= id .or. l /= jd) then
    write (6,*) "ERROR: inconsistent texture dimensions, id= ", n, ', jd= ', l
    stop 111
endif
start= 1 ;  count= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_int (ncid, varid, start, count, itex)

fill_val= -999999
rcode= nf_inq_attid (ncid, varid, '_FillValue', attnum)
if (rcode == 0) rcode= nf_get_att_int (ncid, varid, '_FillValue', fill_val)
write (6,*) 'fill_val= ', fill_val

long_name=''
rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', long_name)
write (6,*) 'long_name= ', trim(long_name)
   
var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
do j= 1,jd
   j1= jdp1-j
   do i= 1,id/2
      if (itex(i+id/2,j1) /= fill_val) then
          tex_fr(i,j)= real(itex(i+id/2,j1))
      endif
   enddo
   do i= id/2+1,id
      if (itex(i-id/2,j1) /= fill_val) then
          tex_fr(i,j)= real(itex(i-id/2,j1))
      endif
   enddo
enddo
   

! close netcdf file
rcode= nf_close (ncid)

write (6,*) 'fao data read'


gtype= 0.
do j= 1,jd
   do i= 1,id
      if (tex_fr(i,j) /= 0.) then
          if (tex_fr(i,j) == 1.) then
              gtype(i,j,igc)= 1.
          else if (tex_fr(i,j) == 2.) then
              gtype(i,j,igm)= 1.
          else if (tex_fr(i,j) == 3.) then
              gtype(i,j,igf)= 1.
          endif
      endif
   enddo
enddo


area_st= 0. ;  sum= 0. ;  sum2= 0. ;  kt= 0
do j= 1,jd
   do i= 1,id
      do n= 1,ngtype
         if (gtype(i,j,n) /= 0.) then
             sum= sum + arlat(i,j)
             sum2= sum2 + gtype(i,j,n)
             kt(n)= kt(n) + 1
             area_st(n)= area_st(n) + arlat(i,j)
         endif
      enddo
   enddo
enddo
write (6,*)
write (6,*) 'ground-type distribution'
write (6,*) 'sum area= ', sum
write (6,*) 'sum points= ', sum2
write (6,*)

sum2= 0. ; gmin= 9999999 ;  gmax= -9999999 ;  kmin= 0 ;  kmax= 0
do n= 1,ngtype
   write (6,'(a15,f12.4,i10)') covertypename(n), area_st(n)/sum, kt(n)
   sum2= sum2 + area_st(n)/sum
   if (n <= ngsoil) then
       gmin= min(gmin,area_st(n))
       gmax= max(gmax,area_st(n))
       if (gmin == area_st(n)) kmin= n
       if (gmax == area_st(n)) kmax= n
   endif
enddo
write (6,*) 'sum2= ', sum2
write (6,*) 'soil type with least coverage= ', covertypename(kmin)
write (6,*) 'soil type with most coverage=  ', covertypename(kmax)

iout= 269 ;  jout= 1


! ----------------------------------------------------------------------
!  read glcc file -- get lats, lons, snow_ice, waterbod, pwetland
! ----------------------------------------------------------------------
write (6,*)
rcode= NF_OPEN (trim(glcc_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open glcc netcdf file"  
    write (6,*) trim(glcc_file)
    stop 100
endif

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lat variable" ; stop 102
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jd_glcc) then
    write (6,*) "ERROR: inconsistent lat dimension, glcc, jd= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jd_glcc
rcode= nf_get_vara_double (ncid, latid, start, count, latl)

if (latl(1) > -89.) then
    write (6,*) "WARNING: glcc lats start at ", latl(1)
    if (latl(1) > 89.) then
        write (6,*) "  -- flip lats to start at 90S -- "
        latinl= latl
        do j= 1,jd_glcc
           j1= jdp1_glcc-j
           latl(j)= latinl(j1)
        enddo
    else
        stop 22
    endif
endif

latlb(1)= lat1
latlb(jdp1_glcc)= -(lat1)
do j= 2,jd_glcc
   latlb(j)= 0.5*(latl(j)+latl(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= id_glcc) then
    write (6,*) "ERROR: inconsistent lon dimension, glcc, id= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= id_glcc
rcode= nf_get_vara_double (ncid, lonid, start, count, lonl)

if (lonl(1) < -179.) then
    write (6,*) "WARNING: glcc lons start at ", lonl(1)
    write (6,*) "  -- shift lons to start at zero -- "
    loninl= lonl
    do i= 1,id_glcc/2
       lonl(i)= loninl(i+id_glcc/2)
    enddo
    do i= id_glcc/2+1,id_glcc
       lonl(i)= lon2 + loninl(i-id_glcc/2)
    enddo
endif

lonlb(1)= lon1
lonlb(idp1_glcc)= lon2
do i= 2,id_glcc
   lonlb(i)= 0.5*(lonl(i)+lonl(i-1))
enddo

! compute area of latitude
arlatl= 0.
do j= 1,jd_glcc
   do i= 1,id_glcc
      arlatl(i,j)= erad*erad*(lonlb(i+1)-lonlb(i))*dtr* &
           (sin(latlb(j+1)*dtr) - sin(latlb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jd_glcc
   do i= 1,id_glcc
      sum= sum + arlatl(i,j)
   enddo
enddo
write (6,*) 'glcc area sum= ', sum

write (10,'(/a)') 'glcc lats'
do j= 1,jd_glcc
   write (10,'(i6,5f10.3)') j, latl(j), latlb(j), latlb(j+1), arlatl(1,j)
enddo

write (10,'(/a)') 'glcc lons'
do i= 1,id_glcc
   write (10,'(i6,3f10.3)') i, lonl(i), lonlb(i), lonlb(i+1)
enddo


write (6,*)
ice_glcc= mval
write (6,*) 'read glcc data: ', trim(vname_glcc)
rcode= nf_inq_varid (ncid, trim(vname_glcc), varid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc variable ",  trim(vname_glcc)
    stop 110
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), i)
rcode= nf_inq_dimlen (ncid, dimids(2), l)
if (i /= id_glcc .or. l /= jd_glcc) then
    write (6,*) "ERROR: inconsistent Snow_Ice dimensions, glcc, id= ", i, ', jd= ', l
    stop 111
endif
start= 1 ;  count= 1 ;  count(1)= id_glcc ;  count(2)= jd_glcc
rcode= nf_get_vara_int (ncid, varid, start, count, icov)

fill_val= -999999
rcode= nf_inq_attid (ncid, varid, '_FillValue', attnum)
if (rcode == 0) rcode= nf_get_att_int (ncid, varid, '_FillValue', fill_val)
write (6,*) 'fill_val= ', fill_val
   
scale= -99999999.
rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
write (6,*) 'scale= ', scale

rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', long_name)
write (6,*) 'long_name= ', trim(long_name)
   
var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
do j= 1,jd_glcc
   j1= jdp1_glcc-j
!   write (6,*) 'j= ', j, j1
   do i= 1,id_glcc/2
      if (icov(i+id_glcc/2,j1) /= fill_val) then
          ice_glcc(i,j)= real(icov(i+id_glcc/2,j1))*scale
      endif
   enddo
   do i= id_glcc/2+1,id_glcc
      if (icov(i-id_glcc/2,j1) /= fill_val) then
          ice_glcc(i,j)= real(icov(i-id_glcc/2,j1))*scale
      endif
   enddo
enddo

rcode= nf_close (ncid)

do j= 1,jd_glcc
   do i= 1,id_glcc
      if (ice_glcc(i,j) == mval) then
          write (6,*) "WARNING: ice_glcc has missing value, ", ice_glcc(i,j)
      endif
   enddo
enddo


! ----------------------------------------------------------------------
! REMOVE some snow_ice points in mexico and islands
!   lake chapala, e. mex coast, arrecife alacran (scorpion reef), 
!   clipperton island
!   also remove island(s) off Antarctica
! ----------------------------------------------------------------------

write (6,'(a)') '*** RESET snow_ice fraction to zero at following (tropical, semi-tropical) locations ***'
do j= 1,jd_glcc
   if (latl(j) > 10. .and. latl(j) < 28.) then
       do i= 1,id_glcc
          if ((lonl(i) > 250. .and. lonl(i) < 271.) .and. &
               ice_glcc(i,j) > 0.) then
              write (6,'(2i6,3f12.4)') j, i, latl(j), lonl(i), ice_glcc(i,j)
              ice_glcc(i,j)= 0.
          endif
       enddo
   endif
   if (latl(j) > -62. .and. latl(j) < -58.) then
       do i= 1,id_glcc
          if ((lonl(i) > 303. .and. lonl(i) < 320.) .and. &
               ice_glcc(i,j) > 0.) then
              write (6,'(2i6,3f12.4)') j, i, latl(j), lonl(i), ice_glcc(i,j)
              ice_glcc(i,j)= 0.
          endif
       enddo
   endif
   if (latl(j) > -68. .and. latl(j) < -64.) then
       do i= 1,id_glcc
          if ((lonl(i) > 160. .and. lonl(i) < 166.) .and. &
               ice_glcc(i,j) > 0.) then
              write (6,'(2i6,3f12.4)') j, i, latl(j), lonl(i), ice_glcc(i,j)
              ice_glcc(i,j)= 0.
          endif
       enddo
   endif
   if (latl(j) > -70. .and. latl(j) < -66.) then
       do i= 1,id_glcc
          if ((lonl(i) > 268. .and. lonl(i) < 270.) .and. &
               ice_glcc(i,j) > 0.) then
              write (6,'(2i6,3f12.4)') j, i, latl(j), lonl(i), ice_glcc(i,j)
              ice_glcc(i,j)= 0.
          endif
       enddo
   endif
enddo


! ----------------------------------------------------------------------
! interpolate ice data from 1/8-degree to 1/12-degree
! ----------------------------------------------------------------------
!allocate (interp_mask(id_glcc,jd_glcc), interp_out(id,jd))
!interp_mask= 1.
!where (ice_glcc == mval) interp_mask= 0.

!write (6,*) 'here i am 1'
!call flush(6)
!call horiz_interp (ice_glcc, lonlb*dtr, latlb*dtr, lonb*dtr, latb*dtr, &
!     ice_fao, verbose=1, mask_in=interp_mask, mask_out=interp_out, &
!     interp_method="bilinear")
!write (6,*) 'here i am 2'
!call flush(6)
!where (interp_out == 0.) ice_fao= mval
!write (6,*) 'here i am 3'
!call flush(6)

!deallocate (interp_mask, interp_out)

ktr= 0 ;  ktr2= 0
do j= 1,jd_glcc
   do i= 1,id_glcc
      if (ice_glcc(i,j) > 0.) then
          if (use_antarc_greenland_ice_only) then
              if (latl(j) > -60. .and. latl(j) <  59.) go to 17
              if (latl(j) >= 59.) then
                  if (lonl(i) < 286. .or.  lonl(i) > 350.) go to 17
                  if (latl(j) <  73. .and. lonl(i) < 300.) go to 17
                  if (latl(j) <  68. .and. lonl(i) > 335.) go to 17
              endif
              if (latl(j) > 79.3 .and. lonl(i) < 292.0) go to 17
              if (latl(j) > 80.3 .and. lonl(i) < 293.7) go to 17
              if (latl(j) > 81.1 .and. lonl(i) < 296.0) go to 17
              if (latl(j) > 81.6 .and. lonl(i) < 298.4) go to 17
              if (latl(j) > 82.0 .and. lonl(i) < 299.0) go to 17
          endif
          if (exclude_mountain_glacier) then
              if ((latl(j) >  20.  .and. latl(j) <  50. .and. &
                   lonl(i) >   0.  .and. lonl(i) < 120.) .or. &
                  (latl(j) >  35.  .and. latl(j) <  65. .and. &
                   lonl(i) > 228.3 .and. lonl(i) < 290.) .or. &
                   latl(j) > -44.  .and. latl(j) < -40. .and. &
                   lonl(i) > 168.  .and. lonl(i) < 172.) then
                     go to 17
              endif
          endif
          ktr= ktr + 1
          do j2= 1,jd
             if (lat(j2) >= latlb(j) .and. lat(j2) < latlb(j+1)) then
                 do i2= 1,id
                    if (lon(i2) >= lonlb(i) .and. lon(i2) < lonlb(i+1)) then
                        ktr2= ktr2 + 1
!                        write (10,'(2i6,2f12.4,) j2, i2, lat(j2), lon(i2), tex_fr(i2,j2)
                        if (tex_fr(i2,j2) == 0.) gtype(i2,j2,iice)= 1.
                    endif
                 enddo
             endif
          enddo
17        continue
      endif
   enddo
enddo
write (6,*) 'ktr= ', ktr
write (6,*) 'ktr2= ', ktr2

do j= 1,jd
   do i= 1,id
      ktr= 0
      do n= 1,ngtype
         if (gtype(i,j,n) > 0.) ktr= ktr + 1
      enddo
      if (ktr > 1) then
          write (6,*) "ERROR: multiple gtypes defined at cell, ", j, i, lat(j), lon(i)
          stop 50
      endif
   enddo
enddo
      
area_st= 0. ;  sum= 0. ;  sum2= 0. ;  kt= 0
do j= 1,jd
   do i= 1,id
      do n= 1,ngtype
         if (gtype(i,j,n) /= 0.) then
             sum= sum + arlat(i,j)*gtype(i,j,n)
             sum2= sum2 + gtype(i,j,n)
             kt(n)= kt(n) + 1
             area_st(n)= area_st(n) + arlat(i,j)*gtype(i,j,n)
         endif
      enddo
   enddo
enddo
write (6,*)
write (6,*) 'ground-type distribution with snow/ice'
write (6,*) 'sum= ', sum
write (6,*) 'sum2= ', sum2
write (6,*)


gmask= 0 ;  ktr= 0
do j= 1,jd
   do i= 1,id
      do n= 1,ngtype
         if (gtype(i,j,n) /= mval .and. gtype(i,j,n) > 0.) then
             gmask(i,j)= 1
             ktr= ktr + 1
             go to 20
         endif
      enddo
20    continue
   enddo
enddo
write (6,*) 'number of points with valid gtype=   ', ktr


! okay now interpolate to areas of missing data, if so specified
!  get nearest cell to each grid cell of missing data (using slm's
!    "nearest" routine); fill in cells of missing data

if ( interp_missing_data ) then

!  open lsea mask netcdf file
    write (6,*)
    rcode= NF_OPEN (trim(lsea_mask_file), NF_NOWRITE, ncid)
    if  (rcode /= 0) then
        write (6,*) "ERROR: cannot open lsea mask netcdf file"  ; stop
    endif

! get latitudes and longitudes of land-sea grid

    start= 1 ; count= 1
    
    rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
    rcode= nf_inq_vardimid (ncid, latid, dimids)
    rcode= nf_inq_dimlen (ncid, dimids(1), jd_ls)
    
    count(1)= jd_ls
    allocate (lat_ls(jd_ls))
    rcode= nf_get_vara_double (ncid, latid, start, count, lat_ls)

    rcode= nf_inq_varid (ncid, 'latb', latid)         ! number of lats
    rcode= nf_inq_vardimid (ncid, latid, dimids)
    rcode= nf_inq_dimlen (ncid, dimids(1), jdp1_ls)
    
    count(1)= jdp1_ls
    allocate (latb_ls(jdp1_ls))
    rcode= nf_get_vara_double (ncid, latid, start, count, latb_ls)

    rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
    rcode= nf_inq_vardimid (ncid, lonid, dimids)
    rcode= nf_inq_dimlen (ncid, dimids(1), id_ls)

    count(1)= id_ls
    allocate (lon_ls(id_ls))
    rcode= nf_get_vara_double (ncid, lonid, start, count, lon_ls)

    rcode= nf_inq_varid (ncid, 'lonb', lonid)         ! number of lons
    rcode= nf_inq_vardimid (ncid, lonid, dimids)
    rcode= nf_inq_dimlen (ncid, dimids(1), idp1_ls)

    count(1)= idp1_ls
    allocate (lonb_ls(idp1_ls))
    rcode= nf_get_vara_double (ncid, lonid, start, count, lonb_ls)

    write (10,*) 'lats, lsea mask'
    do j= 1,jd_ls
       write (10,'(i6,4f12.3)') j, lat_ls(j), latb_ls(j), latb_ls(j+1)
    enddo

    write (10,*) 'lons, lsea mask'
    do i= 1,id_ls
       write (10,'(i6,3f12.3)') i, lon_ls(i), lonb_ls(i), lonb_ls(i+1)
    enddo

    mval_ls= -999999.
    rcode= nf_inq_varid (ncid, trim(ls_var), varid)
    rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
    if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_ls)
    write (6,*) 'rcode= ', rcode
    write (6,*) 'mval_ls= ', mval_ls

    allocate (lsea(id_ls,jd_ls), lsea1(id_ls,jd_ls))
    
    lsea= mval_ls
    count(1)= id_ls ;  count(2)= jd_ls
    rcode= nf_get_vara_double (ncid, varid, start, count, lsea)

! close netcdf file
    rcode= nf_close (ncid)

    write (6,*) 'lsea mask read'

    where (lsea /= mval_ls) lsea= 1.

! screen out land cell at 36.5W, 54.5S
   do j= 1,jd_ls
      if (lat_ls(j) > -55. .and. lat_ls(j) < -54.) then
          do i= 1,id_ls
             if (lon_ls(i) > 323. .and. lon_ls(i) < 324.) then
                 if (lsea(i,j) == 1.) then
                     write (6,*) "NOTE: changing grid point from land to sea"
                     write (6,*) "lat, lon= ", lat_ls(j), lon_ls(i)
                     lsea(i,j)= mval_ls
                 endif
             endif
          enddo
      endif
   enddo

    write (10,*) 'lsea'
    do j= 1,jd_ls
       write (10,'(360i1)') (int(lsea(i,j)), i= 1,id_ls)
    enddo
    
! check to see if there is a gtype cell with valid data within 10 degrees lat,lon
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea(i,j) == 1.) then
              do j2= 1,jd
                 if (abs(lat(j2)-lat_ls(j)) < 10.) then
                     do i2= 1,id
                        if (abs(lon(i2)-lon_ls(i)) < 10.) then
                            if (gmask(i2,j2) /= 0) go to 31
                        endif
                     enddo
                 endif
              enddo
              write (6,*) 'no gtype cell found w/i 10 degrees of ', i, j, lon_ls(i), lat_ls(j)
              lsea(i,j)= mval_ls
31            continue
          endif
       enddo
    enddo
    
! expand land edges 
    lsea1= lsea
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea1(i,j) /= mval_ls) then
              jj1= j-nexpnd ;  jj2= j+nexpnd
              if (jj1 < 1) jj1= 1
              if (jj2 > jd_ls) jj2= jd_ls
              do j2= jj1,jj2
                 ii1= i-nexpnd ;  ii2= i+nexpnd
                 if (ii1 < 1) ii1= 1
                 if (ii2 > id_ls) ii2= id_ls
                 do i2= ii1,ii2
                    lsea(i2,j2)= 1.
                 enddo
              enddo
          endif
       enddo
    enddo
    
! check to see if there is a gtype cell with valid data within 15 degrees lat,lon
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea(i,j) == 1.) then
              do j2= 1,jd
                 if (abs(lat(j2)-lat_ls(j)) < 15.) then
                     do i2= 1,id
                        if (abs(lon(i2)-lon_ls(i)) < 15.) then
                            if (gmask(i2,j2) /= 0) go to 33
                        endif
                     enddo
                 endif
              enddo
              write (6,*) 'no gtype cell found w/i 15 degrees of ', i, j, lon_ls(i), lat_ls(j)
              stop 33
33            continue
          endif
       enddo
    enddo
    
    ls_gmask= 0 ; ktr= 0
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea(i,j) == mval_ls) go to 76
          do j2= 1,jd
             if (lat(j2) >= latb_ls(j) .and. lat(j2) < latb_ls(j+1)) then
                 do i2= 1,id
                    if (lon(i2) >= lonb_ls(i) .and. lon(i2) < lonb_ls(i+1)) then
                        ls_gmask(i2,j2)= 1
                        ktr= ktr + 1
                    endif
                 enddo
             endif
          enddo
76        continue
       enddo
    enddo
    
    write (6,*) 'ls_gmask created'
    write (6,*) 'ktr= ', ktr

    rlat= lat*dtr
    rlon= lon*dtr

    iout= 2160 ;  jout= 904

    write (10,*) 'iout= ', iout, ', jout= ', jout
    igt= 0
    jgt= 0
    do j= 1,jd
       if (mod(j,144) == 0) write (6,*) 'nearest, j= ', j
       do i= 1,id
          if (ls_gmask(i,j) == 0 .or. gmask(i,j) /= 0) go to 86
          r= (8.*pi)**2    ! initialize r to some big value
          ip= 0 ;jp= 0
          if (i == iout .and. j == jout) write (10,*) 'ls_gmask= ', ls_gmask(iout,jout), ', r= ', r
          if (i == iout .and. j == jout) write (6,*) 'ls_gmask= ', ls_gmask(iout,jout), ', r= ', r

!          do j2= 1,jd
          jj1= j-npt ;  jj2= j+npt
          if (jj1 < 1) jj1= 1
          if (jj2 > jd) jj2= jd
          do j2= jj1,jj2
             dlat2= rlat(j) - rlat(j2)
             
!             do i2= i,id
             ii1= i-npt ;  ii2= i+npt
             if (ii1 < 1) ii1= 1
             if (ii2 > id) ii2= id
             do i2= ii1,ii2
                if (gmask(i2,j2) /= 0) then
                    dlon2= rlon(i) - rlon(i2)
                    do while (dlon2 < -pi)
                       dlon2= dlon2+2.*pi
                    enddo
                    do while (dlon2 > pi)
                       dlon2= dlon2-2.*pi
                    enddo
                    r1= dlon2**2 + (prox_ratio*dlat2)**2
                    if (r1 < r) then
                        ip= i2
                        jp= j2
                        r= r1
                        if (i == iout .and. j == jout) write (10,*) j2, i2, dlat2, dlon2, r1
                    endif
                endif
             enddo
          enddo
          igt(i,j)= ip
          jgt(i,j)= jp
86        continue
       enddo
    enddo

!  now fill missing values
    do j= 1,jd
       do i= 1,id
          if (gmask(i,j) == 0 .and. ls_gmask(i,j) == 1) then
              if (igt(i,j) == 0 .or. jgt(i,j) == 0) then
                  write (6,*) "ERROR: igt or jgt= 0, ", i, j, igt(i,j), jgt(i,j)
                  write (6,*) "   lon= ", lon(i), ", lat= ", lat(j)
                  stop 86
              endif
              do n= 1,ngtype
                 if (gtype(igt(i,j),jgt(i,j),n) /= 0.) then
                     gtype(i,j,n)= gtype(igt(i,j),jgt(i,j),n)
                     go to 88
                 endif
              enddo
              write (6,*) "ERROR: gtype fill not found, ", i, j, igt(i,j), jgt(i,j)
              write (6,*) "   lon= ", lon(i), ", lat= ", lat(j)
              stop 88
88            continue
          endif
       enddo
    enddo
    
    write (6,*) 'for point at ', iout, jout, ', iget, jget= ', igt(iout,jout), jgt(iout,jout)
    do n= 1,ngtype
       write (6,*) n, gtype(igt(iout,jout),jgt(iout,jout),n)
    enddo
    
    deallocate (lat_ls, lon_ls, latb_ls, lonb_ls, lsea, lsea1)
    
    do j= 1,jd
       do i= 1,id
          ktr= 0
          do n= 1,ngtype
             if (gtype(i,j,n) > 0.) ktr= ktr + 1
          enddo
          if (ktr == 1) go to 90
          if (ktr == 0) then
              gtype(i,j,kmax)= 1.
          else if (ktr > 1) then
              write (6,*) "ERROR: multiple gtypes defined at cell, ", j, i, lat(j), lon(i)
              stop 150
          endif
90        continue
       enddo
    enddo
      
endif
    
!  create netcdf file
fname_nc= 'ground_type.nc'
rcode= NF_CREATE (trim(fname_nc), NF_NOCLOBBER, ncid)
rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname_nc), &
                              trim(fname_nc))

!  create dimensions (lon, lonb, lat, latb)
rcode= NF_DEF_DIM (ncid, 'lon',  id,   londim)
rcode= NF_DEF_DIM (ncid, 'lonb', idp1, lonbdim)
rcode= NF_DEF_DIM (ncid, 'lat',  jd,   latdim)
rcode= NF_DEF_DIM (ncid, 'latb', jdp1, latbdim)

rcode= NF_DEF_DIM (ncid, 'zax1_10', ngtype, zdim)

!  create coordinate variables
rcode= NF_DEF_VAR (ncid, 'lon',  NF_FLOAT, 1, londim,  lonid)
rcode= NF_DEF_VAR (ncid, 'lonb', NF_FLOAT, 1, lonbdim, lonbid)
rcode= NF_DEF_VAR (ncid, 'lat',  NF_FLOAT, 1, latdim,  latid)
rcode= NF_DEF_VAR (ncid, 'latb', NF_FLOAT, 1, latbdim, latbid)
rcode= NF_DEF_VAR (ncid, 'zax1_10', NF_FLOAT, 1, zdim, zid)

!    longitude:
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 9, 'longitude')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'cartesian_axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')
rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'edges', 4, 'lonb')

!    longitude edges:
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'long_name', 15, 'longitude edges')
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'cartesian_axis', 1, 'X')
rcode= NF_PUT_ATT_TEXT (ncid, lonbid, 'units', 9, 'degrees_E')

!    latitude:
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'cartesian_axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')
rcode= NF_PUT_ATT_TEXT (ncid, latid, 'edges', 4, 'latb')

!    latitude edges:
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'long_name', 14, 'latitude edges')
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'cartesian_axis', 1, 'Y')
rcode= NF_PUT_ATT_TEXT (ncid, latbid, 'units', 9, 'degrees_N')

!    z axis
rcode= NF_PUT_ATT_TEXT (ncid, zid, 'axis', 1, 'Z')
rcode= NF_PUT_ATT_TEXT (ncid, zid, 'point_spacing', 4, 'even')

!  create data variable and attributes
   
mval_out= mval
ndims(1)= londim
ndims(2)= latdim
ndims(3)= zdim

rcode = NF_DEF_VAR (ncid, 'frac', NF_FLOAT, 3, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 21, 'fraction of soil type')
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
rcode= NF_PUT_ATT_REAL (ncid, varid, 'missing_value', NF_FLOAT, 1, mval_out)

!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data

start= 1 ;  count= 1
      
allocate (adat1(idp1))
count(1)= id
adat1(1:id)= lon(1:id)
rcode= NF_PUT_VARA_REAL (ncid, lonid, start, count, adat1(1:id))

count(1)= idp1
adat1(1:idp1)= lonb(1:idp1)
rcode= NF_PUT_VARA_REAL (ncid, lonbid, start, count, adat1(1:idp1))
deallocate (adat1)

allocate (adat1(jdp1))
count(1)= jd
adat1(1:jd)= lat(1:jd)
rcode= NF_PUT_VARA_REAL (ncid, latid, start, count, adat1(1:jd))

count(1)= jdp1
adat1(1:jdp1)= latb(1:jdp1)
rcode= NF_PUT_VARA_REAL (ncid, latbid, start, count, adat1(1:jdp1))
deallocate (adat1)

count(1)= ngtype
rcode= NF_PUT_VARA_REAL (ncid, zid, start(1), count(1), zdat)

allocate (adat3(id,jd,ngtype))
start= 1 ;  count= 1 ;  count(1)= id ;  count(2)= jd ;  count(3)= ngtype
adat3= gtype
rcode= NF_PUT_VARA_REAL (ncid, varid, start, count, adat3)
deallocate (adat3)

rcode= NF_CLOSE (ncid)



close (10)


stop

end
