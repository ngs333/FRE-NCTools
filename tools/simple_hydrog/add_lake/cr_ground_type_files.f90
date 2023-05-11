program cr_ground_type_files

! test groundwater parameterization computations

implicit none

integer, parameter :: id= 2880, jd= 1440
integer, parameter :: idp1= id+1, jdp1= jd+1
integer, parameter :: maxdims= 3
integer, parameter :: nslptex= 10
integer, parameter :: ntex= 4
integer, parameter :: icoarse= 1, imed= 2, ifine= 3, iundef= 4
integer, parameter :: igc= 1, igm= 2, igf= 3, igcm= 4, igcf= 5
integer, parameter :: igmf= 6, igcmf= 7, igo= 8, iice= 9, ilake= 10
integer, parameter :: iwbd= 1, ipwt= 2, isi= 3
integer, parameter :: npt= 160
integer, parameter :: nexpnd= 6
integer, parameter :: ngtype= 10, ngsoil= 7
integer, parameter :: nvar_glcc= 3

real, parameter :: lat1= -90., lon1= 0., lon2= 360.
real, parameter :: res= 0.5
real, parameter :: mval= -9999.
real, parameter :: erad= 6.371e+3   ! earth's radius in km
real, parameter :: prox_ratio= 5.

include '/usr/local/include/netcdf.inc'

integer :: i, j, k, n, l, ns, nz, j1, rcode, ncid, londim, latdim, lonbdim
integer :: latbdim, lonid, latid, lonbid, latbid, vecdim, fill_val, attnum
integer :: ivar, iout, jout, id_ls, idp1_ls, jd_ls, jdp1_ls, i2, j2, varid
integer :: ii1, ii2, jj1, jj2, ktr, ip, jp, zid, zdim, ics, imd, ifn, ktr2
integer :: ichk, ktr3, ktr4, kmin, kmax

real :: scale, pi, dtr, sum, sum2, h1, mval_ls, r, r1
real :: dlon2, dlat2, gmin, gmax

real*4 :: mval_out

logical :: interp_missing_data

character(len=8)   :: ls_var
character(len=60)  :: fname_nc, long_name, var_units
character(len=100) :: fname, fao_slope_file, lsea_mask_file, glcc_file

integer, dimension (maxdims)         :: start, count, ndims, dimids
integer, dimension (id,jd)           :: itex, gmask, igt, jgt, icov, ls_gmask
integer, dimension (nslptex)         :: kt

real, dimension (id)                 :: lon, lonin, rlon, lonl
real, dimension (idp1)               :: lonb, lonlb
real, dimension (jd)                 :: lat, latin, rlat, latl
real, dimension (jdp1)               :: latb, latlb
real, dimension (nslptex)            :: area_st
real, dimension (id,jd)              :: arlat, area_ij, arlatl
real, dimension (id,jd,ntex)         :: tex_fr
real, dimension (id,jd,ngtype)       :: gtype
real, dimension (id,jd,nvar_glcc)    :: wbdat
real, dimension (id,jd,nslptex)      :: slptex

real, allocatable, dimension (:)       :: lat_ls, lon_ls, latb_ls, lonb_ls
real, allocatable, dimension (:,:)     :: lsea, lsea1

real*4, allocatable, dimension (:)     :: adat1
real*4, allocatable, dimension (:,:,:) :: adat3

real*4, dimension (ngtype)     :: zdat= &
(/ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. /)

character(len=11), dimension (nslptex)  :: slopetexname= &
(/ 'CoarseFlat ', 'CoarseHilly', 'CoarseSteep', 'MedFlat    ', 'MedHilly   ', &
   'MedSteep   ', 'FineFlat   ', 'FineHilly  ', 'FineSteep  ', 'Undefined  ' /)

character(len=15), dimension (ngtype)  :: covertypename= &
(/ 'Coarse         ', 'Medium         ', 'Fine           ', 'Coarse/Medium  ', &
   'Coarse/Fine    ', 'Medium/Fine    ', 'Coarse/Med/Fine', 'Organic        ', &
   'Ice            ', 'Lake           ' /)

character(len=8), dimension (nvar_glcc) :: vname_glcc= &
(/ 'WaterBod', 'PWetland', 'Snow_Ice' /)
   

pi= 4.*atan(1.)
dtr= pi/180.

open (10, form= 'formatted')


read (5,'(a)') fao_slope_file
read (5,'(a)') glcc_file
read (5,'(l1)') interp_missing_data
if ( interp_missing_data ) then
    read (5,'(a)') lsea_mask_file
    read (5,'(a)') ls_var
endif

write (6,*) 'interp_missing_data= ', interp_missing_data
write (6,*) 'ls_var= ', trim(ls_var)

! ----------------------------------------------------------------------
!  open fao netcdf file; get slope data
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(fao_slope_file), NF_NOWRITE, ncid)
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
write (6,*) 'sum= ', sum

write (10,*) 'lats'
do j= 1,jd
   write (10,'(i6,5f9.3)') j, lat(j), latb(j), latb(j+1), arlat(1,j)
enddo

write (10,*) 'lons'
do i= 1,id
   write (10,'(i6,3f9.3)') i, lon(i), lonb(i), lonb(i+1)
enddo


tex_fr= 0. ;  slptex= fill_val
do ns= 1,nslptex

   write (6,*) 'read ', trim(slopetexname(ns))
   rcode= nf_inq_varid (ncid, trim(slopetexname(ns)), varid)       
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), n)
   rcode= nf_inq_dimlen (ncid, dimids(2), l)
   if (n /= id .or. l /= jd) then
       write (6,*) "ERROR: inconsistent tex/slope dimensions, id= ", n, ', jd= ', l
       stop 111
   endif
   start= 1 ;  count= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_int (ncid, varid, start, count, itex)

   fill_val= -999999
   rcode= nf_inq_attid (ncid, varid, '_FillValue', attnum)
   if (rcode == 0) rcode= nf_get_att_int (ncid, varid, '_FillValue', fill_val)
   write (6,*) 'fill_val= ', fill_val

   scale= -99999999.
   rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
   if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
   write (6,*) 'scale= ', scale
!   scale= 0.00111111113801599

   long_name=''
   rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
   if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', long_name)
   write (6,*) 'long_name= ', trim(long_name)
   
   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   i= len_trim(slopetexname(ns))
   write (6,*) 'i= ', i
   if (slopetexname(ns)(1:6) == 'Coarse') then
       write (6,*) 'Texture=Coarse'
       ivar= icoarse
   else if (slopetexname(ns)(1:3) == 'Med') then
       write (6,*) 'Texture=Med'
       ivar= imed
   else if (slopetexname(ns)(1:4) == 'Fine') then
       write (6,*) 'Texture=Fine'
       ivar= ifine
   else if (slopetexname(ns)(1:i) == 'Undefined') then
       write (6,*) 'Texture=Undefined'
       ivar= iundef
   else
       write (6,*) 'ERROR: could not determine texture type, ', trim(slopetexname(ns))
       stop 1
   endif

! flip lats and lons, scale by scale factor
   do j= 1,jd
      j1= jdp1-j
      do i= 1,id/2
         if (itex(i+id/2,j1) /= fill_val) then
             tex_fr(i,j,ivar)= tex_fr(i,j,ivar) + real(itex(i+id/2,j1))*scale
             slptex(i,j,ns)= real(itex(i+id/2,j1))*scale
         endif
      enddo
      do i= id/2+1,id
         if (itex(i-id/2,j1) /= fill_val) then
             tex_fr(i,j,ivar)= tex_fr(i,j,ivar) + real(itex(i-id/2,j1))*scale
             slptex(i,j,ns)= real(itex(i-id/2,j1))*scale
         endif
      enddo
   enddo
   
enddo

area_st= 0. ;  sum= 0. ;  sum2= 0.
do j= 1,jd
   do i= 1,id
      do ns= 1,nslptex
         if (slptex(i,j,ns) /= fill_val) then
             sum= sum + arlat(i,j)*slptex(i,j,ns)
             sum2= sum2 + slptex(i,j,ns)
             area_st(ns)= area_st(ns) + arlat(i,j)*slptex(i,j,ns)
         else
             write (6,*) 'slptex= fill_val, ', i, j
             stop 100
         endif
      enddo
   enddo
enddo
write (6,*) 'sum= ', sum
write (6,*) 'sum2= ', sum2
write (6,*)

write (6,'(a9,3a12)')   'Texture', 'Flat ', 'Hilly', 'Steep'
write (6,'(a9,3f12.1)') 'Coarse ', (area_st(l), l= 1,3)
write (6,'(a9,3f12.1)') 'Medium ', (area_st(l), l= 4,6)
write (6,'(a9,3f12.1)') 'Fine   ', (area_st(l), l= 7,9)
write (6,*)
write (6,'(a9,f12.1)') 'Undefined', area_st(nslptex)

write (6,*)
write (6,'(3a15)') (trim(slopetexname(l)), l= 1,3)
write (6,'(3a15)') (trim(slopetexname(l)), l= 4,6)
write (6,'(3a15)') (trim(slopetexname(l)), l= 7,9)
write (6,'(3a15)') (trim(slopetexname(l)), l= 10,10)
write (6,*)

write (6,'(a9,3a9)')   'Texture', 'Flat ', 'Hilly', 'Steep'
write (6,'(a9,3f9.3)') 'Coarse ', (area_st(l)/sum, l= 1,3)
write (6,'(a9,3f9.3)') 'Medium ', (area_st(l)/sum, l= 4,6)
write (6,'(a9,3f9.3)') 'Fine   ', (area_st(l)/sum, l= 7,9)
write (6,*)
write (6,'(a9,f9.3)') 'Undefined', area_st(nslptex)/sum

area_st= 0. ;  sum= 0. ;  sum2= 0. ;  ktr= 0
do j= 1,jd
   do i= 1,id
      ichk= 0
      do ns= 1,ntex
         if (tex_fr(i,j,ns) /= fill_val) then
             if (ichk == 0) then
                 ktr= ktr + 1
                 ichk= 1
             endif
             sum= sum + arlat(i,j)*tex_fr(i,j,ns)
             sum2= sum2 + tex_fr(i,j,ns)
             area_st(ns)= area_st(ns) + arlat(i,j)*tex_fr(i,j,ns)
         endif
      enddo
   enddo
enddo
write (6,*)
write (6,*) 'sum= ', sum
write (6,*) 'sum2= ', sum2
write (6,*) 'ktr= ', ktr
write (6,*)

write (6,'(a9,4a12)')   'Slope ', 'Coarse', 'Med', 'Fine', 'Undef'
write (6,'(a9,4f12.1)') 'All    ', (area_st(l), l= 1,4)
write (6,'(a9,4f12.3)') 'All    ', (area_st(l)/sum, l= 1,4)
write (6,*)

write (10,*) 'slope data'
write (6,*) 'lat(720)= ', lat(720)
do j= 1,jd
   if (j == 720) then
       do i= 1,id
          write (10,'(i6,2f20.7)') i, lon(i), tex_fr(i,j,1)+ tex_fr(i,j,2)+ tex_fr(i,j,3)+ tex_fr(i,j,4)
       enddo
   endif
enddo

! assign undefined to medium texture category

where (tex_fr(:,:,2) /= fill_val .and. tex_fr(:,:,ntex) /= fill_val)
    tex_fr(:,:,2)= tex_fr(:,:,2) + tex_fr(:,:,ntex)
endwhere

area_st= 0. ;  sum= 0. ;  sum2= 0.
do j= 1,jd
   do i= 1,id
      do ns= 1,ntex-1
         if (tex_fr(i,j,ns) /= fill_val) then
             sum= sum + arlat(i,j)*tex_fr(i,j,ns)
             sum2= sum2 + tex_fr(i,j,ns)
             area_st(ns)= area_st(ns) + arlat(i,j)*tex_fr(i,j,ns)
         endif
      enddo
   enddo
enddo
write (6,*)
write (6,*) 'assign undefined to flat'
write (6,*) 'sum= ', sum
write (6,*) 'sum2= ', sum2
write (6,*)

write (6,'(a9,4a12)')   'Slope ', 'Coarse', 'Med', 'Fine', 'Undef'
write (6,'(a9,4f12.1)') 'All    ', (area_st(l), l= 1,4)
write (6,'(a9,4f12.3)') 'All    ', (area_st(l)/sum, l= 1,4)
write (6,*)

! close netcdf file
rcode= nf_close (ncid)

write (6,*) 'fao data read'


gtype= 0.
do j= 1,jd
   do i= 1,id
      if (tex_fr(i,j,icoarse) /= fill_val .or. tex_fr(i,j,imed) /= fill_val &
            .or. tex_fr(i,j,ifine) /= fill_val) then
          ics= 0 ;  imd= 0 ;  ifn= 0
          if (tex_fr(i,j,icoarse) /= fill_val .and. tex_fr(i,j,icoarse) > 0.) ics= 1
          if (tex_fr(i,j,imed) /= fill_val    .and. tex_fr(i,j,imed) > 0.) imd= 1
          if (tex_fr(i,j,ifine) /= fill_val   .and. tex_fr(i,j,ifine) > 0.) ifn= 1
          
          if (ics == 1 .and. imd+ifn == 0) then
              gtype(i,j,igc)= tex_fr(i,j,1)
          else if (ics == 1 .and. imd+ifn == 2) then
              gtype(i,j,igcmf)= tex_fr(i,j,icoarse)+tex_fr(i,j,imed)+tex_fr(i,j,ifine)
          else if (ics == 1 .and. imd == 1) then
              gtype(i,j,igcm)= tex_fr(i,j,icoarse)+tex_fr(i,j,imed)
          else if (ics == 1 .and. ifn == 1) then
              gtype(i,j,igcf)= tex_fr(i,j,icoarse)+tex_fr(i,j,ifine)
          else if (ics == 0 .and. imd+ifn == 2) then
              gtype(i,j,igmf)= tex_fr(i,j,imed)+tex_fr(i,j,ifine)
          else if (ics == 0 .and. imd == 1) then
              gtype(i,j,igm)= tex_fr(i,j,imed)
          else if (ics == 0 .and. ifn == 1) then
              gtype(i,j,igf)= tex_fr(i,j,ifine)
          endif
  !           if ((i == 802 .or. i == 803) .and. (j == 740 .or. j == 741)) then
  !                write (10,'(a,3i6,3f15.4,4e15.6)') '1 ', nz, j, i, arlat(i,j), area_ij(i,j), &
  !                   tex_fr(i,j,nz), h_relief(i,j), h_slope(i,j)
  !           endif
      endif
   enddo
enddo
!write (6,*) 'sum= ', sum

ktr= 0 ;  ktr2= 0 ;  ktr3= 0 ;  ktr4= 0
do n= 1,ngtype
   do j= 1,jd
      do i= 1,id
         if (gtype(i,j,n) /= mval) then
             ktr= ktr + 1
             if (gtype(i,j,n) < 0.99999) then
                 ktr2= ktr2 + 1
!                write (6,'(a,f8.4,3i6,2f10.3)') 'ERROR: gtype < 0, ', gtype(i,j,n), &
!                       n, j, i, lat(j), lon(i)
    !            stop 40
             endif
             if (gtype(i,j,n) < 0.90) ktr3= ktr3 + 1
             if (gtype(i,j,n) == 0.0) ktr4= ktr4 + 1
         endif
      enddo
   enddo
enddo
write (6,*) 'number of points with valid gtype=   ', ktr
write (6,*) 'number of points where gtype < 1=    ', ktr2
write (6,*) 'number of points where gtype < 0.9=  ', ktr3
write (6,*) 'number of points where gtype = 0=    ', ktr4

gmin= 9999999. ;  gmax= -9999999.
do n= 1,ngtype
   do j= 1,jd
      do i= 1,id
         if (gtype(i,j,n) /= mval) then
             gmin= min(gmin,gtype(i,j,n))
             gmax= max(gmax,gtype(i,j,n))
         endif
      enddo
   enddo
enddo
write (6,*) 'gmin= ', gmin, ', gmax= ', gmax

where (gtype > 0.) gtype= 1.

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
write (6,*) 'ground-type distribution'
write (6,*) 'sum= ', sum
write (6,*) 'sum2= ', sum2
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
if (n /= jd) then
    write (6,*) "ERROR: inconsistent lat dimension, glcc, jd= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, latl)

if (latl(1) > -89.) then
    write (6,*) "WARNING: glcc lats start at ", latl(1)
    if (latl(1) > 89.) then
        write (6,*) "  -- flip lats to start at 90S -- "
        latin= latl
        do j= 1,jd
           j1= jdp1-j
           latl(j)= latin(j1)
        enddo
    else
        stop 22
    endif
endif

latlb(1)= lat1
latlb(jdp1)= -(lat1)
do j= 2,jd
   latlb(j)= 0.5*(latl(j)+latl(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= id) then
    write (6,*) "ERROR: inconsistent lon dimension, glcc, id= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lonl)

if (lonl(1) < -179.) then
    write (6,*) "WARNING: glcc lons start at ", lonl(1)
    write (6,*) "  -- shift lons to start at zero -- "
    lonin= lonl
    do i= 1,id/2
       lonl(i)= lonin(i+id/2)
    enddo
    do i= id/2+1,id
       lonl(i)= lon2 + lonin(i-id/2)
    enddo
endif

lonlb(1)= lon1
lonlb(idp1)= lon2
do i= 2,id
   lonlb(i)= 0.5*(lonl(i)+lonl(i-1))
enddo

! compute area of latitude
arlatl= 0.
do j= 1,jd
   do i= 1,id
      arlatl(i,j)= erad*erad*(lonlb(i+1)-lonlb(i))*dtr* &
           (sin(latlb(j+1)*dtr) - sin(latlb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jd
   do i= 1,id
      sum= sum + arlatl(i,j)
   enddo
enddo
write (6,*) 'sum of glcc area= ', sum

write (10,'(/a)') 'glcc lats'
do j= 1,jd
   write (10,'(i6,5f10.3)') j, latl(j), latlb(j), latlb(j+1), arlatl(1,j)
enddo

write (10,'(/a)') 'glcc lons'
do i= 1,id
   write (10,'(i6,3f10.3)') i, lonl(i), lonlb(i), lonlb(i+1)
enddo

do j= 1,jd
   if (latl(j) /= lat(j)) then
       write (6,*) "ERROR: lats differ, ", lat(j), latl(j)
       stop 32
   endif
enddo

do i= 1,id
   if (lonl(i) /= lon(i)) then
       write (6,*) "ERROR: lons differ, ", lon(i), lonl(i)
       stop 33
   endif
enddo


write (6,*)
wbdat= mval
do n= 1,nvar_glcc
   write (6,*) 'read glcc data: ', trim(vname_glcc(n))
   rcode= nf_inq_varid (ncid, trim(vname_glcc(n)), varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find glcc variable ",  trim(vname_glcc(n))
       stop 110
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), i)
   rcode= nf_inq_dimlen (ncid, dimids(2), l)
   if (i /= id .or. l /= jd) then
       write (6,*) "ERROR: inconsistent WaterBod dimensions, glcc, id= ", i, ', jd= ', l
       stop 111
   endif
   start= 1 ;  count= 1 ;  count(1)= id ;  count(2)= jd
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
   do j= 1,jd
      j1= jdp1-j
!      write (6,*) 'j= ', j, j1
      do i= 1,id/2
         if (icov(i+id/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(icov(i+id/2,j1))*scale
         endif
      enddo
      do i= id/2+1,id
         if (icov(i-id/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(icov(i-id/2,j1))*scale
         endif
      enddo
   enddo
enddo

rcode= nf_close (ncid)

do j= 1,jd
   do i= 1,id
      if (wbdat(i,j,1) == mval .or. wbdat(i,j,2) == mval .or. &
          wbdat(i,j,3) == mval) then
          write (6,*) "WARNING: wbdat has missing value, ", wbdat(i,j,1), &
              wbdat(i,j,2), wbdat(i,j,3)
      endif
   enddo
enddo


! ----------------------------------------------------------------------
! REMOVE some snow_ice points in mexico and islands
!   lake chapala, e. mex coast, arrecife alacran (scorpion reef), 
!   clipperton island
! ----------------------------------------------------------------------

write (6,'(a)') '*** RESET snow_ice fraction to zero at following (tropical, semi-tropical) locations ***'
do j= 1,jd
   if (lat(j) > 10. .and. lat(j) < 28.) then
       do i= 1,id
          if ((lon(i) > 250. .and. lon(i) < 271.) .and. &
               wbdat(i,j,isi) > 0.) then
              write (6,'(2i6,3f12.4)') j, i, lat(j), lon(i), wbdat(i,j,isi)
              wbdat(i,j,isi)= 0.
          endif
       enddo
   endif
enddo

!  assign glcc snow/ice data to cover_type field
do j= 1,jd
   do i= 1,id
      gtype(i,j,iice)=  wbdat(i,j,isi)
      if (gtype(i,j,iice) > 0.) then
          gtype(i,j,iice)= 1.
          do n= 1,ngtype
             if (n /= iice) gtype(i,j,n)= 0.
          enddo
      endif
   enddo
enddo

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

!do j= 1,jd
!   do i= 1,id
!      if (gmask(i,j) == 1) then
!          gtype(i,j,ilake)= wbdat(i,j,iwbd)
!          if (gtype(i,j,ilake) > 0.) then
!              do n= 1,ngtype-1
!                 if (gtype(i,j,n)+gtype(i,j,ilake) > 1.) gtype(i,j,n)=  &
!                     gtype(i,j,n)-gtype(i,j,ilake)
!              enddo
!          endif
!      endif
!   enddo
!enddo


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

    lsea= mval_ls
    count(1)= id_ls ;  count(2)= jd_ls
    allocate (lsea(id_ls,jd_ls), lsea1(id_ls,jd_ls))
    rcode= nf_get_vara_double (ncid, varid, start, count, lsea)

! close netcdf file
    rcode= nf_close (ncid)

    write (6,*) 'lsea mask read'

    where (lsea /= mval_ls) lsea= 1.

    write (10,*) 'lsea'
    do j= 1,jd_ls
       write (10,'(360i1)') (int(lsea(i,j)), i= 1,id_ls)
    enddo
    
! check to see if there is a gtype cell with valid data within 15 degrees lat,lon
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea(i,j) == 1.) then
              do j2= 1,jd
                 if (abs(lat(j2)-lat_ls(j)) < 15.) then
                     do i2= 1,id
                        if (abs(lon(i2)-lon_ls(i)) < 15.) then
                            if (gmask(i2,j2) /= 0) go to 31
!                            do n= 1,ngtype
!                               if (gtype(i2,j2,n) /= 0.) go to 31
!                            enddo
                        endif
                     enddo
                 endif
              enddo
              write (6,*) 'no gtype cell found w/i 15 degrees of ', i, j, lon_ls(i), lat_ls(j)
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
    
! check to see if there is a slope cell with valid data within 18 degrees lat,lon
    do j= 1,jd_ls
       do i= 1,id_ls
          if (lsea(i,j) == 1.) then
              do j2= 1,jd
                 if (abs(lat(j2)-lat_ls(j)) < 18.) then
                     do i2= 1,id
                        if (abs(lon(i2)-lon_ls(i)) < 18.) then
                            if (gmask(i2,j2) /= 0) go to 33
                        endif
                     enddo
                 endif
              enddo
              write (6,*) 'no gtype cell found w/i 18 degrees of ', i, j, lon_ls(i), lat_ls(j)
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
             if (lat(j2) > latb_ls(j) .and. lat(j2) <= latb_ls(j+1)) then
                 do i2= 1,id
                    if (lon(i2) > lonb_ls(i) .and. lon(i2) <= lonb_ls(i+1)) then
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
