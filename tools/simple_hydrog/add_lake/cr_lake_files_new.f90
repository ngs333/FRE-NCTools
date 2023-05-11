program cp_lake_files_new

! read big lake frac file and Gascoyne prec and wroff fields

use horiz_interp_mod

implicit none

include "param.h"  !  get ntiles here

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx= 10000
integer, parameter :: nlake_def= 16, nlake_defp1= nlake_def+1
integer, parameter :: idl= 2880, jdl= 1440
integer, parameter :: idlp1= idl+1, jdlp1= jdl+1
integer, parameter :: nvar_glcc= 2
integer, parameter :: iwbd= 1, ipwt= 2
integer, parameter :: ido= 360, jdo= 180
integer, parameter :: idop1= ido+1, jdop1= jdo+1
integer, parameter :: ncell_cnv= 8
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: iaral= 15, ibalk= 13, ichad= 16, icaspian= nlake_def+1
integer, parameter :: imich= 1, ihuron= 2, ierie= 10, iont= 12

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20
real, parameter :: mval_mdl= -9999.
real, parameter :: casp_area= 3.71e11  ! m2
real, parameter :: lake_depth_large= 20.
real, parameter :: lake_depth_small= 2.
real, parameter :: lake_tau_large= 86400.*100.
real, parameter :: lake_tau_small= 1200.
real, parameter :: lake_type_real= 0.
real, parameter :: lake_type_phan= 1.
real, parameter :: sec_per_day= 86400.
real, parameter :: day_per_year= 365.25
real, parameter :: reso= 1.
real, parameter :: ldepth_chad= 2.   ! depth of lake chad in m

include '/usr/local/include/netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, nlake, rcode2, ndm, latgid, longid
integer :: varid2, varid3, varid4, varid5, varid6, varid7, i2, j2
integer :: ip, jp, fill_val, j1, i1, idx, jdx, ii, jj, idp2, jdp2
integer :: n1, np, nwb0, ncasp, varid8
real :: pi, dtr, sum, mval_frac, casp_ratio, mval_tocell, mval_landf
real :: mval_cella, mval_laka, mval_basin, ln1, lf, sum_casp, mval_lad
real :: area_land, scale, lonb1, lonb2, sum1, sum2, asum1, asum2
real :: wmax, wmin
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname, glcc_file
logical :: read_lake_frac, scale_area, read_waterbod_to_zero, cnct_next_casp_chad_outlet0
logical :: cnct_next_all_outlet0, interp_glcc_to_1deg

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, itlake, lake_idx
integer, dimension (nlkmx)             :: iwb0, jwb0, itwb0, icasp, jcasp, itcasp
integer, dimension (idl,jdl)           :: idat, iget, jget
integer, dimension (ntiles)            :: itw, ite, its, itn

real, dimension (nlkmx)                :: lake_area, lfrac, wbfrac
real, dimension (nlake_defp1)          :: sum_lake, ldepth
real, dimension (idl)                  :: lonl, lonin
real, dimension (jdl)                  :: latl, latin
real, dimension (ido)                  :: lono
real, dimension (jdo)                  :: lato
real, dimension (idlp1)                :: lonlb
real, dimension (jdlp1)                :: latlb
real, dimension (idop1)                :: lonob
real, dimension (jdop1)                :: latob
real, dimension (idl,jdl)              :: arlat
real, dimension (ido,jdo)              :: arlato
real, dimension (idl,jdl,nvar_glcc)    :: wbdat
real, dimension (ido,jdo,nvar_glcc)    :: wbd_cnv

!  Lakes (in order) are Michigan, Huron, Superior, Victoria, Tanganyika, Baikal,
!    Great Bear, Malawi, Great Slave, Erie, Winnipeg, Ontario, Balkhash, Ladoga,
!    Aral Sea, Chad, Caspian  ; areas in km2
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

real, dimension (ni_cells,nj_cells)    :: out_flow = &
    (/  8.,   4.,   2., &
       16.,   0.,   1., &
       32.,  64., 128. /)

character(len=8), dimension (nvar_glcc) :: vname_glcc= &
(/ 'WaterBod', 'PWetland' /)
   
character(len=11), dimension (nlake_defp1) :: lake_name = &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ', &
   'Caspian    ' /)

character(len=100), dimension (ntiles) :: lname_glcc, river_input_file, &
           land_month_file

integer, allocatable, dimension (:,:,:) :: ga_mask, lake_int, lake_code

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:)     :: interp_out, interp_mask
real, allocatable, dimension (:,:,:)   :: lat, lon, tocell, land_frac, cell_area, basin, &
                                          lake_frac, lake_type, lake_depth, lake_tau, &
                                          wbd, pwt, cnct_next, whole_lake

real*4, allocatable, dimension (:,:)     :: adat2

pi= 4.*atan(1.)
dtr= pi/180.

read_lake_frac= .false.
cnct_next_casp_chad_outlet0= .false.
cnct_next_all_outlet0= .false.
interp_glcc_to_1deg= .true.

open (10, file= 'out.cr_lake_files_new', form= 'formatted')
open (11, file= 'table.cr_lake_files.wbd_gt_0.1', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') read_lake_frac
read (5,'(l1)') read_waterbod_to_zero
read (5,'(l1)') cnct_next_casp_chad_outlet0
read (5,'(l1)') cnct_next_all_outlet0
read (5,'(l1)') interp_glcc_to_1deg
read (5,'(a)') glcc_file
close (5)

write (6,*) 'fort.5 data read'
write (6,*) 'read_lake_frac= ', read_lake_frac
write (6,*) 'read_waterbod_to_zero= ', read_waterbod_to_zero
write (6,*) 'cnct_next_casp_chad_outlet0= ', cnct_next_casp_chad_outlet0
write (6,*) 'cnct_next_all_outlet0= ', cnct_next_all_outlet0
write (6,*) 'interp_glcc_to_1deg= ', interp_glcc_to_1deg


! ----------------------------------------------------------------------
! read lat and lon of lake to be inserted in cover field
!   lake_area is the area of lake in the grid cell
! ----------------------------------------------------------------------

write (10,*) 'input lake info:'
open (20, form= 'formatted')
if (read_lake_frac) then
    ktr= 0
1   ktr= ktr + 1
      read (20,*,end=2) itlake(ktr), jlake(ktr), ilake(ktr), lfrac(ktr), lake_idx(ktr)
      write (10,'(3i6,f7.1,i6)') itlake(ktr), jlake(ktr), ilake(ktr), lfrac(ktr), lake_idx(ktr)
      go to 1
2   continue
    nlake= ktr - 1
else
    ktr= 0
3   ktr= ktr + 1
      read (20,*,end=4) itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr), lake_idx(ktr)
      write (10,'(3i6,f12.0,i6)') itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr), lake_idx(ktr)
      go to 3
4   continue
    nlake= ktr - 1
endif
write (6,*) 'number of lakes to be inserted = ', nlake
close (20)


! ----------------------------------------------------------------------
! read lat and lon of waterbody fraction to be reset to zero
! ----------------------------------------------------------------------

write (10,*) 'input lake info:'
if (read_waterbod_to_zero) then
    open (22, form= 'formatted')
    ktr= 0
11  ktr= ktr + 1
      read (22,*,end=12) itwb0(ktr), jwb0(ktr), iwb0(ktr), wbfrac(ktr)
      go to 11
12  continue
    nwb0= ktr - 1
endif
write (6,*) 'number of lakes to be inserted = ', nwb0
close (22)


! ----------------------------------------------------------------------
!  read glcc file -- get lats, lons, waterbod, pwetland
! ----------------------------------------------------------------------
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
if (n /= jdl) then
    write (6,*) "ERROR: inconsistent lat dimension, glcc, jdl= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jdl
rcode= nf_get_vara_double (ncid, latid, start, count, latl)

if (latl(1) > -89.) then
    write (6,*) "WARNING: glcc lats start at ", latl(1)
    if (latl(1) > 89.) then
        write (6,*) "  -- flip lats to start at 90S -- "
        latin= latl
        do j= 1,jdl
           j1= jdlp1-j
           latl(j)= latin(j1)
        enddo
    else
        stop 2
    endif
    
endif

latlb(1)= lat1
latlb(jdlp1)= -(lat1)
do j= 2,jdl
   latlb(j)= 0.5*(latl(j)+latl(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= idl) then
    write (6,*) "ERROR: inconsistent lon dimension, glcc, idl= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= idl
rcode= nf_get_vara_double (ncid, lonid, start, count, lonl)

if (lonl(1) < -179.) then
    write (6,*) "WARNING: glcc lons start at ", lonl(1)
    write (6,*) "  -- shift lons to start at zero -- "
    lonin= lonl
    do i= 1,idl/2
       lonl(i)= lonin(i+idl/2)
    enddo
    do i= idl/2+1,idl
       lonl(i)= lon2 + lonin(i-idl/2)
    enddo
endif

lonlb(1)= lon1
lonlb(idlp1)= lon2
do i= 2,idl
   lonlb(i)= 0.5*(lonl(i)+lonl(i-1))
enddo

! compute area of latitude
arlat= 0.
do j= 1,jdl
   do i= 1,idl
      arlat(i,j)= erad*erad*(lonlb(i+1)-lonlb(i))*dtr* &
           (sin(latlb(j+1)*dtr) - sin(latlb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdl
   do i= 1,idl
      sum= sum + arlat(i,j)
   enddo
enddo
write (6,*) 'sum of glcc area= ', sum

write (10,'(/a)') 'glcc lats'
do j= 1,jdl
   write (10,'(i6,5f10.3)') j, latl(j), latlb(j), latlb(j+1), arlat(1,j)
enddo

write (10,'(/a)') 'glcc lons'
do i= 1,idl
   write (10,'(i6,3f10.3)') i, lonl(i), lonlb(i), lonlb(i+1)
enddo


wbdat= mval_mdl ;   lname_glcc=''
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
   if (i /= idl .or. l /= jdl) then
       write (6,*) "ERROR: inconsistent WaterBod dimensions, glcc, id= ", i, ', jd= ', l
       stop 111
   endif
   start= 1 ;  count= 1 ;  count(1)= idl ;  count(2)= jdl
   rcode= nf_get_vara_int (ncid, varid, start, count, idat)

   fill_val= -999999
   rcode= nf_inq_attid (ncid, varid, '_FillValue', attnum)
   if (rcode == 0) rcode= nf_get_att_int (ncid, varid, '_FillValue', fill_val)
   write (6,*) 'fill_val= ', fill_val
   
   scale= -99999999.
   rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
   if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
   write (6,*) 'scale= ', scale

   rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
   if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', lname_glcc(n))
   write (6,*) 'long_name= ', trim(lname_glcc(n))
   
   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
   do j= 1,jdl
      j1= jdlp1-j
!      write (6,*) 'j= ', j, j1
      do i= 1,idl/2
         if (idat(i+idl/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(idat(i+idl/2,j1))*scale
         endif
      enddo
      do i= idl/2+1,idl
         if (idat(i-idl/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(idat(i-idl/2,j1))*scale
         endif
      enddo
   enddo
enddo

rcode= nf_close (ncid)

do j= 1,jdl
   do i= 1,idl
      if (wbdat(i,j,1) == mval_mdl .or. wbdat(i,j,2) == mval_mdl) then
          write (6,*) "WARNING: wbdat has missing value, ", wbdat(i,j,1), wbdat(i,j,2)
      endif
   enddo
enddo

sum= 0.
do j= 1,jdl
   if (latl(j) >= 12. .and. latl(j) <= 14.) then
       do i= 1,idl
          if (lonl(i) >= 13. .and. lonl(i) <= 16.) then
              sum= sum + arlat(i,j)*wbdat(i,j,1)
          endif
       enddo
   endif
enddo
write (6,*)
write (6,*) 'lake chad 0.125-deg grid point area= ', sum

! ----------------------------------------------------------------------
! try averaging 0.125-degree data to 1-degree
! ----------------------------------------------------------------------
do j= 1,jdo
   lato(j)= lat1 + reso*0.5 + real(j-1)*reso
enddo

latob(1)= lat1
latob(jdop1)= -(lat1)
do j= 2,jdo
   latob(j)= 0.5*(lato(j)+lato(j-1))
enddo

do i= 1,ido
   lono(i)= lon1 + reso*0.5 + real(i-1)*reso
enddo

lonob(1)= lon1
lonob(idop1)= lon2
do i= 2,ido
   lonob(i)= 0.5*(lono(i)+lono(i-1))
enddo

! compute area of latitude
arlato= 0.
do j= 1,jdo
   do i= 1,ido
      arlato(i,j)= erad*erad*(lonob(i+1)-lonob(i))*dtr* &
           (sin(latob(j+1)*dtr) - sin(latob(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdo
   do i= 1,ido
      sum= sum + arlato(i,j)
   enddo
enddo
write (6,*) 'sum of 1-degree area= ', sum

write (10,'(/a)') '1-degree lats'
do j= 1,jdo
   write (10,'(i6,5f10.3)') j, lato(j), latob(j), latob(j+1), arlato(1,j)
enddo

write (10,'(/a)') '1-degree lons'
do i= 1,ido
   write (10,'(i6,3f10.3)') i, lono(i), lonob(i), lonob(i+1)
enddo

if ( interp_glcc_to_1deg ) then
    wbd_cnv= mval_mdl ; sum= 0.
    do j= 1,jdo
       j1= (j-1)*ncell_cnv+1
       j2= j*ncell_cnv
       do i= 1,ido
          i1= (i-1)*ncell_cnv+1
          i2= i*ncell_cnv
          asum1= 0. ;  asum2= 0.
          sum1= 0. ;  sum2= 0.
          do jj= j1,j2
             do ii= i1,i2
                if (wbdat(ii,jj,1) /= mval_mdl) then
                    sum= sum + arlat(ii,jj)
                    asum1= asum1 + arlat(ii,jj)
                    sum1= sum1 + wbdat(ii,jj,1)*arlat(ii,jj)
                endif
                if (wbdat(ii,jj,2) /= mval_mdl) then
                    asum2= asum2 + arlat(ii,jj)
                    sum2= sum2 + wbdat(ii,jj,2)*arlat(ii,jj)
                endif
             enddo
          enddo
          if (abs(asum1 - arlato(i,j)) < 1.0e-10) then
              wbd_cnv(i,j,1)= sum1/asum1
          else
              write (6,*) "ERROR: wbd 1-deg areas do not agree, ", asum1, arlato(i,j)
              stop 17
          endif
          if (abs(asum2 - arlato(i,j)) < 1.0e-10) then
              wbd_cnv(i,j,2)= sum2/asum2
          else
              write (6,*) "ERROR: pwt 1-deg areas do not agree, ", asum2, arlato(i,j)
              stop 17
          endif
       enddo
    enddo
    write (6,*) 'sum of conversion area= ', sum
endif

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
!   cellarea, and basin
! ----------------------------------------------------------------------

allocate (lat(idp2,jdp2,ntiles), lon(idp2,jdp2,ntiles))
allocate (tocell(idp2,jdp2,ntiles), land_frac(idp2,jdp2,ntiles), cell_area(idp2,jdp2,ntiles))
allocate (basin(idp2,jdp2,ntiles))

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
   rcode= nf_get_vara_double (ncid, latid, start, count, lat(2:idp1,2:jdp1,n))
       
       
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
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2:jdp1,n))


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
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(2:idp1,2:jdp1,n))

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
   rcode= nf_get_vara_double (ncid, varid, start, count, land_frac(2:idp1,2:jdp1,n))

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
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_area(2:idp1,2:jdp1,n))

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
   

   write (6,*) 'read basin, if available'
   rcode= nf_inq_varid (ncid, 'basin', varid)     ! basin field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find basin, skipping"
       go to 75
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent basin dimension, ", ndm, id ;  stop 75
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent basin dimension, ", ndm, jd ;  stop 75
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, basin(2:idp1,2:jdp1,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_basin= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_basin)
   endif
   write (6,*) 'mval= ', mval_basin

   where (basin(:,:,n) == mval_basin) basin(:,:,n)= mval_mdl
75 continue
   

   rcode= nf_close (ncid)
   
   write (10,'(/"river lats, tile", i4)') n
   do j= 2,jdp1
      write (10,*) 'j= ', j
      write (10,'(10f10.4)') (lat(i,j,n), i= 2,idp1)
   enddo

   write (10,'(/"river lons, tile", i4)') n
   do j= 2,jdp1
      write (10,*) 'j= ', j
      write (10,'(10f10.4)') (lon(i,j,n), i= 2,idp1)
   enddo

   write (10,'(/"cell_area, tile", i4)') n
   do j= 2,jdp1
      write (10,*) 'j= ', j
      write (10,'(10f10.1)') (cell_area(i,j,n)/1.e6, i= 2,idp1)
   enddo

enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         sum= sum + cell_area(i,j,n)
      enddo
   enddo
enddo
write (6,'(/"sum of cell_area= ",f15.1)') sum/1.e+6

scale_area = .true.
if (sum/1.e6 < 510064460.) scale_area = .false.
write (6,*) 'scale_area= ', scale_area

       
!  write lats and lons of input lakes
write (6,'(/"input lake coordinates:")')
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   ln1= lon(i+1,j+1,n)
   if (ln1 > 180.) ln1= ln1 - 360.
!   write (6,'(4i6,2f10.2,f15.0)') l, n, i, j, ln1, lat(i+1,j+1,n), basin(i+1,j+1,n)
enddo


! compute total land area, excluding Greenland and Antarctica
allocate (ga_mask(idp2,jdp2,ntiles))
ga_mask= 0
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if ( lat(i,j,n) > 59.  .and. &
             (lon(i,j,n) > 298. .and. lon(i,j,n) < 350.)) ga_mask(i,j,n)= 1
         if (lat(i,j,n) < -60.) ga_mask(i,j,n)= 1
      enddo
   enddo
enddo

area_land= 0.
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         if (ga_mask(i,j,n) == 0) then
             if (scale_area .and. land_frac(i,j,n) /= mval_mdl) then 
                 area_land= area_land + cell_area(i,j,n)*land_frac(i,j,n)
             else
                 area_land= area_land + cell_area(i,j,n)
             endif
         endif
      enddo
   enddo
enddo
write (6,'(/"sum of land area, no greenland or antarctica= ",f15.1)') area_land/1.e+6


allocate (lake_frac(idp2,jdp2,ntiles), lake_type(idp2,jdp2,ntiles))
allocate (lake_depth(idp2,jdp2,ntiles), lake_tau(idp2,jdp2,ntiles))

lake_frac= mval_mdl
lake_type= mval_mdl

where (land_frac == 1.) lake_frac= 0.
!where (land_frac > 0.) lake_frac= 0.


! ----------------------------------------------------------------------
!  Caspian Sea
! ----------------------------------------------------------------------
allocate (lake_int(idp2,jdp2,ntiles), lake_code(idp2,jdp2,ntiles))

lake_int= 0 ;  lake_code= 0
write (6,'(/"Caspian lake fractions")')
open (21, form= 'formatted')
ktr= 0 ;  sum_casp= 0.
200 ktr= ktr + 1
    read (21,*,end= 202) itcasp(ktr), jcasp(ktr), icasp(ktr), lf
    n= itcasp(ktr) ;  j= jcasp(ktr) ;  i= icasp(ktr)
    if (cell_area(i+1,j+1,n) == mval_mdl) then
        write (6,*) "ERROR: missing cell_area, ", n, j, i ; stop 200
    endif
    lake_frac(i+1,j+1,n)= lf
    lake_int(i+1,j+1,n)= 1
    lake_code(i+1,j+1,n)= icaspian
    if ( cnct_next_casp_chad_outlet0 .or. cnct_next_all_outlet0 ) &
         lake_int(i+1,j+1,n)= 0
    sum_casp= sum_casp + lf*cell_area(i+1,j+1,n)
!    write (6,'(4i6,2f10.2,f12.1,f10.3)') ktr, n, i, j, lon(i+1,j+1,n), lat(i+1,j+1,n), &
!        lf*cell_area(i+1,j+1,n)/1.e6, lake_frac(i+1,j+1,n)
    go to 200
202 continue
write (6,*) 'number of Caspian points= ', ktr-1
ncasp= ktr-1

write (6,*) "obs caspian area= ", casp_area/1.e6
write (6,*) "mdl caspian area= ", sum_casp/1.e6
if (sum_casp > casp_area) then
    write (6,*) "ERROR: model caspian area exceeds obs"
endif
           

! ----------------------------------------------------------------------
!  now insert lakes in lake_frac field
!    do not alter Caspian fractions -- already present in frac field
! ----------------------------------------------------------------------

write (10,'(/"input lake fractions")')
write (6,'(/"input lake fractions")')
if (read_lake_frac) then
    sum_lake= 0.
    do l= 1,nlake
       i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
       if (lake_frac(i+1,j+1,n) > 0.) write (6,'(a,i6,3f10.2)') "lake already present at ", &
           n, lon(i+1,j+1,n), lat(i+1,j+1,n), lake_frac(i+1,j+1,n)
       lake_frac(i+1,j+1,n)= lfrac(l)
       sum_lake(lake_idx(l))= sum_lake(lake_idx(l)) + lfrac(l)*cell_area(i+1,j+1,n)
       if (lake_idx(l) == iaral .or. lake_idx(l) == ibalk .or. lake_idx(l) == ichad) then
           lake_int(i+1,j+1,n)= 1
       endif
       if ( cnct_next_casp_chad_outlet0 .and. lake_idx(l) == ichad) lake_int(i+1,j+1,n)= 0
       if ( cnct_next_all_outlet0 ) lake_int(i+1,j+1,n)= 0
       lake_code(i+1,j+1,n)= lake_idx(l)
   
       write (10,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i+1,j+1,n), lat(i+1,j+1,n), &
            lfrac(l)*cell_area(i+1,j+1,n)/1.e6, cell_area(i+1,j+1,n)/1.e6, lfrac(l)
    enddo
    do l= 1,nlake_def
       write (6,'(i6,3x,a11,2f12.0)') l, lake_name(l), lake_def_area(l), sum_lake(l)/1.e6
    enddo
else
    do l= 1,nlake
       i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
       if (lake_frac(i+1,j+1,n) > 0.) write (6,'(a,i6,3f10.2)') "lake already present at ", &
           n, lon(i+1,j+1,n), lat(i+1,j+1,n), lake_frac(i+1,j+1,n)
       lake_frac(i+1,j+1,n)= lake_area(l)/(cell_area(i+1,j+1,n)/1.e6)
       if (lake_idx(l) == iaral .or. lake_idx(l) == ibalk .or. lake_idx(l) == ichad) then
           lake_int(i+1,j+1,n)= 1
       endif
       if ( cnct_next_casp_chad_outlet0 .and. lake_idx(l) == ichad) lake_int(i+1,j+1,n)= 0
       if ( cnct_next_all_outlet0 ) lake_int(i+1,j+1,n)= 0
       lake_code(i+1,j+1,n)= lake_idx(l)
   
       write (10,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i+1,j+1,n), lat(i+1,j+1,n), &
            lake_area(l), cell_area(i+1,j+1,n)/1.e6, lake_area(l)/(cell_area(i+1,j+1,n)/1.e6)
       if (lake_area(l) > cell_area(i+1,j+1,n)/1.e6) then
           write (6,'(a,2f10.2,2f12.0)') "ERROR: lake area exceeds grid cell area, ", &
             lon(i+1,j+1,n), lat(i+1,j+1,n), lake_area(l), cell_area(i+1,j+1,n)/1.e6
            stop 40
       endif
    enddo
endif

! check to make sure lake does not occur at partial-land cell
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (lake_frac(i,j,n) > 0. .and. land_frac(i,j,n) < 1.) then
             write (6,*) 'ERROR: lake occurs at partial-land cell'
             write (6,*) 'n= ', n, ', i= ', i-1, ', j= ', j-1
             stop 18
         endif
      enddo
   enddo
enddo

!  set lake_frac to zero (currently missing value) at partial-land cells
where (land_frac > 0. .and. land_frac < 1.) lake_frac= 0.


! ----------------------------------------------------------------------
! interpolate glcc water bodies and permanent wetlands to model grid
! ----------------------------------------------------------------------

allocate (wbd(idp2,jdp2,ntiles), pwt(idp2,jdp2,ntiles), interp_out(id,jd))

if (interp_glcc_to_1deg) then
    allocate (interp_mask(ido,jdo))
    interp_mask= 1.
    where (wbd_cnv(:,:,1) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
       write (6,*) 'WaterBod, tile= ', n
       call horiz_interp (wbd_cnv(:,:,1), lonob*dtr, latob*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, wbd(2:idp1,2:jdp1,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       where (interp_out(:,:) == 0.) wbd(2:idp1,2:jdp1,n)= mval_mdl
    enddo
   
    interp_mask= 1.
    where (wbd_cnv(:,:,2) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
       write (6,*) 'PWetland, tile= ', n
       call horiz_interp (wbd_cnv(:,:,2), lonob*dtr, latob*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, pwt(2:idp1,2:jdp1,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       where (interp_out(:,:) == 0.) pwt(2:idp1,2:jdp1,n)= mval_mdl
    enddo

else
    allocate (interp_mask(idl,jdl))
    interp_mask= 1.
    where (wbdat(:,:,1) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
       write (6,*) 'WaterBod, tile= ', n
       call horiz_interp (wbdat(:,:,1), lonlb*dtr, latlb*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, wbd(2:idp1,2:jdp1,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       where (interp_out(:,:) == 0.) wbd(2:idp1,2:jdp1,n)= mval_mdl
    enddo
   
    interp_mask= 1.
    where (wbdat(:,:,2) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
       write (6,*) 'PWetland, tile= ', n
       call horiz_interp (wbdat(:,:,2), lonlb*dtr, latlb*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, pwt(2:idp1,2:jdp1,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       where (interp_out(:,:) == 0.) pwt(2:idp1,2:jdp1,n)= mval_mdl
    enddo
endif

deallocate (interp_out, interp_mask)

! set ocean values to missing
where (land_frac == 0.) 
   wbd= mval_mdl
   pwt= mval_mdl
endwhere

! set part-land values to zero
where (land_frac > 0 .and. land_frac < 1.)
   wbd= 0.
   pwt= 0.
endwhere

! check to see if 'waterbod' occurs at grid cells where lake_frac > 0
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if ( lake_frac(i,j,n) > 0 .and. &
             (wbd(i,j,n) == 0. .or. wbd(i,j,n) == mval_mdl)) then
              write (6,'(a,f9.4,a,f9.4,3i6)') 'WARNING: lake_frac= ', &
                lake_frac(i,j,n), ', waterbod= ', wbd(i,j,n), n, j-1, i-1
         endif
      enddo
   enddo
enddo

! now set specified waterbod fractions (coastal) to zero
if (read_waterbod_to_zero) then
    do l= 1,nwb0
       i= iwb0(l) ;  j= jwb0(l) ;  n= itwb0(l)
       if (abs(wbfrac(l)-wbd(i+1,j+1,n)) > 0.001) then
           write (6,*) "ERROR: check wb0 fraction"
           write (6,'(4i6,2f8.4)') l, n, j, i, wbfrac(l), wbd(i+1,j+1,n)
           stop 140
       endif
       wbd(i+1,j+1,n)= 0.
    enddo
endif

! compute lake area
sum_lake= 0.
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   sum_lake(lake_idx(l))= sum_lake(lake_idx(l)) + wbd(i+1,j+1,n)*cell_area(i+1,j+1,n)
!   write (6,'(4i6,2f10.2,2f12.0,f12.3)') l, i, j, lake_idx(l), lon(i+1,j+1,n), lat(i+1,j+1,n), &
!        wbd(i+1,j+1,n)*cell_area(i+1,j+1,n)/1.e6, cell_area(i+1,j+1,n)/1.e6, wbd(i+1,j+1,n)
enddo
    
sum_casp= 0.
do l= 1,ncasp
   i= icasp(l) ;  j= jcasp(l) ;  n= itcasp(l)
   sum_lake(icaspian)= sum_lake(icaspian) + wbd(i+1,j+1,n)*cell_area(i+1,j+1,n)
!   write (6,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i+1,j+1,n), lat(i+1,j+1,n), &
!        wbd(i+1,j+1,n)*cell_area(i+1,j+1,n)/1.e6, cell_area(i+1,j+1,n)/1.e6, wbd(i+1,j+1,n)
enddo

!  compute lake depth
do l= 1,nlake_defp1
   ldepth(l)= lake_def_vol_wenc(l)/lake_def_area_wenc(l)*1.e3
enddo
ldepth(ichad)= ldepth_chad

do l= 1,nlake_defp1
   write (6,'(i6,3x,a11,2f12.0,f12.2)') l, lake_name(l), lake_def_area_wenc(l), sum_lake(l)/1.e6, ldepth(l)
enddo

!  add michigan and huron
sum= sum_lake(imich)+ sum_lake(ihuron)
sum_lake(imich)= sum
sum_lake(ihuron)= sum

sum=  lake_def_vol_wenc(imich)+lake_def_vol_wenc(ihuron)
sum2= lake_def_area_wenc(imich)+lake_def_area_wenc(ihuron)
ldepth(imich)= sum/sum2*1.e3
ldepth(ihuron)= sum/sum2*1.e3

write (6,*)
write (6,*) 'add michigan and huron'
do l= 1,nlake_defp1
   write (6,'(i6,3x,a11,2f12.0,f12.2)') l, lake_name(l), lake_def_area_wenc(l), sum_lake(l)/1.e6, ldepth(l)
enddo

! define whole_lake_area field and lake depth field
allocate (whole_lake(idp2,jdp2,ntiles))
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         whole_lake(i,j,n)= wbd(i,j,n)*cell_area(i,j,n)
         lake_depth(i,j,n)= lake_depth_small
         lake_tau(i,j,n)= 0.
      enddo
   enddo
enddo

do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   whole_lake(i+1,j+1,n)= sum_lake(lake_idx(l))
   lake_depth(i+1,j+1,n)= ldepth(lake_idx(l))
   lake_tau(i+1,j+1,n)= lake_tau_large
enddo
do l= 1,ncasp
   i= icasp(l) ;  j= jcasp(l) ;  n= itcasp(l)
   whole_lake(i+1,j+1,n)= sum_lake(icaspian)
   lake_depth(i+1,j+1,n)= ldepth(icaspian)
   lake_tau(i+1,j+1,n)= lake_tau_large
enddo

where (land_frac == 0.) whole_lake= mval_mdl
where (land_frac == 0.) lake_depth= mval_mdl
where (land_frac == 0.) lake_tau= mval_mdl

!where (lake_frac > 0.) 
!   lake_tau= lake_tau_large
!endwhere
!where (lake_frac == 0.)
!   lake_tau=    0.
!endwhere

!where (land_frac < 1. .and. land_frac > 0.)
!   lake_tau=    0.
!endwhere


! ----------------------------------------------------------------------
! create 'connected_to_next' field, showing whether grid cells within
!   a large lake are connected
! ----------------------------------------------------------------------

! set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

if (ntiles == 1) then
    itw(n)= 1 ;  ite(n)= 1 ;  its(n)= 1 ;  itn(n)= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_frac)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_area)
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
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_frac)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_area)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
endif

allocate (cnct_next(idp2,jdp2,ntiles))

where (land_frac == 0.) 
   cnct_next= mval_mdl
elsewhere
   cnct_next= 0.
endwhere
where (lake_frac > 0.) cnct_next= 1.

! for computing connected_to_next field, combine lakes michigan and huron in lake_code field
where (lake_code == ihuron) lake_code= imich

do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 170
             ktr= 0
             n1= n
          
             if (n == 3 .and. i == 16 .and. j == 2) then
                 write (6,'(a,3i6,3f10.2)') 'aral, ', n, i, j, tocell(i,j,n), &
                     basin(i,j,n), lake_frac(i,j,n)
             endif
             do jj= 1,nj_cells
                jp= j+jj-2
                
                do ii= 1,ni_cells
                   ip=i+ii-2
                
                if (tocell(i,j,n) == out_flow(ii,jj)) then
                    ktr= ktr + 1
                    
                    if (tocell(i,j,n) == 0.) then
                        if (lake_frac(i,j,n) > 0.) then
                            write (6,*) 'ERROR: lake occurs at partial-land cell'
                            stop 170
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
                        if (n == 3 .and. i == 16 .and. j == 2) then
                            write (6,'(a,6i6,2f10.2,i6)') '1 ', n, i, j, n1, i1, j1, &
                              lake_frac(i,j,n),  lake_frac(i1,j1,n1)
                        endif
!                        if (lake_frac(i,j,n) > 0. .and. lake_frac(i1,j1,n1) == 0. .and.  &
                        if (lake_code(i,j,n) /= lake_code(i1,j1,n1) .and. &
                            lake_int(i,j,n) == 0) then
                            cnct_next(i,j,n)= 0.
                        endif
                        if (n == 3 .and. i == 16 .and. j == 2) then
                            write (6,'(a,6i6,2f10.2,i6)') '2 ', n, i, j, n1, i1, j1, &
                           lake_frac(i,j,n),  lake_frac(i1,j1,n1)
                        endif
                    endif
                endif
                
                enddo
             enddo
      
170      continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop
write (6,*) 'connected_to_next field computed'


! output all full-land grid cells where waterbod exceeds 0.1
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (wbd(i,j,n) > 0.1 .and. land_frac(i,j,n) == 1.) then
             do jj= 1,nj_cells
                jp= j+jj-2
                do ii= 1,ni_cells
                   ip=i+ii-2
                   if (land_frac(ip,jp,n) < 1.) go to 175
                enddo
             enddo
             write (11,'(3i5,f9.4,3x,a2)') n, j-1, i-1, wbd(i,j,n), "**"
             go to 176
175          continue
             write (11,'(3i5,f9.4,3x,a2)') n, j-1, i-1, wbd(i,j,n)
176          continue
         endif
      enddo
   enddo
enddo


! check all waterbod fractions where lake_frac > 0.
wmax= -9999999. ;  wmin= 9999999.
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (lake_frac(i,j,n) > 0. .and. cnct_next(i,j,n) == 0.) then
             wmax= max(wmax,wbd(i,j,n))
             wmin= min(wmin,wbd(i,j,n))
             if (wbd(i,j,n) < 0.1) write (6,'(a,3i5,2f9.4)') &
                "***WARNING***  Outlet cell has small lake fraction, ", &
                   n, j-1, i-1, lake_frac(i,j,n), wbd(i,j,n)
         endif
      enddo
   enddo
enddo
write (6,*) 'wmax= ', wmax, ', wmin= ', wmin


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'lake_frac.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))
   
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (grid_x, grid_y)
   rcode= NF_DEF_DIM (ncid, 'grid_x',  id,   londim)
   rcode= NF_DEF_DIM (ncid, 'grid_y',  jd,   latdim)

!  create coordinate variables
   rcode= NF_DEF_VAR (ncid, 'grid_x',  NF_DOUBLE, 1, londim,  lonid)
   rcode= NF_DEF_VAR (ncid, 'grid_y',  NF_DOUBLE, 1, latdim,  latid)

!  create attributes for coordinate variables
!    longitude:
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 16, 'T-cell longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'axis', 1, 'X')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')

!    latitude:
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 15, 'T-cell latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'axis', 1, 'Y')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')

!    create data variable and attributes
   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 9, 'degrees_E')
 
   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 9, 'degrees_N')
 
   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'lake_frac', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 13, 'lake_fraction')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'lake_depth_sill', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 15, 'lake_depth_sill')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'lake_tau', NF_DOUBLE, 2, ndims, varid4)
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 8, 'lake_tau')
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 1, 's')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
   rcode= NF_DEF_VAR (ncid, trim(vname_glcc(1)), NF_DOUBLE, 2, ndims, varid5)
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', len_trim(lname_glcc(1)), &
          trim(lname_glcc(1)))
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, trim(vname_glcc(2)), NF_DOUBLE, 2, ndims, varid6)
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', len_trim(lname_glcc(2)), &
          trim(lname_glcc(2)))
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
   rcode= NF_DEF_VAR (ncid, 'connected_to_next', NF_DOUBLE, 2, ndims, varid7)
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'long_name', 20, 'lake connection flag')
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid7, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
   rcode= NF_DEF_VAR (ncid, 'whole_lake_area', NF_DOUBLE, 2, ndims, varid8)
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'long_name', 18, 'total area of lake')
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'units', 2, 'm2')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid8, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
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

!    lake fraction data
   start= 1 ;  count(1)= id ;  count(2)= jd 
!   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, lake_frac(2:idp1,2:jdp1,n))
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, wbd(2:idp1,2:jdp1,n))

!    lake depth data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, lake_depth(2:idp1,2:jdp1,n))

!    lake tau data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, lake_tau(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, wbd(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, pwt(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, cnct_next(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid8, start, count, whole_lake(2:idp1,2:jdp1,n))

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

close (10)
close (11)

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area, basin)
deallocate (lake_frac, lake_depth, lake_type, lake_tau)
deallocate (wbd, pwt, cnct_next, whole_lake, ga_mask, lake_int)

   
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



