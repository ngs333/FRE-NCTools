program cp_small_lake_frac

! read big lake frac file and Gascoyne prec and wroff fields

use horiz_interp_mod

implicit none

include "param.h"  !  get ntiles here

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx= 1000
integer, parameter :: nlake_def= 16
integer, parameter :: idl= 360, jdl= 180
integer, parameter :: idlp1= idl+1, jdlp1= jdl+1

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20
real, parameter :: mval_mdl= -9999.
real, parameter :: casp_area= 3.71e11  ! m2
real, parameter :: lake_depth_real= 20.
real, parameter :: lake_depth_phan= 2.
real, parameter :: lake_tau_real= 86400.*100.
real, parameter :: lake_tau_phan= 1200.
real, parameter :: lake_type_real= 0.
real, parameter :: lake_type_phan= 1.
real, parameter :: sec_per_day= 86400.
real, parameter :: day_per_year= 365.25
real, parameter :: F_lake= 0.03
real, parameter :: prox_ratio= 5.

include '/usr/local/include/netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, nlake, rcode2, ndm, latgid, longid
integer :: varid2, varid3, varid4, varid5, varid6, varid7, i2, j2
integer :: ip, jp
real :: pi, dtr, sum, mval_frac, casp_ratio, mval_tocell, mval_landf
real :: mval_cella, mval_laka, mval_basin, ln1, lf, sum_casp, mval_lad
real :: area_land, area_lakecell_b, area_lakecell_s, area_lake_b
real :: area_lake_s, F_lake_big, F_lake_small, rr_mn_small, r, r1
real :: dlat2, dlon2, p_const, fmin, fmax
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: long_name, fname, ladworld_file
logical :: read_lake_frac, scale_area

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, itlake, lake_idx
integer, dimension (idl,jdl)           :: iget, jget

real, dimension (nlkmx)                :: lake_area, lfrac
real, dimension (nlake_def)            :: sum_lake
real, dimension (idl)                  :: lonl, rlon
real, dimension (jdl)                  :: latl, rlat
real, dimension (idlp1)                :: lonlb
real, dimension (jdlp1)                :: latlb
real, dimension (idl,jdl)              :: precl, wroffl, rratiol, interp_mask

!  Lakes (in order) are Michigan, Huron, Superior, Victoria, Tanganyika, Baikal,
!    Great Bear, Malawi, Great Slave, Erie, Winnipeg, Ontario, Balkhash, Ladoga,
!    Aral Sea, Chad
real, dimension (nlake_def)              :: lake_def_area = &
(/ 58020., 60700., 83270., 68800., 32900., 30500., 31790., 28900., 28440., &
   25680., 24510., 19230., 17400., 18390., 68000., 26000. /)
   
character(len=11), dimension (nlake_def) :: lake_name = &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ' /)

character(len=100), dimension (ntiles) :: river_input_file, land_month_file

integer, allocatable, dimension (:,:,:) :: ga_mask

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:)     :: interp_out
real, allocatable, dimension (:,:,:)   :: lat, lon, tocell, land_frac, cell_area, basin, &
                                          lake_frac, lake_type, lake_depth, lake_tau, &
                                          prec, wroff, rratio

real*4, allocatable, dimension (:,:)     :: adat2

pi= 4.*atan(1.)
dtr= pi/180.

read_lake_frac= .false.

open (10, file= 'out.cp_small_lake_frac', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') read_lake_frac
read (5,'(a)') ladworld_file
close (5)

write (6,*) 'fort.5 data read'
write (6,*) 'read_lake_frac= ', read_lake_frac


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
      read (20,*,end=4) itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr)
      write (10,'(3i6,f12.0)') itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr)
      go to 3
4   continue
    nlake= ktr - 1
endif
write (6,*) 'number of lakes to be inserted = ', nlake
close (20)


! ----------------------------------------------------------------------
!  read ladworld file -- get lats, lons, edges, precip, wroff
! ----------------------------------------------------------------------
rcode= NF_OPEN (trim(ladworld_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open ladworld netcdf file"  
    write (6,*) trim(ladworld_file)
    stop 100
endif

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld lat variable" ; stop 102
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jdl) then
    write (6,*) "ERROR: inconsistent lat dimension, ladworld, jdl= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jdl
rcode= nf_get_vara_double (ncid, latid, start, count, latl)

rcode= nf_inq_varid (ncid, 'latb', latid)         ! number of lat edges
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld latb variable" ; stop 102
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jdlp1) then
    write (6,*) "ERROR: inconsistent latb dimension, ladworld, jdlp1= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jdlp1
rcode= nf_get_vara_double (ncid, latid, start, count, latlb)

rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= idl) then
    write (6,*) "ERROR: inconsistent lon dimension, ladworld, idl= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= idl
rcode= nf_get_vara_double (ncid, lonid, start, count, lonl)

rcode= nf_inq_varid (ncid, 'lonb', lonid)         ! number of lon edges
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld lonb variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= idlp1) then
    write (6,*) "ERROR: inconsistent lonb dimension, ladworld, idlp1= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= idlp1
rcode= nf_get_vara_double (ncid, lonid, start, count, lonlb)

write (10,'(/a)') 'ladworld lats:'
do j= 1,jdl
   write (10,'(i6,3f8.1)') j, latl(j), latlb(j), latlb(j+1)
enddo

write (10,'(/a)') 'ladworld lons:'
do i= 1,idl
   write (10,'(i6,3f8.1)') i, lonl(i), lonlb(i), lonlb(i+1)
enddo


write (6,*) 'read ladworld precip'
rcode= nf_inq_varid (ncid, 'precip', varid)         ! precipitation
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld precip variable" ; stop 110
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
rcode= nf_inq_dimlen (ncid, dimids(2), l)
if (n /= idl .or. l /= jdl) then
    write (6,*) "ERROR: inconsistent precip dimensions, ladworld, id= ", n, ', jd= ', l
    stop 111
endif
start= 1 ;  count= 1 ;  count(1)= idl ;  count(2)= jdl
rcode= nf_get_vara_double (ncid, varid, start, count, precl)

var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
write (6,*) 'units= ', var_units

mval_lad= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_lad)
endif
write (6,*) 'mval= ', mval_lad
where (precl == mval_lad) 
    precl=  mval_mdl
elsewhere
    precl= precl*sec_per_day*day_per_year
endwhere
   
write (6,*) 'read ladworld wroff'
rcode= nf_inq_varid (ncid, 'wroff', varid)         ! wroff
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find ladworld wroff variable" ; stop 112
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
rcode= nf_inq_dimlen (ncid, dimids(2), l)
if (n /= idl .or. l /= jdl) then
    write (6,*) "ERROR: inconsistent wroff dimensions, ladworld, id= ", n, ', jd= ', l
    stop 113
endif
start= 1 ;  count= 1 ;  count(1)= idl ;  count(2)= jdl
rcode= nf_get_vara_double (ncid, varid, start, count, wroffl)

var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
write (6,*) 'units= ', var_units

mval_lad= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_lad)
endif
write (6,*) 'mval= ', mval_lad
where (wroffl == mval_lad) 
    wroffl=  mval_mdl
elsewhere
    wroffl= wroffl*sec_per_day*day_per_year
endwhere
   
rcode= nf_close (ncid)
   
write (10,'(/a)') 'ladworld precip:'
do j= 1,jdl
   write (10,*) 'j= ', j
   write (10,'(12f9.2)') (precl(i,j), i= 1,idl)
enddo
   
write (10,'(/a)') 'ladworld wroff:'
do j= 1,jdl
   write (10,*) 'j= ', j
   write (10,'(12f9.2)') (wroffl(i,j), i= 1,idl)
enddo



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

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)
  
rcode= nf_close (ncid)


! ----------------------------------------------------------------------
! now open river files -- read lat,lon grids, tocell, land_frac, 
!   cellarea, and basin
! ----------------------------------------------------------------------

allocate (lat(id,jd,ntiles), lon(id,jd,ntiles))
allocate (tocell(id,jd,ntiles), land_frac(id,jd,ntiles), cell_area(id,jd,ntiles))
allocate (basin(id,jd,ntiles))

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
       write (6,*) "ERROR: inconsistent lon dimension, ", ndm, jd ;  stop 35
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
   rcode= nf_get_vara_double (ncid, varid, start, count, basin(:,:,n))

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
   do j= 1,jd
      write (10,*) 'j= ', j
      write (10,'(10f10.4)') (lat(i,j,n), i= 1,id)
   enddo

   write (10,'(/"river lons, tile", i4)') n
   do j= 1,jd
      write (10,*) 'j= ', j
      write (10,'(10f10.4)') (lon(i,j,n), i= 1,id)
   enddo

   write (10,'(/"cell_area, tile", i4)') n
   do j= 1,jd
      write (10,*) 'j= ', j
      write (10,'(10f10.1)') (cell_area(i,j,n)/1.e6, i= 1,id)
   enddo

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

scale_area = .true.
if (sum/1.e6 < 510064460.) scale_area = .false.
write (6,*) 'scale_area= ', scale_area


!  write lats and lons of input lakes
write (6,'(/"input lake coordinates:")')
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   ln1= lon(i,j,n)
   if (ln1 > 180.) ln1= ln1 - 360.
   write (6,'(4i6,2f10.2,f15.0)') l, n, i, j, ln1, lat(i,j,n), basin(i,j,n)
enddo


! ----------------------------------------------------------------------
!  get prec and wroff at ocean points -- this is to have values where
!    interpolation leaves missing values at model coastal points
! ----------------------------------------------------------------------

do j= 1,jdl
   do i= 1,idl
      if ((precl(i,j) /= mval_mdl .and. wroffl(i,j) == mval_mdl) .or. &
          (precl(i,j) == mval_mdl .and. wroffl(i,j) /= mval_mdl)) then
          write (6,*) 'prec and wroff fields inconsistent w.r.t. missing data', &
                       i, j, precl(i,j), wroffl(i,j)
      endif
   enddo
enddo


rlat= latl*dtr
rlon= lonl*dtr

iget= 0 ;  jget= 0
do j= 1,jdl
   do i= 1,idl
      if (precl(i,j) /= mval_mdl) go to 86
      r= (8.*pi)**2    ! initialize r to some big value

      do j2= 1,jdl
         dlat2= rlat(j) - rlat(j2)
             
         do i2= i,idl
            if (precl(i2,j2) /= mval_mdl) then
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
                endif
            endif
         enddo
      enddo
      iget(i,j)= ip
      jget(i,j)= jp
86    continue
   enddo
enddo

!  now fill missing values
do j= 1,jdl
   do i= 1,idl
      if (precl(i,j) == mval_mdl) then
          if (iget(i,j) == 0 .or. jget(i,j) == 0) then
              write (6,*) "ERROR: iget or jget= 0, prec, ", i, j, iget(i,j), jget(i,j)
              stop 86
          endif
          precl(i,j)= precl(iget(i,j),jget(i,j))
      endif
      if (wroffl(i,j) == mval_mdl) then
          if (iget(i,j) == 0 .or. jget(i,j) == 0) then
              write (6,*) "ERROR: iget or jget= 0, wroff, ", i, j, iget(i,j), jget(i,j)
              stop 86
          endif
          wroffl(i,j)= wroffl(iget(i,j),jget(i,j))
      endif
   enddo
enddo
    
! for now, don't allow precip= zero
where (precl == 0.) precl= 0.00001

! for now, don't allow negative wroff -- set to zero
where (wroffl < 0.) wroffl= 0.

! compute runoff ratio
rratiol= mval_mdl
do j= 1,jdl
   do i= 1,idl
      if (precl(i,j) /= mval_mdl .and. wroffl(i,j) /= mval_mdl) then
          rratiol(i,j)= wroffl(i,j)/precl(i,j)
      endif
   enddo
enddo

! ----------------------------------------------------------------------
! interpolate precip and wroff to model grid
! ----------------------------------------------------------------------

interp_mask= 1.
where (rratiol == mval_mdl) interp_mask= 0.

allocate (prec(id,jd,ntiles), wroff(id,jd,ntiles), rratio(id,jd,ntiles), interp_out(id,jd))

do n= 1,ntiles
   write (6,*) 'tile= ', n
   call horiz_interp (rratiol, lonlb*dtr, latlb*dtr, lon(:,:,n)*dtr, lat(:,:,n)*dtr, &
        rratio(:,:,n), verbose=1, mask_in=interp_mask, mask_out=interp_out, &
        interp_method="bilinear")
   where (interp_out(:,:) == 0.) prec(:,:,n)= mval_mdl
   
!   call horiz_interp (precl, lonlb*dtr, latlb*dtr, lon(:,:,n)*dtr, lat(:,:,n)*dtr, &
!        prec(:,:,n), verbose=1, mask_in=interp_mask, mask_out=interp_out, &
!        interp_method="bilinear")
!   where (interp_out(:,:) == 0.) wroff(:,:,n)= mval_mdl
   
!   call horiz_interp (wroffl, lonlb*dtr, latlb*dtr, lon(:,:,n)*dtr, lat(:,:,n)*dtr, &
!        wroff(:,:,n), verbose=1, mask_in=interp_mask, mask_out=interp_out, &
!        interp_method="bilinear")
!   where (interp_out(:,:) == 0.) wroff(:,:,n)= mval_mdl
enddo

deallocate (interp_out)

! for now, don't allow precip= zero
!where (prec == 0.) prec= 0.00001

! set ocean values to missing
where (land_frac == 0)
   prec= mval_mdl
   wroff= mval_mdl
   rratio= mval_mdl
endwhere

! compute runoff ratio
!rratio= mval_mdl
!do n= 1,ntiles
!   do j= 1,jd
!      do i= 1,id
!         if (prec(i,j,n) /= mval_mdl .and. wroff(i,j,n) /= mval_mdl) then
!             rratio(i,j,n)= wroff(i,j,n)/prec(i,j,n)
!         endif
!      enddo
!   enddo
!enddo

! compute total land area, excluding Greenland and Antarctica
allocate (ga_mask(id,jd,ntiles))
ga_mask= 0
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if ( lat(i,j,n) > 59.  .and. &
             (lon(i,j,n) > 298. .and. lon(i,j,n) < 350.)) ga_mask(i,j,n)= 1
         if (lat(i,j,n) < -60.) ga_mask(i,j,n)= 1
      enddo
   enddo
enddo

area_land= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
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


! ----------------------------------------------------------------------
!  insert lake at every internal drain point
! ----------------------------------------------------------------------

allocate (lake_frac(id,jd,ntiles), lake_type(id,jd,ntiles), lake_depth(id,jd,ntiles))
allocate (lake_tau(id,jd,ntiles))

lake_frac= mval_mdl
lake_type= mval_mdl
lake_depth= mval_mdl
lake_tau= mval_mdl

where (land_frac == 1.) lake_frac= 0.
!where (land_frac > 0.) lake_frac= 0.

ktr= 0
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (tocell(i,j,n) == 0. .and. land_frac(i,j,n) == 1.) then
             ktr= ktr + 1
             lake_frac(i,j,n)= lfrac_const
         endif
      enddo
   enddo
enddo

write (6,*)
write (6,*) 'number of internal drain points= ', ktr
write (6,*) 'this number includes Caspian drain point(s), lake_frac will be overwritten'
write (6,*) 'internal drain points assigned lake fraction= ', lfrac_const


! ----------------------------------------------------------------------
!  use model-output lake fraction for Caspian Sea
! ----------------------------------------------------------------------
write (6,'(/"Caspian lake fractions")')
open (21, form= 'formatted')
ktr= 0 ;  sum_casp= 0.
200 ktr= ktr + 1
    read (21,*,end= 202) n, j, i, lf
    if (cell_area(i,j,n) == mval_mdl) then
        write (6,*) "ERROR: missing cell_area, ", n, j, i ; stop 200
    endif
    lake_frac(i,j,n)= lf
    sum_casp= sum_casp + lf*cell_area(i,j,n)
    write (6,'(4i6,2f10.2,f12.1,f10.3)') ktr, n, i, j, lon(i,j,n), lat(i,j,n), &
        lf*cell_area(i,j,n)/1.e6, lake_frac(i,j,n)
    go to 200
202 continue
write (6,*) 'number of Caspian points= ', ktr-1

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
       if (lake_frac(i,j,n) > 0.) write (6,'(a,i6,3f10.2)') "lake already present at ", &
           n, lon(i,j,n), lat(i,j,n), lake_frac(i,j,n)
       lake_frac(i,j,n)= lfrac(l)
       sum_lake(lake_idx(l))= sum_lake(lake_idx(l)) + lfrac(l)*cell_area(i,j,n)
   
       write (10,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i,j,n), lat(i,j,n), &
            lfrac(l)*cell_area(i,j,n)/1.e6, cell_area(i,j,n)/1.e6, lfrac(l)
    enddo
    do l= 1,nlake_def
       write (6,'(i6,3x,a11,2f12.0)') l, lake_name(l), lake_def_area(l), sum_lake(l)/1.e6
    enddo
else
    do l= 1,nlake
       i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
       if (lake_frac(i,j,n) > 0.) write (6,'(a,i6,3f10.2)') "lake already present at ", &
           n, lon(i,j,n), lat(i,j,n), lake_frac(i,j,n)
       lake_frac(i,j,n)= lake_area(l)/(cell_area(i,j,n)/1.e6)
   
       write (10,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i,j,n), lat(i,j,n), &
            lake_area(l), cell_area(i,j,n)/1.e6, lake_area(l)/(cell_area(i,j,n)/1.e6)
       if (lake_area(l) > cell_area(i,j,n)/1.e6) then
           write (6,'(a,2f10.2,2f12.0)') "ERROR: lake area exceeds grid cell area, ", &
             lon(i,j,n), lat(i,j,n), lake_area(l), cell_area(i,j,n)/1.e6
            stop 40
       endif
    enddo
endif

! check to make sure lake does not occur at partial-land cell
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (lake_frac(i,j,n) > 0. .and. land_frac(i,j,n) < 1.) then
             write (6,*) 'ERROR: lake occurs at partial-land cell'
             write (6,*) 'n= ', n, ', i= ', i, ', j= ', j
             stop 18
         endif
      enddo
   enddo
enddo

!  set lake_frac to zero (currently missing value) at partial-land cells
where (land_frac > 0. .and. land_frac < 1.) lake_frac= 0.

! ----------------------------------------------------------------------
! now compute lakes at small-lake cells
! ----------------------------------------------------------------------

! compute total land area, excluding Greenland and Antarctica
area_lakecell_b= 0. ;  area_lakecell_s= 0. ; area_lake_b= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (ga_mask(i,j,n) == 1 .or. land_frac(i,j,n) == 0.) go to 250
         if (lake_frac(i,j,n) > 0.) then
             area_lakecell_b= area_lakecell_b + cell_area(i,j,n)
             area_lake_b= area_lake_b + cell_area(i,j,n)*lake_frac(i,j,n)
         else
             if (scale_area .and. land_frac(i,j,n) /= mval_mdl) then 
                 area_lakecell_s= area_lakecell_s + cell_area(i,j,n)*land_frac(i,j,n)
             else
                 area_lakecell_s= area_lakecell_s + cell_area(i,j,n)
             endif
         endif
250      continue
      enddo
   enddo
enddo
write (6,*) 'area of large lake grid cells= ', area_lakecell_b/1.e+6
write (6,*) 'area of small lake grid cells= ', area_lakecell_s/1.e+6
write (6,*) 'area of large lakes= ', area_lake_b/1.e+6

F_lake_big= area_lake_b/area_lakecell_b
write (6,*) 'F_lake_big= ', F_lake_big

F_lake_small= (F_lake*area_land - F_lake_big*area_lakecell_b)/area_lakecell_s
write (6,*) 'F_lake_small= ', F_lake_small

! compute area-weighted average of runoff ratio over small-lake cells
rr_mn_small= 0. ;  sum= 0. ;  l= 0
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (ga_mask(i,j,n) == 1 .or. land_frac(i,j,n) == 0.) go to 300
         if (rratio(i,j,n) == mval_mdl) write (6,*) &
             "WARNING: land grid cell has no runoff ratio, ", n, i, j
             if (lake_frac(i,j,n) == 0.) then
                 if (scale_area .and. land_frac(i,j,n) /= mval_mdl) then 
                     rr_mn_small= rr_mn_small + rratio(i,j,n)*cell_area(i,j,n)*land_frac(i,j,n)
                     sum= sum + cell_area(i,j,n)*land_frac(i,j,n)
                 else
                     rr_mn_small= rr_mn_small + rratio(i,j,n)*cell_area(i,j,n)
                     sum= sum + cell_area(i,j,n)
                 endif
                 if (rratio(i,j,n) == 0.) l= l + 1
             endif
300      continue
      enddo
   enddo
enddo

write (6,*) 'sum= ', sum/1.e+6
write (6,*) 'rr_mn_small= ', rr_mn_small
rr_mn_small= rr_mn_small/sum
p_const= F_lake_small/rr_mn_small
write (6,*) 'rr_mn_small= ', rr_mn_small
write (6,*) 'k= ', p_const
write (6,*) 'number of small lake cells with runoff ratio equal to zero= ', l

! now compute small lake fractions
area_lake_s= 0.
fmax= -99999. ;  fmin= 99999.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (ga_mask(i,j,n) == 1 .or. land_frac(i,j,n) == 0.) go to 350
         if (rratio(i,j,n) == mval_mdl) write (6,*) &
             "WARNING: land grid cell has no runoff ratio, ", n, i, j
             if (lake_frac(i,j,n) == 0.) then
                 lake_frac(i,j,n)= p_const*rratio(i,j,n)
                 fmin= min(fmin,p_const*rratio(i,j,n))
                 fmax= max(fmax,p_const*rratio(i,j,n))
                 if (scale_area .and. land_frac(i,j,n) /= mval_mdl) then 
                     area_lake_s= area_lake_s + cell_area(i,j,n)*p_const*rratio(i,j,n)*land_frac(i,j,n)
                 else
                     area_lake_s= area_lake_s + cell_area(i,j,n)*p_const*rratio(i,j,n)
                 endif
             endif
350      continue
      enddo
   enddo
enddo
write (6,*) 'area of small lakes= ', area_lake_s/1.e+6
write (6,*) 'fmin= ', fmin, ', fmax= ', fmax

where (lake_frac > 0.) 
   lake_type= 0.
   lake_depth= lake_depth_real
   lake_tau= lake_tau_real
endwhere
where (lake_frac == 0.)
   lake_depth=  0.
   lake_tau=    0.
endwhere

!where (land_frac < 1. .and. land_frac > 0.)
!   lake_depth=  0.
!   lake_tau=    0.
!endwhere




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
   
!   rcode= NF_DEF_VAR (ncid, 'prec', NF_DOUBLE, 2, ndims, varid5)
!   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', 13, 'precipitation')
!   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 4, 'mm/y')
!   rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
!   rcode= NF_DEF_VAR (ncid, 'wroff', NF_DOUBLE, 2, ndims, varid6)
!   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', 23, 'surface runoff of water')
!   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'mm/y')
!   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
   rcode= NF_DEF_VAR (ncid, 'R', NF_DOUBLE, 2, ndims, varid7)
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'long_name', 12, 'runoff ratio')
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid7, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
!  leave define mode
   rcode= NF_ENDDEF (ncid)

!  write coordinate data
   start= 1 ;  count= 1
      
   count(1)= id
   rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

   count(1)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)
   
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lon(:,:,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, lat(:,:,n))

!    lake fraction data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, lake_frac(:,:,n))

!    lake depth data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, lake_depth(:,:,n))

!    lake tau data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, lake_tau(:,:,n))

!   start= 1 ;  count(1)= id ;  count(2)= jd 
!   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, prec(:,:,n))

!   start= 1 ;  count(1)= id ;  count(2)= jd 
!   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, wroff(:,:,n))

   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, rratio(:,:,n))

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area, basin)
deallocate (lake_frac, lake_depth, lake_type, lake_tau)

   
stop

end



