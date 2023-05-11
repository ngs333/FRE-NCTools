program cr_lake_nc_from_idx

use horiz_interp_mod

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: nlkmx= 500000
integer, parameter :: nlake_def= 16, nlake_defp1= nlake_def+1
integer, parameter :: idl= 2880, jdl= 1440
integer, parameter :: idlp1= idl+1, jdlp1= jdl+1
integer, parameter :: iaral= 15, ibalk= 13, ichad= 16, icaspian= nlake_def+1
integer, parameter :: imich= 1, ihuron= 2, ierie= 10, iont= 12

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, l, rcode, varid, ncid, attnum, id, jd, latid
integer :: lonid, latdim, londim, k, ktr, nlake, rcode2, ndm
integer :: j1, fill_val, varid2, idp1, jdp1
real :: pi, dtr, sum, mval_in, lf, scale
character(len=8)   :: var
character(len=40)  :: var_units, latlon_var
character(len=100) :: fname, latlon_nc_file, glcc_file, lname_glcc

character(len=8)   :: vname_glcc= "WaterBod"

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, itlake, lake_idx
integer, dimension (nlkmx)             :: icasp, jcasp, itcasp
integer, dimension (idl,jdl)           :: idat

real, dimension (nlkmx)                :: lake_area, lfrac
real, dimension (idl)                  :: lonl, lonin
real, dimension (jdl)                  :: latl, latin
real, dimension (idlp1)                :: lonlb
real, dimension (jdlp1)                :: latlb
real, dimension (idl,jdl)              :: arlatl, wbdat, interp_mask


!  Lakes (in order) are Michigan, Huron, Superior, Victoria, Tanganyika, Baikal,
!    Great Bear, Malawi, Great Slave, Erie, Winnipeg, Ontario, Balkhash, Ladoga,
!    Aral Sea, Chad, Caspian  ; areas in km2

character(len=11), dimension (nlake_defp1) :: lake_name = &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ', &
   'Caspian    ' /)

real, allocatable, dimension (:)        :: lat_idx, lon_idx, latb, lonb
real, allocatable, dimension (:,:)      :: lake_var, lake_code, arlat
real, allocatable, dimension (:,:)      :: wbd, interp_out, data_out

pi= 4.*atan(1.)
dtr= pi/180.


open (10, file= 'out.cr_lake_nc_from_idx', form= 'formatted')

read (5,'(a)') latlon_var
read (5,'(a)') latlon_nc_file
read (5,'(a)') glcc_file
close (5)

write (6,*) 'input file: ', trim(latlon_nc_file)

! ----------------------------------------------------------------------
! read lat and lon of lake to be inserted in cover field
!   lake_area is the area of lake in the grid cell
! ----------------------------------------------------------------------

write (10,*) 'input lake info:'
open (20, form= 'formatted')
ktr= 0
1 ktr= ktr + 1
   read (20,*,end=2) itlake(ktr), jlake(ktr), ilake(ktr), lfrac(ktr), lake_idx(ktr)
   write (10,'(3i6,f7.1,i6)') itlake(ktr), jlake(ktr), ilake(ktr), lfrac(ktr), lake_idx(ktr)
   go to 1
2 continue
nlake= ktr - 1
write (6,*) 'number of big-lake points to be inserted = ', nlake
write (6,*)
close (20)


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
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= jdl) then
    write (6,*) "ERROR: inconsistent lat dimension, glcc, jdl= ", ndm
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
latlb(jdlp1)= lat2
do j= 2,jdl
   latlb(j)= 0.5*(latl(j)+latl(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= idl) then
    write (6,*) "ERROR: inconsistent lon dimension, glcc, idl= ", ndm
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
arlatl= 0.
do j= 1,jdl
   do i= 1,idl
      arlatl(i,j)= erad*erad*(lonlb(i+1)-lonlb(i))*dtr* &
           (sin(latlb(j+1)*dtr) - sin(latlb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdl
   do i= 1,idl
      sum= sum + arlatl(i,j)
   enddo
enddo
write (6,*) 'sum of glcc area= ', sum

write (10,'(/a)') 'glcc lats'
do j= 1,jdl
   write (10,'(i6,5f10.3)') j, latl(j), latlb(j), latlb(j+1), arlatl(1,j)
enddo

write (10,'(/a)') 'glcc lons'
do i= 1,idl
   write (10,'(i6,3f10.3)') i, lonl(i), lonlb(i), lonlb(i+1)
enddo


wbdat= mval_mdl ;   lname_glcc=''
write (6,*) 'read glcc data: ', trim(vname_glcc)
rcode= nf_inq_varid (ncid, trim(vname_glcc), varid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc variable ",  trim(vname_glcc)
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
!write (6,*) 'fill_val= ', fill_val
   
scale= -99999999.
rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
!write (6,*) 'scale= ', scale

rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', lname_glcc)
!write (6,*) 'long_name= ', trim(lname_glcc)
   
var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
!write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
do j= 1,jdl
   j1= jdlp1-j
!   write (6,*) 'j= ', j, j1
   do i= 1,idl/2
      if (idat(i+idl/2,j1) /= fill_val) then
          wbdat(i,j)= real(idat(i+idl/2,j1))*scale
      endif
   enddo
   do i= idl/2+1,idl
      if (idat(i-idl/2,j1) /= fill_val) then
          wbdat(i,j)= real(idat(i-idl/2,j1))*scale
      endif
   enddo
enddo

rcode= nf_close (ncid)

do j= 1,jdl
   do i= 1,idl
      if (wbdat(i,j) == mval_mdl) then
          write (6,*) "WARNING: wbdat has missing value, ", wbdat(i,j)
      endif
   enddo
enddo



! ----------------------------------------------------------------------
!  get lon and lat dims from first river file -- should be identical
!    for all files
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(latlon_nc_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open river netcdf file"  
    write (6,*) trim(latlon_nc_file)
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
jdp1= jd+1

allocate (lat_idx(jd), latb(jdp1))
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
idp1= id+1

allocate (lon_idx(id), lonb(idp1))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)


! get lon and lat edges
latb(1)= lat1
latb(jdp1)= lat2
do j= 2,jd
   latb(j)= 0.5*(lat_idx(j)+lat_idx(j-1))
enddo

lonb(1)= lon1
lonb(idp1)= lon2
do i= 2,id
   lonb(i)= 0.5*(lon_idx(i)+lon_idx(i-1))
enddo

! compute area of latitude
allocate (arlat(id,jd))

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
write (6,*) 'sum of area (lake file)= ', sum
  

write (10,'(/a)') 'lake lats'
do j= 1,jd
   write (10,'(i6,5f10.3)') j, lat_idx(j), latb(j), latb(j+1), arlat(1,j)
enddo

write (10,'(/a)') 'lake lons'
do i= 1,id
   write (10,'(i6,3f10.3)') i, lon_idx(i), lonb(i), lonb(i+1)
enddo


! read lake-variable grid

allocate (lake_var(id,jd))

start= 1 ;  count= 1

!   write (6,*) 'read lake_var'
rcode= nf_inq_varid (ncid, trim(latlon_var), varid)     ! lake_var field
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find trim(latlon_var)" ; stop 40
endif
rcode= nf_inq_vardimid (ncid, varid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
if (ndm /= id) then
    write (6,*) "ERROR: inconsistent lake_var dimension, ", ndm, id ;  stop 45
endif
rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
if (ndm /= jd) then
    write (6,*) "ERROR: inconsistent lake_var dimension, ", ndm, jd ;  stop 45
endif

start= 1 ;  count(1)= id ;  count(2)= jd
rcode= nf_get_vara_double (ncid, varid, start, count, lake_var)

var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
!write (6,*) 'units= ', var_units

mval_in= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) then
    rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
endif
!write (6,*) 'mval= ', mval_in
   
where (lake_var == mval_in) lake_var= mval_mdl

rcode= nf_close (ncid)
   
       
allocate (lake_code(id,jd))

lake_code= mval_mdl

where (lake_var /= mval_mdl) lake_code= 0.


! ----------------------------------------------------------------------
!  Caspian Sea
! ----------------------------------------------------------------------

lake_code= 0.
write (6,'(/"Caspian lake fractions")')
open (21, form= 'formatted')
ktr= 0
200 ktr= ktr + 1
    read (21,*,end= 202) itcasp(ktr), jcasp(ktr), icasp(ktr), lf
    j= jcasp(ktr) ;  i= icasp(ktr)
    lake_code(i,j)= real(icaspian)
!    write (6,'(4i6,2f10.2,f12.1,f10.3)') ktr, n, i, j, lon_idx(i), lat_idx(j), &
!        lake_code(i,j)
    go to 200
202 continue
write (6,*) 'number of Caspian points= ', ktr-1


! ----------------------------------------------------------------------
!  now insert lakes in lake_code field
!    do not alter Caspian fractions -- already present in frac field
! ----------------------------------------------------------------------

write (10,'(/"input lake fractions")')
write (6,'(/"input lake fractions")')
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l)
   if (lake_code(i,j) > 0.) write (6,'(a,2i5,3f12.4)') "lake already present at ", &
       j, i, lat_idx(j), lon_idx(i), lake_code(i,j)
   lake_code(i,j)= real(lake_idx(l))
   
   write (10,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon_idx(i), lat_idx(j), lfrac(l)
enddo


! set ocean values to missing
where (lake_var == mval_mdl) lake_code= mval_mdl


! ----------------------------------------------------------------------
! interpolate glcc water bodies to lake grid
! ----------------------------------------------------------------------

allocate (wbd(id,jd), interp_out(id,jd), data_out(id,jd))

wbd= mval_mdl

call horiz_interp_init

interp_mask= 1.
where (wbdat == mval_mdl) interp_mask= 0.

!write (6,*) 'WaterBod, tile= ', n
call horiz_interp (wbdat, lonlb*dtr, latlb*dtr, lonb*dtr, &
     latb*dtr, data_out, verbose=0, mask_in=interp_mask, &
     mask_out=interp_out, interp_method="bilinear")
wbd= data_out
where (interp_out == 0.) wbd= mval_mdl
   
deallocate (interp_out, data_out)


! make sure waterbod > 0 at all lake cells
do j= 1,jd
   do i= 1,id
      if (lake_code(i,j) > 0.) then
          if (wbd(i,j) < 0.05) then
              j1= lake_code(i,j)
              write (6,'(a,3x,a,3x,2i6,f6.0,f8.4)') "WARNING: waterbod<0.05 at lake cell, ", &
                   trim(lake_name(j1)), j, i, lake_code(i,j), wbd(i,j)
          endif
          if (wbd(i,j) == 0.) then
              write (6,'(a,3x,a,3x,2i6,f6.0,f8.4)') "ERROR: waterbod=0 at lake cell, ", &
                   trim(lake_name(j1)), j, i, lake_code(i,j), wbd(i,j)
              write (6,*) "exiting..."
              stop
          endif
      endif
   enddo
enddo


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

write (fname, '(a,i1,a)') 'lake_code.nc'
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
rcode= NF_DEF_VAR (ncid, 'lake_code', NF_DOUBLE, 2, ndims, varid)
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'lake_code')
rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
rcode= NF_DEF_VAR (ncid, 'waterbod', NF_DOUBLE, 2, ndims, varid2)
rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 8, 'waterbod')
rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 4, 'none')
rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
!  leave define mode
rcode= NF_ENDDEF (ncid)

!  write coordinate data
start= 1 ;  count= 1
      
count(1)= id
rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

count(1)= jd
rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)
   
!    lake fraction data
start= 1 ;  count(1)= id ;  count(2)= jd 
rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, lake_code)
   
!    waterbod data
start= 1 ;  count(1)= id ;  count(2)= jd 
rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, wbd)


!  close netcdf file
rcode= NF_CLOSE (ncid)

close (10)

deallocate (lat_idx, lon_idx, latb, lonb, arlat)
deallocate (lake_var, lake_code, wbd)

   
stop

end



