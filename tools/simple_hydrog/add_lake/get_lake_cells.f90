program cr_lake_frac_files

! program creates netcdf files of lake fraction on C48 grid 

implicit none

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx= 1000
integer, parameter :: ntiles= 6

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20
real, parameter :: mval_mdl= -9999.

include '/usr/local/include/netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, nlake, rcode2, ndm, latgid, longid
integer :: idm, jdm
real :: pi, dtr, sum, mval_frac, casp_ratio, mval_tocell, mval_landf
real :: mval_cella, mval_laka, mval_basin, ln1
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: long_name, mask_file, fname

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlkmx)             :: ilake, jlake, itlake, ncell_lake, lake_type

real, dimension (nlkmx)                :: lake_area

character(len=100), dimension (ntiles) :: river_input_file, land_month_file

real, allocatable, dimension (:)       :: lat_idx, lon_idx, lat_msk, lon_msk
real, allocatable, dimension (:,:,:)   :: lat, lon, tocell, land_frac, cell_area, basin, &
                                          lake_area_g, area_g, lake_frac

real*4, allocatable, dimension (:,:)     :: adat2

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.cr_lake_frac_files_C48', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(a)') mask_file
close (5)

write (6,*) 'fort.5 data read'


! ----------------------------------------------------------------------
! read lat and lon of lake to be inserted in cover field
!   lake_area is the area of lake in the grid cell
! ----------------------------------------------------------------------

write (10,*) 'input lake info:'
open (20, form= 'formatted')
ktr= 0
1 ktr= ktr + 1
    read (20,*,end=2) itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr)
    write (10,'(3i6,f12.0)') itlake(ktr), jlake(ktr), ilake(ktr), lake_area(ktr)
    go to 1
2 continue
nlake= ktr - 1
write (6,*) 'number of lakes to be inserted = ', nlake
close (20)


! read mask file

rcode= NF_OPEN (trim(mask_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"  
    write (6,*) trim(mask_file)
    stop 1
endif

rcode= nf_inq_varid (ncid, 'latitude', latid)         ! number of lats
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find lat variable (mask)" ; stop 2
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jdm)
write (6,*) 'jdm= ', jdm

allocate (lat_msk(jdm))
start= 1 ;  count= 1 ;  count(1)= jdm
rcode= nf_get_vara_double (ncid, latid, start, count, lat_msk)

       
rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find lon variable (mask)" ; stop 3
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), idm)
write (6,*) 'idm= ', idm

allocate (lon_msk(idm))
start= 1 ;  count(1)= idm
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_msk)
  
rcode= nf_close (ncid)


! ----------------------------------------------------------------------
!  get lon and lat dims from first river file -- should be identical
!    for all files
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(river_input_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"  
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


! ----------------------------------------------------------------------
!  now read lake data -- to grab Caspian lake fractions
! ----------------------------------------------------------------------

allocate (lake_area_g(id,jd,ntiles), area_g(id,jd,ntiles))

do n= 1,ntiles
   write (6,'(i6,3x,a)')  n, trim(land_month_file(n))
   write (10,'(i6,3x,a)') n, trim(land_month_file(n))

   rcode= NF_OPEN (trim(land_month_file(n)), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file"  ; stop 112
   endif

   start= 1 ;  count= 1

   write (6,*) 'read lake_area'
   rcode= nf_inq_varid (ncid, 'lake_area', varid)     ! lake_area field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find lake_area" ; stop 140
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent lake_area dimension, ", ndm, id ;  stop 145
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent lake_area dimension, ", ndm, jd ;  stop 145
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, lake_area_g(:,:,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_laka= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_laka)
   endif
   write (6,*) 'mval= ', mval_laka

   where (lake_area_g(:,:,n) == mval_laka) lake_area_g(:,:,n)= mval_mdl


   write (6,*) 'read area'
   rcode= nf_inq_varid (ncid, 'area', varid)     ! area field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find area" ; stop 150
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent area dimension, ", ndm, id ;  stop 155
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent area dimension, ", ndm, jd ;  stop 155
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, area_g(:,:,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_laka= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_laka)
   endif
   write (6,*) 'mval= ', mval_laka

   where (area_g(:,:,n) == mval_laka) area_g(:,:,n)= mval_mdl

   rcode= nf_close (ncid)
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


sum= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         sum= sum + area_g(i,j,n)
      enddo
   enddo
enddo
write (6,'("sum of area_g= ",f15.1)') sum/1.e+6

!  write lats and lons of input lakes
write (6,'(/"input lake coordinates:")')
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   ln1= lon(i,j,n)
   if (ln1 > 180.) ln1= ln1 - 360.
   write (6,'(2i6,2f10.2,f15.0)') l, n, ln1, lat(i,j,n), basin(i,j,n)
enddo


! ----------------------------------------------------------------------
!  insert lake at every internal drain point
! ----------------------------------------------------------------------

allocate (lake_frac(id,jd,ntiles))

ktr= 0 ;  lake_frac= 0.
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

write (6,*) lake_area_g(10,44,2)/1.e6


! ----------------------------------------------------------------------
!  use model-output lake fraction for Caspian Sea
! ----------------------------------------------------------------------
write (6,'(/"Caspian lake fractions")')
open (21, form= 'formatted')
ktr= 0
200 ktr= ktr + 1
    read (21,*,end= 202) n, j, i
    if (cell_area(i,j,n) == mval_mdl) then
        write (6,*) "ERROR: missing cell_area, ", n, j, i ; stop 200
    endif
    if (lake_area_g(i,j,n) /= mval_mdl) then
        lake_frac(i,j,n)= lake_area_g(i,j,n)/cell_area(i,j,n)
    else
        lake_frac(i,j,n)= 0.
    endif
    write (6,'(4i6,2f10.2,2f12.1,f10.3)') ktr, n, i, j, lon(i,j,n), lat(i,j,n), &
        lake_area_g(i,j,n)/1.e6, cell_area(i,j,n)/1.e6, lake_frac(i,j,n)
    go to 200
202 continue
write (6,*) 'number of Caspian points= ', ktr-1
             

! ----------------------------------------------------------------------
!  now insert lakes in lake_frac field
!    do not alter Caspian fractions -- already present in frac field
! ----------------------------------------------------------------------

!allocate (frac_new(id,jd,ntype))

!write (10,'(/"insert lakes:")')
!frac_new= frac

!write (6,*) 'Caspian drain point, frac= ', frac(21,67,idx_lake)

! multiply all Caspian lake-area by pre-determined factor; do this
!  only for Caspian, the only lake present in the frac field
!do j= 1,jd
!   do i= 1,id
!      if (frac(i,j,idx_lake) > 0.) then
!          frac_new(i,j,idx_lake)= frac(i,j,idx_lake)*casp_ratio
!          if (frac_new(i,j,idx_lake) > 1.) then
!              write (6,'(a,f15.6,2i6,2f12.2)') 'ERROR: Casp, frac > 1, ', &
!                 frac_new(i,j,idx_lake), j, i, lat(j), lon(i)
!          endif
!      endif
!   enddo
!enddo

!write (6,*) 'Caspian drain point, frac_new= ', frac_new(21,67,idx_lake)

write (6,'(/"input lake fractions")')
do l= 1,nlake
   i= ilake(l) ;  j= jlake(l) ;  n= itlake(l)
   if (lake_frac(i,j,n) > 0.) write (6,'(a,i6,3f10.2)') "lake already present at ", &
       n, lon(i,j,n), lat(i,j,n), lake_frac(i,j,n)
   lake_frac(i,j,n)= lake_area(l)/(cell_area(i,j,n)/1.e6)
   
   write (6,'(3i6,2f10.2,2f12.0,f12.3)') l, i, j, lon(i,j,n), lat(i,j,n), &
        lake_area(l), cell_area(i,j,n)/1.e6, lake_area(l)/(cell_area(i,j,n)/1.e6)
   if (lake_area(l) > cell_area(i,j,n)/1.e6) then
       write (6,'(a,2f10.2,2f12.0)') "ERROR: lake area exceeds grid cell area, ", &
         lon(i,j,n), lat(i,j,n), lake_area(l), cell_area(i,j,n)/1.e6
        stop 40
   endif
enddo

!write (6,*) 'Caspian drain point, frac_new= ', frac_new(21,67,idx_lake)


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

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area, basin)
deallocate (lake_area_g, area_g)
deallocate (lake_frac)

   
stop

end



