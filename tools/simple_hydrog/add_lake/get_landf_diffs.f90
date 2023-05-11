program get_landf_diffs

! get diffs in land_fracs from 2 river regrid files

implicit none

include "param.h"  !  get ntiles here

integer, parameter :: maxdims= 3
integer, parameter :: nfile= 2

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

!!TODO:
!!include '/usr/local/include/netcdf.inc'
include 'netcdf.inc'


integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, latid
integer :: lonid, k, rcode2, ndm, kta, nf
real :: pi, dtr, sum, mval_tocell, mval_landf, mval_cella
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ntiles)            :: ktr

character(len=100), dimension (ntiles,nfile) :: river_input_file

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:,:,:) :: lat, lon, tocell, land_frac, cell_area

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.get_landf_diffs', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n,1)
enddo
do n= 1,ntiles
   read (5,'(a)') river_input_file(n,2)
enddo
close (5)

write (6,*) 'fort.5 data read'


! ----------------------------------------------------------------------
!  get lon and lat dims from first river file -- should be identical
!    for all files
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(river_input_file(1,1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open river netcdf file"
    write (6,*) trim(river_input_file(1,1))
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
!   cellarea
! ----------------------------------------------------------------------

allocate (lat(id,jd,ntiles,nfile), lon(id,jd,ntiles,nfile))
allocate (tocell(id,jd,ntiles,nfile), land_frac(id,jd,ntiles,nfile))
allocate (cell_area(id,jd,ntiles,nfile))

do nf = 1,nfile
   do n= 1,ntiles
      write (6,'(i6,3x,a)')  n, trim(river_input_file(n,nf))
      write (10,'(/i6,3x,a)') n, trim(river_input_file(n,nf))

      rcode= NF_OPEN (trim(river_input_file(n,nf)), NF_NOWRITE, ncid)
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
      rcode= nf_get_vara_double (ncid, latid, start, count, lat(:,:,n,nf))


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
      rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,:,n,nf))


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
      rcode= nf_get_vara_double (ncid, varid, start, count, tocell(:,:,n,nf))

      var_units= ' '
      rcode= nf_get_att_text (ncid, varid, "units", var_units)
      write (6,*) 'units= ', var_units

      mval_tocell= 1.e+20
      rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
      if (rcode == 0) then
          rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_tocell)
      endif
      write (6,*) 'mval= ', mval_tocell

      where (tocell(:,:,n,nf) == mval_tocell) tocell(:,:,n,nf)= mval_mdl


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
      rcode= nf_get_vara_double (ncid, varid, start, count, land_frac(:,:,n,nf))

      var_units= ' '
      rcode= nf_get_att_text (ncid, varid, "units", var_units)
      write (6,*) 'units= ', var_units

      mval_landf= 1.e+20
      rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
      if (rcode == 0) then
          rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
      endif
      write (6,*) 'mval= ', mval_landf

      where (land_frac(:,:,n,nf) == mval_landf) land_frac(:,:,n,nf)= mval_mdl


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
      rcode= nf_get_vara_double (ncid, varid, start, count, cell_area(:,:,n,nf))

      var_units= ' '
      rcode= nf_get_att_text (ncid, varid, "units", var_units)
      write (6,*) 'units= ', var_units

      mval_cella= 1.e+20
      rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
      if (rcode == 0) then
          rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_cella)
      endif
      write (6,*) 'mval= ', mval_cella

      where (cell_area(:,:,n,nf) == mval_cella) cell_area(:,:,n,nf)= 0.

      rcode= nf_close (ncid)

   enddo
enddo


ktr= 0 ; kta= 0
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (land_frac(i,j,n,1) /= land_frac(i,j,n,2)) then
             ktr(n)= ktr(n) + 1
             kta= kta + 1
             write (6,'(4i6,2f15.9,2f15.1)') ktr(n), n, j, i, land_frac(i,j,n,1), &
                 land_frac(i,j,n,2), cell_area(i,j,n,1), cell_area(i,j,n,2)
         endif
      enddo
   enddo
enddo
write (6,*) 'total number of affected grid cells= ', kta

open (23, file= 'table.land_frac0_cells', form= 'formatted')
! write i-1 and j-1 for input to ncvarput1 (indices start at 0)
ktr= 0 ; kta= 0
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if ((land_frac(i,j,n,1) /= land_frac(i,j,n,2)) .and. land_frac(i,j,n,1) == 0.) then
             ktr(n)= ktr(n) + 1
             kta= kta + 1
             write (6,'(4i6,2f15.9,2f15.1)') ktr(n), n, j, i, land_frac(i,j,n,1), &
                 land_frac(i,j,n,2), cell_area(i,j,n,1), cell_area(i,j,n,2)
             write (23,'(3i6,f20.15,f15.1)') n, j-1, i-1, land_frac(i,j,n,2), cell_area(i,j,n,2)
         endif
      enddo
   enddo
enddo
close (23)
write (6,*) 'total number of affected grid cells= ', kta


close (10)

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area)


stop

end



