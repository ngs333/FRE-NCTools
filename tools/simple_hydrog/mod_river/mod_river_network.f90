program mod_river_network

! =========================================================================
!   program modifies river-network tocell fields
!   reads and writes lat, lon, tocell, cellarea, and land_frac fields
! =========================================================================

implicit none

integer, parameter :: maxdims= 2
!integer, parameter :: nindmx= 50000
integer, parameter :: ntilmx= 6

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, latid, lonid
integer :: latdim, londim, varid3, varid6, i1, j1, k, ktr, ipts, jpts
integer :: rcode2, nindx, ntiles, longid, latgid, ndm, nindmx
real :: pi, dtr, sum, mval_cella, mval_tocell, mval_landf, reslon
real :: bmin, bmax, b1
logical :: read_casp_network
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname

integer, dimension (maxdims)           :: start, count, dimids, ndims
!integer, dimension (nindmx)            :: idx, jdx, ndx

!real, dimension (nindmx)               :: tocell_new

character(len=100), dimension (ntilmx) :: river_input_file

integer, allocatable, dimension (:)    :: idx, jdx, ndx
real, allocatable, dimension (:)       :: lat_idx, lon_idx, tocell_new
real, allocatable, dimension (:,:)     :: tocell_casp
real, allocatable, dimension (:,:,:)   :: lat, lon, cell_a, tocell, land_fr

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.mod_river_network', form= 'formatted')

read (5,*) nindmx
read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') read_casp_network
close (5)

allocate (idx(nindmx), jdx(nindmx), ndx(nindmx))
allocate (tocell_new(nindmx))


open (21, form= 'formatted')
if (ntiles == 1) then
    nindx= 0 ; ndx= 1
1     nindx= nindx + 1
      read (21,*,end=2) jdx(nindx), idx(nindx), tocell_new(nindx)
      go to 1
2   continue
    close (21)
    nindx= nindx - 1
else
    nindx= 0
11    nindx= nindx + 1
      read (21,*,end=12) ndx(nindx), jdx(nindx), idx(nindx), tocell_new(nindx)
      go to 11
12  continue
    close (21)
    nindx= nindx - 1
endif
write (6,*) 'nindx= ', nindx

do k= 1,nindx
   write (6,'(3i6,f7.0)') ndx(k), jdx(k), idx(k), tocell_new(k)
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

write (6,*) 'jd= ', jd

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
  
write (6,*) 'id= ', id

rcode= nf_close (ncid)


allocate (lat(id,jd,ntiles), lon(id,jd,ntiles))
allocate (cell_a(id,jd,ntiles), land_fr(id,jd,ntiles), tocell(id,jd,ntiles))

! ----------------------------------------------------------------------
!  now get lons and lats from all input river files
! ----------------------------------------------------------------------
do n= 1,ntiles

   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
!   write (10,'(i6,3x,a)') n, trim(river_input_file(n))

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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(1,:,n))
       do i= 2,id
          lat(i,:,:)= lat(1,:,:)
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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,1,n))
       do j= 2,jd
          lon(:,j,:)= lon(:,1,:)
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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(:,:,n))
       
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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(:,:,n))

   endif

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

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_cella= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_cella)
   endif
   write (6,*) 'mval= ', mval_cella

   cell_a(:,:,n)= mval_cella
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_a(:,:,n))

   where (cell_a(:,:,n) == mval_cella) cell_a(:,:,n)= mval_mdl

!   write (10,'(/"cellarea, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f10.1)') (cell_a(i,j,n)/1.e6, i= 1,id)
!   enddo


   write (6,*) 'read land_frac'
   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_landf= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
   endif
   write (6,*) 'mval= ', mval_landf

   land_fr(:,:,n)= mval_landf
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_fr(:,:,n))

   where (land_fr(:,:,n) == mval_landf) land_fr(:,:,n)= mval_mdl

!   write (10,'(/"land_frac, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f6.2)') (land_fr(i,j,n), i= 1,id)
!   enddo


   write (6,*) 'read tocell'
   rcode= nf_inq_varid (ncid, 'tocell', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_tocell= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_tocell)
   endif
   write (6,*) 'mval= ', mval_tocell

   tocell(:,:,n)= mval_tocell
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(:,:,n))

   where (tocell(:,:,n) == mval_tocell) tocell(:,:,n)= mval_mdl

!   write (10,'(/"tocell, tile", i4)') n
!   do j= 1,jd
!      write (10,*) 'j= ', j
!      write (10,'(10f8.1)') (tocell(i,j,n), i= 1,id)
!   enddo

   rcode= nf_close (ncid)
enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jd
      do i= 1,id
         if (cell_a(i,j,n) /= mval_mdl) sum= sum + cell_a(i,j,n)
      enddo
   enddo
   write (6,*) 'sum cell_a (km2)= ', sum/1.e6
enddo


! ----------------------------------------------------------------------
! change tocell values for input j and i indices
! ----------------------------------------------------------------------

write (6,*)
write (6,*) 'tocell changes:'
do k= 1,nindx
   write (6,'(4i6,a,f6.0,a,f6.0)') k, ndx(k), jdx(k), idx(k), ', old= ', &
            tocell(idx(k),jdx(k),ndx(k)), ', new= ', tocell_new(k)
   tocell(idx(k),jdx(k),ndx(k))= tocell_new(k)
enddo


! ----------------------------------------------------------------------
! change caspian drainage so that there is only one value of tocell=0
!   set tocell= 0 at 42.25N, 50.25E (0.5-degree)
! ----------------------------------------------------------------------

allocate (tocell_casp(id,jd))

if ( read_casp_network ) then
! read new caspian network
    open (20, form= 'formatted')
    read (20,'(2i4)') i1, ipts
    read (20,'(2i4)') j1, jpts
    do l= 1,4
       read (20,'(10x)')
    enddo

    tocell_casp= -1.
    do j= j1+jpts-1,j1,-1
       read (20,'(6x,25f8.0)') (tocell_casp(i,j), i= i1,i1+ipts-1)
    enddo
    close (20)

    where (tocell_casp /= -1.)
       tocell(:,:,1)= tocell_casp(:,:)
    endwhere
endif


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'river_network_mod.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))
      
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lat)
   rcode= NF_DEF_DIM (ncid, 'grid_x', id, londim)
   rcode= NF_DEF_DIM (ncid, 'grid_y', jd, latdim)

!  create coordinate variables
   rcode= NF_DEF_VAR (ncid, 'grid_x', NF_DOUBLE, 1, londim, lonid)
   rcode= NF_DEF_VAR (ncid, 'grid_y', NF_DOUBLE, 1, latdim, latid)

!  create attributes for coordinate variables
!    longitude:
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 9, 'longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 12, 'degrees_east')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'cartesian_axis', 1, 'X')

!    latitude:
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'cartesian_axis', 1, 'Y')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 13, 'degrees_north')

!    create data variable and attributes
   ndims(1)= londim ; ndims(2)= latdim
 
   rcode= NF_DEF_VAR (ncid, 'cellarea', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'cell area')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 2, 'm2')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'tocell', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 28, 'direction to downstream cell')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'land_frac', NF_DOUBLE, 2, ndims, varid6)
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', 23, 'land/sea mask(land = 1)')
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 12, 'degrees_east')
 
   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 13, 'degrees_north')
 
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

!    cell area data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, cell_a(:,:,n))

!    tocell data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, tocell(:,:,n))

!    land fraction data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, land_fr(:,:,n))

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

deallocate (lat_idx, lon_idx)
deallocate (lat, lon)
deallocate (cell_a, tocell, land_fr, tocell_casp)
deallocate (idx, jdx, ndx, tocell_new)

   
stop

end



