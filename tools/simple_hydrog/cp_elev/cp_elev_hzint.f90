program cp_elev

   use horiz_interp_mod

   implicit none

   include 'netcdf.inc'

   type (horiz_interp_type) :: Interp

!!include "param.h"  !  get ntiles here

   integer :: ntiles= 6
   integer :: ide= 43200, jde= 21600
   real    :: rese= 1./120.
   integer :: ido= 1440, jdo= 720
   integer :: ncell_cnv= 30
   real    :: reso= 0.25


   integer, parameter :: maxdims= 3

   integer, parameter  :: ni_cells= 3, nj_cells= 3

   integer, parameter  :: npoints= 5
   integer, parameter :: itile1= 1, ilat1= 1, ilon1= 1
   real, parameter  :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
   real, parameter  :: erad= 6.371e+3
   real, parameter  :: mval_mdl= -9999.

   integer :: idep1, jdep1
   integer :: idop1, jdop1
   integer :: npts


   integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
   integer :: latid, lonid, latdim, londim, k, rcode2, ndm, latgid, longid
   integer :: i2, j2, ip, jp, j1, i1, idx, jdx, ii, jj, np, isum, m
   integer :: iwt, ifail, ktr, varid2, varid3, varid4, nxp, nyp
   integer :: idst, jdst, isrc, jsrc, varid5, elev_res
   real :: pi, dtr, sum, mval_landf, sum1, asum1, el_min, el_max, sum2, asum2
   real :: xmean, s2, s3, s4, xmin, xmax, wtsum, mval_elev, interp_res
   character(len=8)   :: var
   character(len=40)  :: var_units, lname_in
   character(len=100) :: fname
   logical :: interp_elev
   logical :: mx, nx, my, ny

   integer, allocatable, dimension (:)       :: start, count, dimids, ndims !!(maxdims)
   integer*2, allocatable, dimension (:,:)   :: elev1 !!ide,jde
   real, allocatable, dimension (:)          :: xdat, wt !!(npts)
   real, allocatable, dimension (:)          :: xpt, ypt !!(npoints)


   real, allocatable, dimension (:,:)    :: out_flow  !!(ni_cells,nj_cells)
!!TODO: Initialization moved below to compile.
!!    (/  8.,   4.,   2., &
!!       16.,   0.,   1., &
!!       32.,  64., 128. /)

   character(len=100), allocatable, dimension (:) :: river_input_file, grid_file
   integer, allocatable, dimension (:,:)  :: kcell
   real, allocatable, dimension (:)       :: lono, lato, lonob, latob
   real, allocatable, dimension (:)       :: lat_idx, lon_idx
   real, allocatable, dimension (:)       :: lone, late, loneb, lateb
   real, allocatable, dimension (:,:)     :: interp_out, interp_mask, data_in, data_out
   real, allocatable, dimension (:,:)     :: xgrid, ygrid, sumi1, sumi2, asumi1, asumi2
   real, allocatable, dimension (:,:)     :: arlatm, elev_chk, elev_in, arlat
   real, allocatable, dimension (:,:)     :: arlato, elv_cnv, elv_sd, elv_min, elv_max
   real, allocatable, dimension (:,:,:)   :: lat, lon, land_frac, elev
   real, allocatable, dimension (:,:,:)   :: elev_sd, elev_min, elev_max
   real, allocatable, dimension (:,:,:)   :: xcrn, ycrn, xcen, ycen
   real, allocatable, dimension (:,:,:)   :: xpoly, ypoly

   character(len=*), parameter :: namelist_file_fp = './cp_elev.nml'

   pi= 4.*atan(1.)
   dtr= pi/180.

   interp_elev= .true.
   interp_res= 1.0

   call read_namelist(namelist_file_fp, ntiles, ide, jde, rese, ido, jdo, ncell_cnv, reso)

   idep1= ide+1
   jdep1= jde+1
   idop1= ido+1
   jdop1= jdo+1
   npts= ncell_cnv*ncell_cnv

   allocate(river_input_file (ntiles), grid_file(ntiles))

   allocate( start(maxdims), count(maxdims), dimids(maxdims), ndims(maxdims) )
   allocate( elev1 (ide,jde))
   allocate( xdat  (npts) , wt (npts) )
   allocate( xpt (npoints) , ypt (npoints)  )
   allocate( out_flow (ni_cells,nj_cells) )

   out_flow = reshape ((/  8.,   4.,   2., 16.,   0.,   1.,  32.,  64., 128. /), shape (out_flow))

   open (10, file= 'out.cp_elev', form= 'formatted')

   read (5,*) elev_res
   do n= 1,ntiles
      read (5,'(a)') river_input_file(n)
   enddo
   read (5,'(l1)') interp_elev

   if ( .not. interp_elev) then
      do n= 1,ntiles
         read (5,'(a)') grid_file(n)
      enddo
   else
      read (5,*) interp_res
   endif
   close (5)

   write (6,*) 'fort.5 data read'
   write (6,*) 'interp_elev= ', interp_elev
   write (6,*) 'interp_res=  ', interp_res


! ----------------------------------------------------------------------
!  read elevation data
! ----------------------------------------------------------------------

   allocate (lone(ide), late(jde), loneb(idep1), lateb(jdep1))
   allocate (elev_in(ide,jde), arlat(ide,jde))

   if (elev_res == 5) then

!   -- ETOPO5 (5-minute, 1/12-degree)

! flip latitudes on read
      open (37, form= 'unformatted', access= 'direct', recl=8640)
      write (10,*) 'north pole elevations'
      do j= 1,jde
         j1= jde - j + 1
         read (37,rec=j) (elev1(i,j1), i= 1,ide)
         if (j == 1) write (10,*) j, j1
         if (j == 1) write (10,*) (elev1(i,j1), i= 1,ide)
      enddo
      close (37)

      write (10,*) 'north pole elevations'
      write (10,*) j, j1
      write (10,*) (elev1(i,j1), i= 1,ide)

      elev_in= real(elev1)

! compute lons and lats for elevation data set
      do j= 1,jdep1
         lateb(j)= lat1 + real(j-1)*rese
      enddo

      do j= 1,jde
         late(j)= 0.5*(lateb(j)+lateb(j+1))
      enddo

      do i= 1,idep1
         loneb(i)= lon1 + real(i-1)*rese
      enddo

      do i= 1,ide
         lone(i)= 0.5*(loneb(i)+loneb(i+1))
      enddo

      el_min=  9999999.
      el_max= -9999999.
      do j= 1,jde
         do i= 1,ide
            el_min= min(el_min,elev_in(i,j))
            el_max= max(el_max,elev_in(i,j))
         enddo
      enddo
      write (6,*) 'el_min= ', el_min, ', el_max= ', el_max

! set negative elevations in Caspian Sea area to 0
      do j= 1,jde
         if (late(j) > 36. .and. late(j) < 48.) then
            do i= 1,ide
               if (lone(i) > 45. .and. lone(i) < 56.) then
                  if (elev_in(i,j) < 0.) elev_in(i,j)= 0.
               endif
            enddo
         endif
      enddo

! set other negative elevations to missing value
      where (elev_in < 0.) elev_in= mval_mdl

! write input elevation data as a check
      write (fname, '(a)') 'elev_data_input.nc'
      rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
      rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

!  create dimensions (grid_x, grid_y)
      rcode= NF_DEF_DIM (ncid, 'lon',  ide,   londim)
      rcode= NF_DEF_DIM (ncid, 'lat',  jde,   latdim)

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
      rcode= NF_DEF_VAR (ncid, 'elev', NF_DOUBLE, 2, ndims, varid)
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

!  leave define mode
      rcode= NF_ENDDEF (ncid)

!  write coordinate data
      start= 1 ;  count= 1 ;  count(1)= ide
      rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lone)

      start= 1 ;  count= 1 ;  count(1)= jde
      rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, late)

!    elevation data
      start= 1 ;  count(1)= ide ;  count(2)= jde
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, elev_in)

!  close netcdf file
      rcode= NF_CLOSE (ncid)

   else if (elev_res == 30) then

!   -- GTOPO30 (30 arc-second, 1/120 degree)

      rcode= NF_OPEN ('gtopo30.nc', NF_NOWRITE, ncid)
      if  (rcode /= 0) then
         write (6,*) "ERROR: cannot open netcdf file"  ; stop
      endif

! get latitudes and longitudes of elevation grid

      rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
      rcode= nf_inq_vardimid (ncid, latid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), n)
      if (n /= jde) then
         write (6,*) "ERROR: inconsistent jde, ", n, jde
         stop 20
      endif

      start= 1 ;  count= 1 ;  count(1)= jde
      rcode= nf_get_vara_double (ncid, latid, start, count, late)

      rcode= nf_inq_varid (ncid, 'latb', latid)        ! number of lat edges
      rcode= nf_inq_vardimid (ncid, latid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), n)
      if (n /= jdep1) then
         write (6,*) "ERROR: inconsistent jdep1, ", n, jdep1
         stop 20
      endif

      start= 1 ;  count= 1 ;  count(1)= jdep1
      rcode= nf_get_vara_double (ncid, latid, start, count, lateb)

      rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
      rcode= nf_inq_vardimid (ncid, lonid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), n)
      if (n /= ide) then
         write (6,*) "ERROR: inconsistent ide, ", n, ide
         stop 20
      endif

      start= 1 ;  count= 1 ;  count(1)= ide
      rcode= nf_get_vara_double (ncid, lonid, start, count, lone)

      rcode= nf_inq_varid (ncid, 'lonb', lonid)        ! number of lon edges
      rcode= nf_inq_vardimid (ncid, lonid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), n)
      if (n /= idep1) then
         write (6,*) "ERROR: inconsistent idep1, ", n, idep1
         stop 20
      endif

      start= 1 ;  count= 1 ;  count(1)= idep1
      rcode= nf_get_vara_double (ncid, lonid, start, count, loneb)

!  read elevation data and attributes

      mval_elev= -999999.
      rcode= nf_inq_varid (ncid, 'elev', varid)
      rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
      if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_elev)
      write (6,*) 'mval= ', mval_elev

      rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
      if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', lname_in)
      write (6,*) 'long name= ', trim(lname_in)

      rcode= nf_inq_attid (ncid, varid, 'units', attnum)
      if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'units', var_units)
      write (6,*) 'units= ', trim(var_units)

      rcode= nf_inq_vardimid (ncid, varid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), n)
      rcode= nf_inq_dimlen (ncid, dimids(2), i)
      if (n /= ide .or. i /= jde) then
         write (6,*) "ERROR: inconsistent elev dims, ", n, ide, i, jde
         stop 25
      endif

      start= 1 ; count= 1
      count(1)= ide  ;  count(2)= jde
      rcode= nf_get_vara_double (ncid, varid, start, count, elev_in)

! close netcdf file
      rcode= nf_close (ncid)

      el_min=  9999999.
      el_max= -9999999.
      do j= 1,jde
         do i= 1,ide
            if (elev_in(i,j) /= mval_mdl) then
               el_min= min(el_min,elev_in(i,j))
               el_max= max(el_max,elev_in(i,j))
            endif
         enddo
      enddo
      write (6,*) 'el_min= ', el_min, ', el_max= ', el_max

! set all negative elevations to 0 (gtopo30 provides elevation only over land)
      where (elev_in /= mval_mdl .and. elev_in < 0.) elev_in= 0.
! set all missing elevations (water) to 0   12 mar 2014
! (some land points near ob/yenisei mouth have undefined elevation on c384)
      where (elev_in == mval_mdl) elev_in= 0.
   else
      write (6,*) "ERROR: invalid elevation resolution ", elev_res
      stop 30
   endif


! compute area
   arlat= 0.
   do j= 1,jde
      do i= 1,ide
         arlat(i,j)= erad*erad*(loneb(i+1)-loneb(i))*dtr* &
            (sin(lateb(j+1)*dtr) - sin(lateb(j)*dtr))
      enddo
   enddo

   sum= 0.
   do j= 1,jde
      do i= 1,ide
         sum= sum + arlat(i,j)
      enddo
   enddo
   write (6,*) 'sum of elev area= ', sum

   write (10,'(/a)') 'elev lats'
   do j= 1,jde
      write (10,'(i6,5f10.3)') j, late(j), lateb(j), lateb(j+1), arlat(1,j)
   enddo

   write (10,'(/a)') 'elev lons'
   do i= 1,ide
      write (10,'(i6,3f10.3)') i, lone(i), loneb(i), loneb(i+1)
   enddo



! ----------------------------------------------------------------------
! try averaging high-res data to 1-degree
! ----------------------------------------------------------------------
   if ( interp_elev ) then
      allocate (lato(jdo), lono(ido), latob(jdop1), lonob(idop1))
      allocate (arlato(ido,jdo), elv_cnv(ido,jdo), elv_sd(ido,jdo))
      allocate (elv_min(ido,jdo), elv_max(ido,jdo))

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

      write (10,'(/a)') 'lats for manual interpolation of input elevation'
      do j= 1,jdo
         write (10,'(i6,5f10.3)') j, lato(j), latob(j), latob(j+1), arlato(1,j)
      enddo

      write (10,'(/a)') 'lons for manual interpolation of input elevation'
      do i= 1,ido
         write (10,'(i6,3f10.3)') i, lono(i), lonob(i), lonob(i+1)
      enddo

      sum= 0.
      do j= 1,jde
         do i= 1,ide
            if (elev_in(i,j) >= 0.) sum= sum+arlat(i,j)
         enddo
      enddo
      write (6,*) 'area, elev>=0 = ', sum

!  use a manual interpolation of elevation data; compute max, min, std deviation
      elv_cnv= mval_mdl  ;  elv_sd= mval_mdl  ; sum= 0.
      elv_min= 99999999. ;  elv_max= -99999999.
      do j= 1,jdo
         j1= (j-1)*ncell_cnv+1
         j2= j*ncell_cnv
         do i= 1,ido
            i1= (i-1)*ncell_cnv+1
            i2= i*ncell_cnv
            asum1= 0. ;  sum1= 0. ;  ktr= 0 ;  asum2= 0.
            do jj= j1,j2
               do ii= i1,i2
                  if (elev_in(ii,jj) /= mval_mdl) then
                     ktr= ktr + 1
                     sum= sum + arlat(ii,jj)
                     asum1= asum1 + arlat(ii,jj)
                     asum2= asum2 + arlat(ii,jj)**2.
                     sum1= sum1 + elev_in(ii,jj)*arlat(ii,jj)
                     elv_min(i,j)= min(elv_min(i,j),elev_in(ii,jj))
                     elv_max(i,j)= max(elv_max(i,j),elev_in(ii,jj))
                     xdat(ktr)= elev_in(ii,jj)
                     wt(ktr) = arlat(ii,jj)
                  endif
               enddo
            enddo
!          if (ktr /= npts) then
!              write (6,*) "ERROR: wrong number of points, 1-deg"
!              write (6,*) npts, ktr
!              stop 10
!          endif
            if (ktr < 3) go to 15

            sum2= 0.
            do jj= j1,j2
               do ii= i1,i2
                  if (elev_in(ii,jj) /= mval_mdl) then
                     sum2= sum2 + ((elev_in(ii,jj)-(sum1/asum1))**2.)*arlat(ii,jj)
                  endif
               enddo
            enddo

            elv_cnv(i,j)= sum1/asum1
            elv_sd(i,j)= (sum2/(asum1-asum2/asum1))**0.5
            if (ktr == npts .and. abs(asum1 - arlato(i,j)) > 1.0e-8) then
               write (6,*) "ERROR: elev 1-deg areas do not agree, ", asum1, arlato(i,j)
               stop 17
            endif

            iwt= 1 ;  ifail= 0
!          call G01AAF (ktr, xdat, iwt, wt, xmean, s2, s3, s4, xmin, xmax, wtsum, ifail)

!          if (abs(xmean-elv_cnv(i,j)) > 1.e-10) then
!              write (6,*) "ERROR: manual/NAG values disagree, mean"
!              write (6,*) elv_cnv(i,j), xmean
!              stop 21
!          endif
!          if (abs(s2-elv_sd(i,j)) > 1.e-11) then
!              write (6,*) "ERROR: manual/NAG values disagree, sd"
!              write (6,*) elv_sd(i,j), s2
!              stop 22
!          endif
!          if (xmin /= elv_min(i,j)) then
!              write (6,*) "ERROR: manual/NAG values disagree, min"
!              write (6,*) elv_min(i,j), xmin
!              stop 23
!          endif
!          if (xmax /= elv_max(i,j)) then
!              write (6,*) "ERROR: manual/NAG values disagree, max"
!              write (6,*) elv_max(i,j), xmax
!              stop 24
!          endif

15          continue

         enddo
      enddo
      write (6,*) 'sum of conversion area= ', sum

      where (elv_cnv == mval_mdl)
         elv_sd= mval_mdl
         elv_min= mval_mdl
         elv_max= mval_mdl
      endwhere
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
!write (6,*) 'jd= ', jd
   jdp1= jd + 1

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
!write (6,*) 'id= ', id
   idp1= id + 1

   allocate (lon_idx(id))
   start= 1 ;  count(1)= id
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)

   rcode= nf_close (ncid)



! ----------------------------------------------------------------------
! now open river files -- read lat,lon grids, land_frac
! ----------------------------------------------------------------------

   allocate (lat(id,jd,ntiles), lon(id,jd,ntiles))
   allocate (land_frac(id,jd,ntiles))

   do n= 1,ntiles
!   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
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


!   write (6,*) 'read land_frac'
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
!   write (6,*) 'units= ', var_units

      mval_landf= 1.e+20
      rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
      if (rcode == 0) then
         rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
      endif
!   write (6,*) 'mval= ', mval_landf

      where (land_frac(:,:,n) == mval_landf) land_frac(:,:,n)= mval_mdl


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

   enddo


   if (.not. interp_elev) then

! get lat and lon data from grid file

      rcode= NF_OPEN (trim(grid_file(1)), NF_NOWRITE, ncid)
      if (rcode /= 0) then
         write (6,*) "ERROR: cannot open grid netcdf file"
         write (6,*) trim(grid_file(1))
         stop 100
      endif

      rcode= nf_inq_varid (ncid, 'y', latid)         ! number of lats
      if (rcode /= 0) then
         write (6,*) "ERROR: cannot find y" ; stop 20
      endif
      rcode= nf_inq_vardimid (ncid, latid, dimids)
      rcode= nf_inq_dimlen (ncid, dimids(1), nxp)
      rcode= nf_inq_dimlen (ncid, dimids(2), nyp)
      write (6,*) 'nxp= ', nxp, ', nyp= ', nyp

      rcode= nf_close (ncid)

      if (nxp /= id*2+1 .or. nyp /= jd*2+1) then
         write (6,*) "ERROR: id/jd and nxp/nyp disagree"
         write (6,*) id, jd, nxp, nyp
         stop 101
      endif


      allocate (ygrid(nxp,nyp), xgrid(nxp,nyp))
      allocate (xcrn(idp1,jdp1,ntiles), ycrn(idp1,jdp1,ntiles))
      allocate (xcen(id,jd,ntiles), ycen(id,jd,ntiles))
      allocate (xpoly(npoints,id,jd), ypoly(npoints,id,jd))

      do n= 1,ntiles

         write (6,'(i6,3x,a)')  n, trim(grid_file(n))
         write (10,'(/i6,3x,a)') n, trim(grid_file(n))

         rcode= NF_OPEN (trim(grid_file(n)), NF_NOWRITE, ncid)
         if (rcode /= 0) then
            write (6,*) "ERROR: cannot open grid netcdf file"
            write (6,*) trim(grid_file(1))
            stop 200
         endif

         start= 1 ;  count= 1 ;  count(1)= nxp ;  count(2)= nyp
         rcode= nf_get_vara_double (ncid, latid, start, count, ygrid)

         rcode= nf_inq_varid (ncid, 'x', lonid)         ! number of lons
         if (rcode /= 0) then
            write (6,*) "ERROR: cannot find x" ; stop 30
         endif
         rcode= nf_inq_vardimid (ncid, lonid, dimids)
         rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
         if (ndm /= nxp) then
            write (6,*) "ERROR: inconsistent lon dimension, ", ndm, nxp ;  stop 250
         endif
         rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
         if (ndm /= nyp) then
            write (6,*) "ERROR: inconsistent lat dimension, ", ndm, nyp ;  stop 250
         endif

         start= 1 ;  count(1)= nxp ;  count(2)= nyp
         rcode= nf_get_vara_double (ncid, lonid, start, count, xgrid)

         do j= 1,jdp1
            do i= 1,idp1
               ycrn(i,j,n)= ygrid(2*i-1,2*j-1)
               xcrn(i,j,n)= xgrid(2*i-1,2*j-1)
            enddo
         enddo

         do j= 1,jd
            do i= 1,id
               ycen(i,j,n)= ygrid(2*i,2*j)
               xcen(i,j,n)= xgrid(2*i,2*j)
            enddo
         enddo

!       write (10,'(/"grid lat centers, tile", i4)') n
!       do j= 1,jd
!          write (10,*) 'j= ', j
!          write (10,'(10f10.4)') (ycen(i,j,n), i= 1,id)
!       enddo

!       write (10,'(/"grid lon centers, tile", i4)') n
!       do j= 1,jd
!          write (10,*) 'j= ', j
!          write (10,'(10f10.4)') (xcen(i,j,n), i= 1,id)
!       enddo

!       write (10,'(/"grid lat corners, tile", i4)') n
!       do j= 1,jdp1
!          write (10,*) 'j= ', j
!          write (10,'(10f10.4)') (ycrn(i,j,n), i= 1,idp1)
!       enddo

!       write (10,'(/"grid lon corners, tile", i4)') n
!       do j= 1,jdp1
!          write (10,*) 'j= ', j
!          write (10,'(10f10.4)') (xcrn(i,j,n), i= 1,idp1)
!       enddo

         do j= 1,jd
            do i= 1,id
               if (ycen(i,j,n) /= lat(i,j,n)) then
                  write (6,*) "ERROR: lats differ river/grid"
                  write (6,*) j, i, lat(i,j,n), ycen(i,j,n)
                  stop 300
               endif
               if (xcen(i,j,n) /= lon(i,j,n)) then
                  write (6,*) "ERROR: lons differ river/grid"
                  write (6,*) j, i, lon(i,j,n), xcen(i,j,n)
                  stop 301
               endif
            enddo
         enddo

         do j= 1,jd
            do i= 1,id
               xpoly(1,i,j)= xcrn(i,j,n)
               xpoly(2,i,j)= xcrn(i,j+1,n)
               xpoly(3,i,j)= xcrn(i+1,j+1,n)
               xpoly(4,i,j)= xcrn(i+1,j,n)
               xpoly(5,i,j)= xcrn(i,j,n)

               ypoly(1,i,j)= ycrn(i,j,n)
               ypoly(2,i,j)= ycrn(i,j+1,n)
               ypoly(3,i,j)= ycrn(i+1,j+1,n)
               ypoly(4,i,j)= ycrn(i+1,j,n)
               ypoly(5,i,j)= ycrn(i,j,n)
            enddo
         enddo

!       write (10,'(/"grid cell polygons, tile", i4)') n
!       do j= 1,jd
!          do i= 1,id
!             write (10,*) 'j= ', j, ', i= ', i
!             write (10,'(a,5f12.3)') '  xpoly= ', (xpoly(l,i,j), l= 1,npoints)
!             write (10,'(a,5f12.3)') '  ypoly= ', (ypoly(l,i,j), l= 1,npoints)
!          enddo
!       enddo

         rcode= nf_close (ncid)
      enddo

   endif

! ----------------------------------------------------------------------
! interpolate high-res elevation data to model grid
! ----------------------------------------------------------------------

   allocate (elev(id,jd,ntiles), elev_sd(id,jd,ntiles))
   allocate (elev_min(id,jd,ntiles), elev_max(id,jd,ntiles))

   elev= mval_mdl ;  elev_sd= mval_mdl

   call horiz_interp_init

   if (interp_elev) then
      elev_min= mval_mdl ;  elev_max= mval_mdl

      allocate (interp_mask(ido,jdo), data_in(ido,jdo))
      allocate (interp_out(id,jd), data_out(id,jd))

      interp_mask= 1.
      where (elv_cnv == mval_mdl) interp_mask= 0.

      do n= 1,ntiles
!       write (6,'(/a,i4)') 'interpolation, tile= ', n
         data_in= elv_cnv
         call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(:,:,n)*dtr, &
            lat(:,:,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         elev(:,:,n)= data_out
         where (interp_out == 0.) elev(:,:,n)= mval_mdl
      enddo

      interp_mask= 1.
      where (elv_sd == mval_mdl) interp_mask= 0.
      do n= 1,ntiles
!       write (6,'(/a,i4)') 'interpolation, tile= ', n
         data_in= elv_sd
         call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(:,:,n)*dtr, &
            lat(:,:,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         elev_sd(:,:,n)= data_out
         where (interp_out == 0.) elev_sd(:,:,n)= mval_mdl
      enddo

      interp_mask= 1.
      where (elv_min == mval_mdl) interp_mask= 0.
      do n= 1,ntiles
!       write (6,'(/a,i4)') 'interpolation, tile= ', n
         data_in= elv_min
         call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(:,:,n)*dtr, &
            lat(:,:,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         elev_min(:,:,n)= data_out
         where (interp_out == 0.) elev_min(:,:,n)= mval_mdl
      enddo

      interp_mask= 1.
      where (elv_max == mval_mdl) interp_mask= 0.
      do n= 1,ntiles
!       write (6,'(/a,i4)') 'interpolation, tile= ', n
         data_in= elv_max
         call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(:,:,n)*dtr, &
            lat(:,:,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         elev_max(:,:,n)= data_out
         where (interp_out == 0.) elev_max(:,:,n)= mval_mdl
      enddo

      deallocate (data_in)

   else

      allocate (interp_mask(ide,jde))
      allocate (interp_out(id,jd), data_out(id,jd))

      elev_min= 99999999. ;  elev_max= -99999999.

!  if just one tile, then compute elev, elev_min, elev_max, elev_sd fields manually
!    elev and model data are ordered south to north from -90 and east to west from 0
      if (ntiles == 1) then

! compute area of latitude
         allocate (arlatm(id,jd), elev_chk(id,jd))
         arlatm= 0.
         do j= 1,jd
            do i= 1,id
               arlatm(i,j)= erad*erad*(xcrn(i+1,ilat1,itile1)-xcrn(i,ilat1,itile1))*dtr* &
                  (sin(ycrn(ilon1,j+1,itile1)*dtr) - sin(ycrn(ilon1,j,  itile1)*dtr))
            enddo
         enddo

         sum= 0.
         do j= 1,jd
            do i= 1,id
               sum= sum + arlatm(i,j)
            enddo
         enddo
         write (6,*) 'sum of elev area (model)= ', sum

         interp_mask= 1.
         where (elev_in == mval_mdl) interp_mask= 0.

         call horiz_interp_init
         call horiz_interp_new (Interp, loneb*dtr, lateb*dtr, xcrn(:,:,itile1)*dtr, &
            ycrn(:,:,itile1)*dtr, verbose=2,  &
            interp_method="conservative")

         write (6,'(/a,i4)') 'interpolation, tile= 1'
         write (10,'(/a,i4)') 'interpolation, tile= 1'

         call horiz_interp (Interp, elev_in, data_out, mask_in=interp_mask, mask_out=interp_out)
         write (6,*) 'in main, nxgrid= ', Interp%nxgrid

!        call horiz_interp_init
!        call horiz_interp_new (Interp, loneb*dtr, lateb*dtr, xcrn(:,:,itile1)*dtr,  &
!                ycrn(:,:,itile1)*dtr, verbose=2, interp_method="conservative")
!        call horiz_interp (Interp, elev_in, data_out, verbose=2, mask_in=interp_mask, mask_out=interp_out)

         elev_chk= data_out
         where (interp_out == 0.) elev_chk= mval_mdl

         elev= mval_mdl  ;  elev_sd= mval_mdl  ; sum= 0.
         do j= 1,jd
            do i= 1,id
               asum1= 0. ;  sum1= 0. ;  ktr= 0 ;  asum2= 0.
               do jj= 1,jde
                  if (late(jj) <= ycrn(ilon1,j,  itile1)) go to 110
                  if (late(jj) >  ycrn(ilon1,j+1,itile1)) go to 115
                  do ii= 1,ide
                     if (lone(ii) <= xcrn(i,  ilat1,itile1)) go to 105
                     if (lone(ii) >  xcrn(i+1,ilat1,itile1)) go to 110
                     if (elev_in(ii,jj) /= mval_mdl) then
                        ktr= ktr + 1
!                            if (i == 33 .and. j == 52) then
!                                write (6,'(5i6,f12.5)') ktr, i, j, jj, ii, elev_in(ii,jj)
!                            endif
                        sum= sum + arlat(ii,jj)
                        asum1= asum1 + arlat(ii,jj)
                        asum2= asum2 + arlat(ii,jj)**2.
                        sum1= sum1 + elev_in(ii,jj)*arlat(ii,jj)
                        elev_min(i,j,itile1)= min(elev_min(i,j,itile1),elev_in(ii,jj))
                        elev_max(i,j,itile1)= max(elev_max(i,j,itile1),elev_in(ii,jj))
                        xdat(ktr)= elev_in(ii,jj)
                        wt(ktr) = arlat(ii,jj)
                     endif
105                  continue
                  enddo           ! end of ii loop
110               continue
               enddo              ! end of jj loop
115            continue
               if (ktr < 3) then
                  elev_chk(i,j)= mval_mdl
                  if (ktr > 0) write (6,'(a,3i6,2f12.5)') 'ktr= ', ktr, j, i, &
                     xcrn(i,ilat1,itile1), ycrn(ilon1,j,itile1)
                  go to 120
               endif

               sum2= 0.
               do k= 1,ktr
                  sum2= sum2 + ((xdat(k)-(sum1/asum1))**2.)*wt(k)
               enddo

               elev(i,j,itile1)= sum1/asum1
               elev_sd(i,j,itile1)= (sum2/(asum1-asum2/asum1))**0.5

!               if (i == 33 .and. j == 52) then
!                   write (6,*) elev(i,j,itile1), elev_sd(i,j,itile1)
!               endif
!!              if (abs(asum1 - arlatm(i,j)) > 1.0e-10) then
!!                  write (6,*) "ERROR: elev model areas do not agree, ", asum1, arlatm(i,j), i, j
!!                  stop 117
!!              endif

120            continue

            enddo                 ! end of i loop
         enddo                    ! end of j loop
         write (6,*) 'sum of conversion area= ', sum

         where (elev == mval_mdl)
            elev_sd= mval_mdl
            elev_min= mval_mdl
            elev_max= mval_mdl
         endwhere

         deallocate (arlatm)

      else

         allocate (sumi1(id,jd), sumi2(id,jd), asumi1(id,jd), asumi2(id,jd))
         allocate (kcell(ide,jde))

         interp_mask= 1.
         where (elev_in == mval_mdl) interp_mask= 0.

         write (10,*) "Interp data"
         kcell= 0

         do n= 1,ntiles
            call horiz_interp_init
            call horiz_interp_new (Interp, loneb*dtr, lateb*dtr, xcrn(:,:,n)*dtr, ycrn(:,:,n)*dtr, &
               verbose=2, mask_in=interp_mask, mask_out=interp_out, interp_method="conservative")

            write (6,'(/a,i4)') 'interpolation, tile= ', n
            write (10,'(/a,i4)') 'interpolation, tile= ', n
!           call horiz_interp (elev_in, loneb*dtr, lateb*dtr, xcrn(:,:,n)*dtr, &
!                ycrn(:,:,n)*dtr, data_out, verbose=2, mask_in=interp_mask, &
!                mask_out=interp_out, interp_method="conservative")

!           call horiz_interp (Interp, elev_in, data_out)
            write (6,*) 'in main, nxgrid= ', Interp%nxgrid
!           elev(:,:,n)= data_out
!           where (interp_out == 0.) elev(:,:,n)= mval_mdl

            sumi1= 0. ;  sumi2= 0. ;  asumi1= 0. ;  asumi2= 0.
            do i = 1, Interp%nxgrid
               idst= Interp%i_dst(i) ;  jdst= Interp%j_dst(i)
               isrc= Interp%i_src(i) ;  jsrc= Interp%j_src(i)
               kcell(isrc,jsrc)= kcell(isrc,jsrc)+1

               if (elev_in(isrc,jsrc) == mval_mdl) then
                  write (6,*) "ERROR: missing input elevation"
                  write (6,*) "  i= ", isrc, ", j= ", jsrc
                  stop 100
               endif

               asumi1(idst,jdst)= asumi1(idst,jdst) + Interp%area_frac_dst(i)
               asumi2(idst,jdst)= asumi2(idst,jdst) + Interp%area_frac_dst(i)**2.
               sumi1(idst,jdst)=  sumi1(idst,jdst)  + Interp%area_frac_dst(i)*elev_in(isrc,jsrc)

               elev_min(idst,jdst,n)= min(elev_min(idst,jdst,n),elev_in(isrc,jsrc))
               elev_max(idst,jdst,n)= max(elev_max(idst,jdst,n),elev_in(isrc,jsrc))

               if (Interp%i_src(i) == 1423 .and. Interp%j_src(i) == 1389) &
                  write (10,'(a,7i8,f20.10)') 'chk point: ', i, kcell(Interp%i_src(i),Interp%j_src(i)), &
                  Interp%i_src(i), Interp%j_src(i), Interp%i_dst(i), Interp%j_dst(i), n, &
                  Interp%area_frac_dst(i)
            enddo

            do i= 1, Interp%nxgrid
               idst= Interp%i_dst(i) ;  jdst= Interp%j_dst(i)
               isrc= Interp%i_src(i) ;  jsrc= Interp%j_src(i)
               sumi2(idst,jdst)= sumi2(idst,jdst) + ((elev_in(isrc,jsrc)- &
                  sumi1(idst,jdst)/asumi1(idst,jdst))**2.)*Interp%area_frac_dst(i)
            enddo

            write (10,*) "CHECK grid point data"
            do i= 1, Interp%nxgrid
               idst= Interp%i_dst(i) ;  jdst= Interp%j_dst(i)
               isrc= Interp%i_src(i) ;  jsrc= Interp%j_src(i)
               elev(idst,jdst,n)= sumi1(idst,jdst)
               if ((asumi1(idst,jdst)-asumi2(idst,jdst)/asumi1(idst,jdst)) /= 0.) then
                  elev_sd(idst,jdst,n)= (sumi2(idst,jdst)/ &
                     (asumi1(idst,jdst)-asumi2(idst,jdst)/asumi1(idst,jdst)))**0.5
               endif
!              if (idst == 25 .and. jdst == 48 .and. n == 1) then
!              if (idst == 7 .and. jdst == 10 .and. n == 1) then
               if (idst == 1 .and. jdst == 5 .and. n == 1) then
                  write (10,'(3i8,13e20.8)') i, isrc, jsrc, elev_in(isrc,jsrc), arlat(isrc,jsrc), &
                     interp_mask(isrc,jsrc), Interp%area_frac_dst(i), sumi1(idst,jdst), sumi2(idst,jdst), &
                     asumi1(idst,jdst), asumi2(idst,jdst), elev(idst,jdst,n), elev_sd(idst,jdst,n), &
                     elev_min(idst,jdst,n), elev_max(idst,jdst,n), interp_out(idst,jdst)
               endif
            enddo

!    for consistency, where sd can't be computed, set elev, max, and min to missing value
            where (elev_sd(:,:,n) == mval_mdl)
               elev(:,:,n)= mval_mdl
               elev_min(:,:,n)= mval_mdl
               elev_max(:,:,n)= mval_mdl
            endwhere

!    use output interp mask to screen data
            where (interp_out == 0.)
               elev(:,:,n)= mval_mdl
               elev_sd(:,:,n)= mval_mdl
               elev_min(:,:,n)= mval_mdl
               elev_max(:,:,n)= mval_mdl
            endwhere
            call horiz_interp_del (Interp)

         enddo      ! end of ntile loop

!        write (10,*) "kcell > 1"
!        ktr= 0
!        do j= 1,jde
!           do i= 1,ide
!              if (kcell(i,j) > 1) then
!                  ktr= ktr + 1
!                  write (10,'(4i8,2f12.3)') ktr, kcell(i,j), j, i, late(j), lone(i)
!              endif
!           enddo
!        enddo

         deallocate (sumi1, sumi2, asumi1, asumi2)
         deallocate (kcell)

      endif

   endif

   deallocate (interp_out, interp_mask)
   deallocate (data_out)
   deallocate (late, lateb, lone, loneb, elev_in, arlat)

! set ocean values to missing
   where (land_frac == 0.)
      elev= mval_mdl
      elev_sd= mval_mdl
      elev_min= mval_mdl
      elev_max= mval_mdl
   endwhere

   if (ntiles == 1 .and. (.not. interp_elev)) then
      do j= 1,jd
         do i= 1,id
            if (elev(i,j,itile1) == mval_mdl) elev_chk(i,j)= mval_mdl
         enddo
      enddo
   endif

! set part-land values to zero
!where (land_frac > 0 .and. land_frac < 1.)
!   elev= 0.
!endwhere



! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

   do n= 1,ntiles
      write (fname, '(a,i1,a)') 'elev.tile', n, '.nc'
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
      rcode= NF_DEF_VAR (ncid, 'elev', NF_DOUBLE, 2, ndims, varid)
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'mean elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      if (ntiles == 1 .and. (.not. interp_elev)) then
         rcode= NF_DEF_VAR (ncid, 'elev_chk', NF_DOUBLE, 2, ndims, varid5)
         rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', 9, 'mean elevation')
         rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 1, 'm')
         rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)
      endif

!   if (interp_elev) then
      ndims(1)= londim ;  ndims(2)= latdim
      rcode= NF_DEF_VAR (ncid, 'elev_sd', NF_DOUBLE, 2, ndims, varid2)
      rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 31, 'standard deviation of elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      rcode= NF_DEF_VAR (ncid, 'elev_min', NF_DOUBLE, 2, ndims, varid3)
      rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 17, 'minimum elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      rcode= NF_DEF_VAR (ncid, 'elev_max', NF_DOUBLE, 2, ndims, varid4)
      rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 17, 'maximum elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)
!   endif

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

!    elevation data
      start= 1 ;  count(1)= id ;  count(2)= jd
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, elev(:,:,n))

      if (ntiles == 1 .and. (.not. interp_elev)) then
         start= 1 ;  count(1)= id ;  count(2)= jd
         rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, elev_chk)
      endif

!   if (interp_elev) then
      start= 1 ;  count(1)= id ;  count(2)= jd
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, elev_sd(:,:,n))

      start= 1 ;  count(1)= id ;  count(2)= jd
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, elev_min(:,:,n))

      start= 1 ;  count(1)= id ;  count(2)= jd
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, elev_max(:,:,n))
!   endif

!  close netcdf file
      rcode= NF_CLOSE (ncid)
   enddo


   if (interp_elev) then
! write 1-deg elevation data as a check
      write (fname, '(a)') 'elev_data_1deg.nc'
      rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
      rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (grid_x, grid_y)
      rcode= NF_DEF_DIM (ncid, 'lon',  ido,   londim)
      rcode= NF_DEF_DIM (ncid, 'lat',  jdo,   latdim)

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

      ndims(1)= londim ;  ndims(2)= latdim
      rcode= NF_DEF_VAR (ncid, 'elev', NF_DOUBLE, 2, ndims, varid)
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      ndims(1)= londim ;  ndims(2)= latdim
      rcode= NF_DEF_VAR (ncid, 'elev_sd', NF_DOUBLE, 2, ndims, varid2)
      rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 31, 'standard deviation of elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      ndims(1)= londim ;  ndims(2)= latdim
      rcode= NF_DEF_VAR (ncid, 'elev_min', NF_DOUBLE, 2, ndims, varid3)
      rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 17, 'minimum elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)

      ndims(1)= londim ;  ndims(2)= latdim
      rcode= NF_DEF_VAR (ncid, 'elev_max', NF_DOUBLE, 2, ndims, varid4)
      rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 17, 'maximum elevation')
      rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 1, 'm')
      rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)

!  leave define mode
      rcode= NF_ENDDEF (ncid)

!  write coordinate data
      start= 1 ;  count= 1

      count(1)= ido
      rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lono)

      count(1)= jdo
      rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lato)

!    elevation data
      start= 1 ;  count(1)= ido ;  count(2)= jdo
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, elv_cnv)

      start= 1 ;  count(1)= ido ;  count(2)= jdo
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, elv_sd)

      start= 1 ;  count(1)= ido ;  count(2)= jdo
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, elv_min)

      start= 1 ;  count(1)= ido ;  count(2)= jdo
      rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, elv_max)

!  close netcdf file
      rcode= NF_CLOSE (ncid)

      deallocate (lato, lono, latob, lonob)
      deallocate (arlato, elv_cnv, elv_sd, elv_min, elv_max)
   endif

   close (10)

   deallocate (lat_idx, lon_idx)
   deallocate (lat, lon, land_frac)

   stop

contains

   subroutine read_namelist(file_path, ntiles, ide, jde, rese,ido, jdo, ncell_cnv, reso)
      character(len=*),  intent(in)    :: file_path
      integer ,intent(inout) :: ntiles
      integer ,intent(inout) :: ide, jde
      real ,intent(inout)   :: rese
      integer ,intent(inout) :: ido, jdo
      integer ,intent(inout) :: ncell_cnv
      real , intent(inout)   :: reso

      integer :: fu, rc
      logical :: there


      !! Namelist definition.
      namelist /cp_elev_nml/ ntiles, ide, jde, rese,ido, jdo, ncell_cnv, reso

      ! Check whether file exists.
      inquire (file=file_path, EXIST=there )

      if (.not. there) then
         write (6, '("Error: input file ", a, " does not exist")') file_path
         stop
         !return
      else
        write (6, '("About to read namelist file ", a)') file_path
      endif

      ! Open and read Namelist file.
      open (action='read', file=file_path, iostat=rc, newunit=fu)
      read (nml=cp_elev_nml, iostat=rc, unit=fu)
      if (rc /= 0) then
        write (6, '("Error: invalid namelist format for file ", a)') file_path
        close (fu)
        stop
      else
        write (6, '("The file namelist ", a," was read. The namelist is:")')
        write(6,cp_elev_nml)
      endif
      close (fu)
   end subroutine read_namelist
end



