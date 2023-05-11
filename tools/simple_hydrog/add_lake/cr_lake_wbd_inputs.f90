program cr_lake_wbd_inputs

! code is the same as cr_lake_files_new except that gage data is
!   read for each tile, and included in the output netcdf file

use horiz_interp_mod

implicit none

include "param.h"  !  get ntiles here

integer, parameter :: maxdims= 3
integer, parameter :: ni_cells= 5, nj_cells= 5
integer, parameter :: nlkmx_all= 1000000
integer, parameter :: nlkmx= 100000
integer, parameter :: nlake_def= 16, nlake_defp1= nlake_def+1
integer, parameter :: iaral= 15, ibalk= 13, ichad= 16, icaspian= nlake_def+1
integer, parameter :: imich= 1, ihuron= 2, ierie= 10, iont= 12
integer, parameter :: idl= 2880, jdl= 1440
integer, parameter :: idlp1= idl+1, jdlp1= jdl+1
integer, parameter :: id_720= 720, jd_720= 720, ntile_720= 6
integer, parameter :: ido= 360, jdo= 180
integer, parameter :: idop1= ido+1, jdop1= jdo+1
integer, parameter :: ncell_cnv= 8
integer, parameter :: nreg= 8
integer, parameter :: k_ice= 9

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: lfrac_const= 0.20
real, parameter :: mval_mdl= -9999.
real, parameter :: casp_area= 3.71e11  ! m2
real, parameter :: lake_depth_small= 2.
real, parameter :: reso= 1.
real, parameter :: ldepth_chad= 2.   ! depth of lake chad in m
real, parameter :: ldepth_max= 50.
real, parameter :: res_fine= 0.0625

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, k, rcode2, ndm,  ii, jj, idp2, jdp2, ktr2
integer :: londim, latdim, fill_val, j1, kt_wbd, latgid, longid
integer :: varid2, varid3, varid4, varid5, idp3, idp4, jdp3, jdp4
integer :: ip, jp, j2, i1, i2, id_ice, jd_ice, idp1_ice, jdp1_ice
integer :: zid, ntype, id_fine, jd_fine, idp1_fine, jdp1_fine
real :: pi, dtr, sum, mval_ice, mval_tocell, mval_landf, mval_cella
real :: mval_ctn, scale, mval_travel, mval_ctn2, travel_thresh
real :: sum1, asum1
character(len=8)   :: var
character(len=40)  :: var_units, zax_name, lname_ice, lake_mask_var
character(len=200) :: fname, ctn_latlon_field, glcc_file, lname_glcc, cover_type_file
logical :: create_large_lakes, screen_small_wbd, interp_glcc_to_1deg, lgm_grid

character(len=8)   :: vname_glcc= 'WaterBod'
   
integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ntiles)            :: itw, ite, its, itn
integer, dimension (nlake_defp1)       :: ktr
integer, dimension (nlkmx_all)         :: il1, jl1, nl1, idx1
integer, dimension (nlkmx_all)         :: iwscr, jwscr, nwscr
integer, dimension (nlkmx,nlake_defp1) :: ilake, jlake, nlake
integer, dimension (idl,jdl)           :: idat

real, dimension (idl)                  :: lonl, lonin
real, dimension (jdl)                  :: latl, latin
real, dimension (ido)                  :: lono
real, dimension (jdo)                  :: lato
real, dimension (idlp1)                :: lonlb
real, dimension (jdlp1)                :: latlb
real, dimension (idop1)                :: lonob
real, dimension (jdop1)                :: latob
real, dimension (ido,jdo)              :: arlato
real, dimension (nlkmx_all)            :: fr1, wbscr
real, dimension (idl,jdl)              :: arlatl, wbdat
real, dimension (ido,jdo)              :: wbd_cnv
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
   
! for c720, lots of land points that are really ocean and so show
!  up as lake --> get rid of these lake points (Sea of Azov, 
!  northern Canada, all on tile 3)

integer, dimension (nreg,ntile_720)          :: ib1 = &
(/ 695, 599,  -1,  -1,  -1,  -1,  -1,  -1, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
    93, 277, 351, 381, 393, 411, 314, 496, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   278,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
    83,  42,  -1,  -1,  -1,  -1,  -1,  -1  /)

integer, dimension (nreg,ntile_720)          :: ib2 = &
(/ 720, 607,  -1,  -1,  -1,  -1,  -1,  -1, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   127, 350, 380, 392, 410, 418, 359, 510, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   287,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   130,  45,  -1,  -1,  -1,  -1,  -1,  -1  /)

integer, dimension (nreg,ntile_720)          :: jb1 = &
(/ 646, 697,  -1,  -1,  -1,  -1,  -1,  -1, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
    83, 497, 497, 497, 497, 497, 185, 637, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   570,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   504, 527,  -1,  -1,  -1,  -1,  -1,  -1  /)

integer, dimension (nreg,ntile_720)          :: jb2 = &
(/ 664, 704,  -1,  -1,  -1,  -1,  -1,  -1, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   103, 540, 527, 540, 526, 509, 219, 681, &
    -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   580,  -1,  -1,  -1,  -1,  -1,  -1,  -1, &
   516, 531,  -1,  -1,  -1,  -1,  -1,  -1  /)


character(len=100), dimension (ntiles) :: river_input_file

integer, allocatable, dimension (:,:) :: big_lake
real, allocatable, dimension (:)      :: lat_idx, lon_idx, zdat
real, allocatable, dimension (:)      :: lonf, latf, lonfb, latfb
real, allocatable, dimension (:)      :: lat_ice, latb_ice, lon_ice, lonb_ice
real, allocatable, dimension (:,:)    :: interp_out, interp_mask, data_in, ice_cover
real, allocatable, dimension (:,:)    :: arlatf, ctn
real, allocatable, dimension (:,:,:)  :: lat, lon, tocell, land_frac, wbd, &
                                         cell_area, travel, ctn2, lake_frac, &
                                         lake_idx

pi= 4.*atan(1.)
dtr= pi/180.
interp_glcc_to_1deg= .false.
lgm_grid= .false.
cover_type_file= ""

open (10, file= 'out.cr_lake_wbd_inputs', form= 'formatted')

do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') create_large_lakes
read (5,'(l1)') screen_small_wbd
read (5,'(a)') ctn_latlon_field
read (5,'(a)') glcc_file
read (5,*) travel_thresh
read (5,'(l1)') interp_glcc_to_1deg
read (5,'(l1)') lgm_grid
if (lgm_grid) then
    read (5,'(a)') cover_type_file
endif
read (5,'(a)') lake_mask_var
close (5)

do n= 1,ntile_720
   write (6,*) 'tile ', n
   do k= 1,nreg
      write (6,'(5i6)') k, ib1(k,n), ib2(k,n), jb1(k,n), jb2(k,n)
   enddo
enddo

write (6,*) 'fort.5 data read'
write (6,*) 'interp_glcc_to_1deg= ', interp_glcc_to_1deg

if (create_large_lakes .or. (screen_small_wbd .and. lgm_grid)) then

!  read connected_to_next field that has been fregridded from
!   c360 to 1/8-deg resolution

    rcode= NF_OPEN (trim(ctn_latlon_field), NF_NOWRITE, ncid)
    if (rcode /= 0) then
        write (6,*) "ERROR: cannot open lat/lon ctn netcdf file"  
        write (6,*) trim(ctn_latlon_field)
        stop 100
    endif

    rcode= nf_inq_varid (ncid, 'grid_y', latid)         ! number of lats
    if (rcode /= 0) then
        write (6,*) "ERROR: cannot find connected_to_next lat variable" ; stop 102
    endif
    rcode= nf_inq_vardimid (ncid, latid, dimids)
    rcode= nf_inq_dimlen (ncid, dimids(1), jd_fine)
    jdp1_fine= jd_fine+1
    allocate (latf(jd_fine), latfb(jdp1_fine))
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
    rcode= nf_inq_dimlen (ncid, dimids(1), id_fine)
    idp1_fine= id_fine+1
    allocate (lonf(id_fine), lonfb(idp1_fine))
    start= 1 ;  count= 1 ;  count(1)= id_fine
    rcode= nf_get_vara_double (ncid, lonid, start, count, lonf)
    
    lonfb(1)= lon1
    lonfb(idp1_fine)= lon2
    do i= 2,id_fine
       lonfb(i)= 0.5*(lonf(i)+lonf(i-1))
    enddo

! compute area of latitude
    allocate (arlatf(id_fine,jd_fine))
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
    
    allocate (ctn(id_fine,jd_fine))

    ctn= mval_mdl
    write (6,*) 'read field: ', lake_mask_var
    rcode= nf_inq_varid (ncid, trim(lake_mask_var), varid)
    if (rcode /= 0) then
        write (6,*) "ERROR: cannot find ", trim(lake_mask_var)
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
    
    deallocate (arlatf)

endif


if (create_large_lakes .or. screen_small_wbd) then
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
    write (6,*) 'fill_val= ', fill_val
   
    scale= -99999999.
    rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
    if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
    write (6,*) 'scale= ', scale

    rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
    if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', lname_glcc)
    write (6,*) 'long_name= ', trim(lname_glcc)
   
    var_units= ' '
    rcode= nf_get_att_text (ncid, varid, "units", var_units)
    write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
    do j= 1,jdl
       j1= jdlp1-j
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

endif


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
          asum1= 0. ;  sum1= 0.
          do jj= j1,j2
             do ii= i1,i2
                if (wbdat(ii,jj) /= mval_mdl) then
                    sum= sum + arlatl(ii,jj)
                    asum1= asum1 + arlatl(ii,jj)
                    sum1= sum1 + wbdat(ii,jj)*arlatl(ii,jj)
                endif
             enddo
          enddo
          if (abs(asum1 - arlato(i,j)) < 1.0e-10) then
              wbd_cnv(i,j)= sum1/asum1
          else
              write (6,*) "ERROR: wbd 1-deg areas do not agree, ", asum1, arlato(i,j)
              stop 17
          endif
       enddo
    enddo
    write (6,*) 'sum of conversion area= ', sum
endif


! ----------------------------------------------------------------------
! if lgm grid, then read ice data
! ----------------------------------------------------------------------

if (screen_small_wbd .and. lgm_grid) then
   rcode= NF_OPEN (trim(cover_type_file), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
        write (6,*) "ERROR: cannot open netcdf file"  ; stop 12
   endif

   start= 1 ; count= 1

! ----------------------------------------------------------------------
! get latitudes and longitudes of input cover file
! ----------------------------------------------------------------------

   rcode= nf_inq_varid (ncid, 'LAT', latid)         ! number of lats
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'lat', latid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find lat variable" ; stop 22
       endif
   endif
   rcode= nf_inq_vardimid (ncid, latid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), jd_ice)
   write (6,*) 'jd_ice= ', jd_ice

   allocate (lat_ice(jd_ice))

   count(1)= jd_ice
   rcode= nf_get_vara_double (ncid, latid, start, count, lat_ice)

   if (lat_ice(1) > 0.) then
       write (6,*) "ERROR: lats in wrong direction"
       stop 10
   endif

   rcode= nf_inq_varid (ncid, 'LATB', latid)         ! number of lats
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'latb', latid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find latb variable" ; stop 22
       endif
   endif
   rcode= nf_inq_vardimid (ncid, latid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), jdp1_ice)
    
   allocate (latb_ice(jdp1_ice))
   count(1)= jdp1_ice
   rcode= nf_get_vara_double (ncid, latid, start, count, latb_ice)
   write (6,*) 'jdp1_ice= ', jdp1_ice

   rcode= nf_inq_varid (ncid, 'LON', lonid)         ! number of lons
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'lon', lonid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find lon variable" ; stop 23
       endif
   endif
   rcode= nf_inq_vardimid (ncid, lonid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), id_ice)
   write (6,*) 'id_ice= ', id_ice

   allocate (lon_ice(id_ice))

   count(1)= id_ice
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon_ice)

   if (lon_ice(1) < 0.) then
       write (6,*) "ERROR: lons in wrong direction"
       stop 11
   endif

   rcode= nf_inq_varid (ncid, 'LONB', lonid)         ! number of lats
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'lonb', lonid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find lonb variable" ; stop 23
       endif
   endif
   rcode= nf_inq_vardimid (ncid, lonid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), idp1_ice)
    
   allocate (lonb_ice(idp1_ice))
   count(1)= idp1_ice
   rcode= nf_get_vara_double (ncid, lonid, start, count, lonb_ice)
   write (6,*) 'idp1_ice= ', idp1_ice


   write (10,'(/"input file lats")')
   do j= 1,jd_ice
      write (10,'(i5,6f12.4)') j, lat_ice(j), latb_ice(j), latb_ice(j + 1)
   enddo

   write (10,'(/"input file lons")')
   do i= 1,id_ice
      write (10,'(i5,4f12.4)') i, lon_ice(i), lonb_ice(i), lonb_ice(i + 1)
   enddo
   write (10,'(/)')


! ----------------------------------------------------------------------
! get z values (cover types)
! ----------------------------------------------------------------------

   write (zax_name,'(a,i2)') 'zax1_10'
   rcode= nf_inq_varid (ncid, trim(zax_name), zid)         ! number of lats
   if (rcode /= 0) then
       write (zax_name,'(a,i2)') 'ZAX1_10'
       rcode2 = nf_inq_varid (ncid, trim(zax_name), zid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find z variable" ; stop 26
       endif
   endif
   rcode= nf_inq_vardimid (ncid, zid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ntype)

   allocate (zdat(ntype))

   count(1)= ntype
   rcode= nf_get_vara_double (ncid, zid, start, count, zdat)

   write (6,*) 'zdat= ', zdat


! ----------------------------------------------------------------------
!  now read cover field
! ----------------------------------------------------------------------

   write (6,*) 'read frac'
   rcode= nf_inq_varid (ncid, 'frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 30
   endif

   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), n)
   if (n /= id_ice) then
       write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', id_ice= ', id_ice
       stop 31
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), n)
   if (n /= jd_ice) then
       write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', jd_ice= ', jd_ice
       stop 32
   endif
   rcode= nf_inq_dimlen (ncid, dimids(3), n)
   if (n /= ntype) then
       write (6,*) 'ERROR: inconsistent dimensions, n= ', n, ', ntype= ', ntype
       stop 33
   endif


   var_units= ' '
   if (rcode == 0) then
       rcode= nf_get_att_text (ncid, varid, "units", var_units)
   endif
   write (6,*) 'units= ', trim(var_units)

   lname_ice= ' '
   rcode= nf_get_att_text (ncid, varid, "long_name", lname_ice)
   write (6,*) 'long_name= ', trim(lname_ice)

   mval_ice= -99999.
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_ice)
   endif
   write (6,*) 'mval= ', mval_ice

   allocate (ice_cover(id_ice,jd_ice))

   ice_cover= mval_ice
   start= 1 ; count= 1
   start(3)= k_ice; count(1)= id_ice ; count(2)= jd_ice ; count(3)= 1
   rcode= nf_get_vara_double (ncid, varid, start, count, ice_cover)

   rcode= nf_close (ncid)

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
jdp3= jd + 3
jdp4= jd + 4

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
idp3= id + 3
idp4= id + 4

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)
  
rcode= nf_close (ncid)



! ----------------------------------------------------------------------
! now open river files -- read lat,lon grids, tocell, land_frac, 
!   cellarea
! ----------------------------------------------------------------------

allocate (lat(idp4,jdp4,ntiles), lon(idp4,jdp4,ntiles), travel(idp4,jdp4,ntiles))
allocate (tocell(idp4,jdp4,ntiles), land_frac(idp4,jdp4,ntiles), cell_area(idp4,jdp4,ntiles))
allocate (ctn2(idp4,jdp4,ntiles), lake_frac(idp4,jdp4,ntiles))

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
   rcode= nf_get_vara_double (ncid, latid, start, count, lat(3:idp2,3:jdp2,n))
       
       
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
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon(3:idp2,3:jdp2,n))


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
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(3:idp2,3:jdp2,n))

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
   rcode= nf_get_vara_double (ncid, varid, start, count, land_frac(3:idp2,3:jdp2,n))

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
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_area(3:idp2,3:jdp2,n))

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
   

   write (6,*) 'read travel'
   rcode= nf_inq_varid (ncid, 'travel', varid)     ! travel field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find travel" ; stop 70
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent travel dimension, ", ndm, id ;  stop 75
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent travel dimension, ", ndm, jd ;  stop 75
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, travel(3:idp2,3:jdp2,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
   write (6,*) 'units= ', var_units

   mval_travel= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_travel)
   endif
   write (6,*) 'mval= ', mval_travel

   where (travel(:,:,n) == mval_travel) travel(:,:,n)= mval_mdl
   
   if (.not. create_large_lakes) then  
       write (6,*) 'read connected_to_next'
       rcode= nf_inq_varid (ncid, 'connected_to_next', varid)     ! connected_to_next field
       if (rcode /= 0) then
           write (6,*) "ERROR: cannot find connected_to_next" ; stop 80
       endif
       rcode= nf_inq_vardimid (ncid, varid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= id) then
           write (6,*) "ERROR: inconsistent connected_to_next dimension, ", ndm, id ;  stop 85
       endif
       rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
       if (ndm /= jd) then
           write (6,*) "ERROR: inconsistent connected_to_next dimension, ", ndm, jd ;  stop 85
       endif

       start= 1 ;  count(1)= id ;  count(2)= jd
       rcode= nf_get_vara_double (ncid, varid, start, count, ctn2(3:idp2,3:jdp2,n))

       var_units= ' '
       rcode= nf_get_att_text (ncid, varid, "units", var_units)
       write (6,*) 'units= ', var_units

       mval_ctn2= 1.e+20
       rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
       if (rcode == 0) then
           rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_ctn2)
       endif
       write (6,*) 'mval= ', mval_ctn2

       where (ctn2(:,:,n) == mval_ctn2) ctn2(:,:,n)= mval_mdl
       lake_frac(:,:,n)= ctn2(:,:,n)
   endif
   
   rcode= nf_close (ncid)
   
!   write (10,'(/"river lats, tile", i4)') n
!   do j= 3,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lat(i,j,n), i= 3,idp2)
!   enddo

!   write (10,'(/"river lons, tile", i4)') n
!   do j= 3,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.4)') (lon(i,j,n), i= 3,idp2)
!   enddo

!   write (10,'(/"cell_area, tile", i4)') n
!   do j= 3,jdp2
!      write (10,*) 'j= ', j
!      write (10,'(10f10.1)') (cell_area(i,j,n)/1.e6, i= 3,idp2)
!   enddo

enddo

sum= 0.
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
         sum= sum + cell_area(i,j,n)
      enddo
   enddo
enddo
write (6,'(/"sum of cell_area= ",f15.1)') sum/1.e+6


if (create_large_lakes .or. screen_small_wbd) then
! ----------------------------------------------------------------------
! interpolate glcc water bodies and permanent wetlands to model grid
! ----------------------------------------------------------------------

  allocate (wbd(idp4,jdp4,ntiles), interp_out(id,jd))

  if (interp_glcc_to_1deg) then
      allocate (interp_mask(ido,jdo))
      allocate (data_in(ido,jdo))
      interp_mask= 1.
      where (wbd_cnv(:,:) == mval_mdl) interp_mask= 0.

      do n= 1,ntiles
         write (6,*) 'WaterBod, tile= ', n
         data_in= wbd_cnv(:,:)
         call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(3:idp2,3:jdp2,n)*dtr, &
            lat(3:idp2,3:jdp2,n)*dtr, wbd(3:idp2,3:jdp2,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         where (interp_out(:,:) == 0.) wbd(3:idp2,3:jdp2,n)= mval_mdl
      enddo
  else
      allocate (interp_mask(idl,jdl))
      allocate (data_in(idl,jdl))

      interp_mask= 1.
      where (wbdat(:,:) == mval_mdl) interp_mask= 0.

      do n= 1,ntiles
         write (6,*) 'WaterBod, tile= ', n
         data_in= wbdat(:,:)
         call horiz_interp (data_in, lonlb*dtr, latlb*dtr, lon(3:idp2,3:jdp2,n)*dtr, &
            lat(3:idp2,3:jdp2,n)*dtr, wbd(3:idp2,3:jdp2,n), verbose=1, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
         where (interp_out(:,:) == 0.) wbd(3:idp2,3:jdp2,n)= mval_mdl
      enddo
   endif
    
    where (tocell == mval_mdl)
       wbd= mval_mdl
    endwhere
! set ocean values to missing
    where (land_frac == 0.)
       wbd= mval_mdl
    endwhere
! set part-land values to zero
    where (land_frac > 0 .and. land_frac < 1.)
       wbd= 0.
    endwhere
    
    deallocate (interp_out, interp_mask)
    deallocate (data_in)
endif


if (create_large_lakes) then

   allocate (big_lake(id_fine,jd_fine))

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

    write (6,*) 43, 89, lat(89,43,1), lon(89,43,1), wbd(89,43,1)
    write (6,*) 43, 90, lat(90,43,1), lon(90,43,1), wbd(90,43,1)
! now if ctn lat/lon field has value > 0, assign as lake in c720 field/tile
    ktr2= 0
    do j= 1,jd_fine
       if (mod(j,100) == 0) write (6,*) 'j= ', j
       do i= 1,id_fine
          if (ctn(i,j) > 0.3) then
              do n= 1,ntiles
                 do jj= 3,jdp2
                    do ii= 3,idp2
                       if (lat(ii,jj,n) > latfb(j) .and. lat(ii,jj,n) <= latfb(j+1) .and. &
                           lon(ii,jj,n) > lonfb(i) .and. lon(ii,jj,n) <= lonfb(i+1) .and. &
                           wbd(ii,jj,n) > 0. .and. land_frac(ii,jj,n) == 1.) then
                               ktr2= ktr2 + 1
                               il1(ktr2)= ii
                               jl1(ktr2)= jj
                               nl1(ktr2)= n
                               idx1(ktr2)= big_lake(i,j)
                               fr1(ktr2)= 1.
                       endif
                    enddo
                 enddo
              enddo
          endif
       enddo
    enddo
    write (6,*) 'ktr2= ', ktr2

    allocate (lake_idx(idp4,jdp4,ntiles))
    ktr= 0 ;  ilake= 0 ;  jlake= 0 ;  nlake= 0
    ctn2= 0. ;  lake_frac= 0. ;  lake_idx= 0.
    do l= 1,ktr2
       do n= 1,nlake_defp1
          if (idx1(l) == n) then
              ktr(n)= ktr(n) + 1
              ilake(ktr(n),n)= il1(l)-2
              jlake(ktr(n),n)= jl1(l)-2
              nlake(ktr(n),n)= nl1(l)
              frlake(ktr(n),n)= fr1(l)
              ctn2(il1(l),jl1(l),nl1(l))= 1.
              lake_frac(il1(l),jl1(l),nl1(l))= fr1(l)
              lake_idx(il1(l),jl1(l),nl1(l))= real(n)
          endif
       enddo
    enddo
    where (tocell == mval_mdl)
       ctn2= mval_mdl
       lake_frac= mval_mdl
       lake_idx= mval_mdl
    endwhere

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
    
endif

if (screen_small_wbd) then
! okay, now set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

    if (ntiles == 1) then
        itw(n)= 1 ;  ite(n)= 1 ;  its(n)= 1 ;  itn(n)= 1
        call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, ctn2)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, wbd)
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
        enddo
       
        call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, ctn2)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, wbd)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
        call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
    endif

    kt_wbd= 0
    if (lgm_grid) then
! if ctn lat/lon field is not defined, allow no waterbod fraction
        do j= 1,jd_fine
           if (mod(j,100) == 0) write (6,*) 'j= ', j
           do i= 1,id_fine
              if (ctn(i,j) == mval_ctn) then
                  do n= 1,ntiles
                     do jj= 3,jdp2
                        do ii= 3,idp2
                           if (lat(ii,jj,n) > latfb(j) .and. lat(ii,jj,n) <= latfb(j+1) .and. &
                               lon(ii,jj,n) > lonfb(i) .and. lon(ii,jj,n) <= lonfb(i+1) .and. &
                               wbd(ii,jj,n) > 0. .and. land_frac(ii,jj,n) > 0.) then
                                   kt_wbd= kt_wbd + 1
                                   nwscr(kt_wbd)= n
                                   jwscr(kt_wbd)= jj-2
                                   iwscr(kt_wbd)= ii-2
                                   wbscr(kt_wbd)= wbd(ii,jj,n)
                                   wbd(ii,jj,n)= 0.
                           endif
                        enddo
                     enddo
                  enddo
              endif
           enddo
        enddo
        write (6,*) 'screen ocean points, kt_wbd= ', kt_wbd
        
! also allow no waterbod fraction over ice field
        ktr2= 0
        do j= 1,jd_ice
           do i= 1,id_ice
              if (ice_cover(i,j) /= mval_ice .and. ice_cover(i,j) > 0.) then 
                  do n= 1,ntiles
                     do jj= 3,jdp2
                        do ii= 3,idp2
                           if (lat(ii,jj,n) > latb_ice(j) .and. lat(ii,jj,n) <= latb_ice(j+1) .and. &
                               lon(ii,jj,n) > lonb_ice(i) .and. lon(ii,jj,n) <= lonb_ice(i+1) .and. &
                               wbd(ii,jj,n) > 0. ) then
                                   kt_wbd= kt_wbd + 1
                                   ktr2= ktr2 + 1
                                   nwscr(kt_wbd)= n
                                   jwscr(kt_wbd)= jj-2
                                   iwscr(kt_wbd)= ii-2
                                   wbscr(kt_wbd)= wbd(ii,jj,n)
                                   wbd(ii,jj,n)= 0.
                           endif
                        enddo
                     enddo
                  enddo
              endif
           enddo
        enddo
        write (6,*) 'screen ice points, ktr2= ', ktr2, ', kt_wbd= ', kt_wbd
    endif

    do n= 1,ntiles
       do j= 3,jdp2
          do i= 3,idp2
             if (tocell(i,j,n) == mval_mdl .or. wbd(i,j,n) == mval_mdl) go to 109
             if (lat(i,j,n) < -60.) then                    ! antarctica
                 if (wbd(i,j,n) == 0.) go to 109
             else if (lat(i,j,n) > 70. ) then               ! arctic
                 if (wbd(i,j,n) == 0.) go to 109
             else if (lat(i,j,n) > 58. .and. &              ! rest of greenland
                 (lon(i,j,n) > 302. .and. lon(i,j,n) < 348.)) then
                 if (wbd(i,j,n) == 0.) go to 109
             else if (id == id_720 .and. jd == jd_720 .and. ntiles == ntile_720) then
                 do k= 1,nreg
                    if (i-2 >= ib1(k,n) .and. i-2 <= ib2(k,n) .and.   &
                        j-2 >= jb1(k,n) .and. j-2 <= jb2(k,n)) then
                        if (wbd(i,j,n) == 0.) then
                            go to 109
                        else
                            go to 108
                        endif
                    endif
                 enddo
                 if (wbd(i,j,n) < 0.1 .or. travel(i,j,n) > 10. .or. &
                     ctn2(i,j,n) > 0.3) then
                     go to 109
                 endif
             else
              ! insert this next section for global hydrosheds, 0.125-degree, 14 nov 2011
              ! ----------
             !    do jj= 1,nj_cells
             !       jp= j+jj-3
             !       do ii= 1,ni_cells
             !          ip=i+ii-3
             !          if (wbd(ip,jp,n) == mval_mdl) go to 108
             !       enddo
             !    enddo
              ! ----------
                 if (wbd(i,j,n) < 0.1 .or. travel(i,j,n) > travel_thresh .or. &
                     ctn2(i,j,n) > 0.3) then
                     go to 109
                 endif
             endif
108          continue
             kt_wbd= kt_wbd + 1
             nwscr(kt_wbd)= n
             jwscr(kt_wbd)= j-2
             iwscr(kt_wbd)= i-2
             wbscr(kt_wbd)= wbd(i,j,n)
             wbd(i,j,n)= 0.
109          continue
          enddo
       enddo
    enddo
    write (6,*) 'kt_wbd = ', kt_wbd
    
    open (11, file= 'wbd_screen_list', form= 'formatted')
    do n= 1,kt_wbd
       write (11,'(3i5,f9.4)') nwscr(n), jwscr(n), iwscr(n), wbscr(n)
    enddo
    close (11)
    
endif

close (10)


if (create_large_lakes) then
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
    
    deallocate (big_lake)
endif

if (create_large_lakes .or. (screen_small_wbd .and. lgm_grid)) then
    deallocate (latf, lonf, latfb, lonfb, ctn)
endif

! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'lake_wbd_inputs.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))
   
! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (grid_x, grid_y, ngage)
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
   rcode= NF_DEF_VAR (ncid, 'tocell', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 13, 'tocell')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   if (create_large_lakes) then
       rcode= NF_DEF_VAR (ncid, 'lake_idx', NF_DOUBLE, 2, ndims, varid2)
       rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 8, 'lake_idx')
       rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 4, 'none')
       rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   endif
   
   rcode= NF_DEF_VAR (ncid, 'travel', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 15, 'travel')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)
 
   rcode= NF_DEF_VAR (ncid, 'lake_frac', NF_DOUBLE, 2, ndims, varid4)
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 9, 'lake_frac')
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   
   if (screen_small_wbd) then
       rcode= NF_DEF_VAR (ncid, 'waterbod', NF_DOUBLE, 2, ndims, varid5)
       rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', 8, 'waterbod')
       rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 4, 'none')
       rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   endif

!  leave define mode
   rcode= NF_ENDDEF (ncid)

!  write coordinate data
   start= 1 ;  count= 1
      
   count(1)= id
   rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

   count(1)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)
   
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lon(3:idp2,3:jdp2,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, lat(3:idp2,3:jdp2,n))

!    tocell data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, tocell(3:idp2,3:jdp2,n))

!    travel data
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, travel(3:idp2,3:jdp2,n))

!    big lakes
   start= 1 ;  count(1)= id ;  count(2)= jd 
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, lake_frac(3:idp2,3:jdp2,n))

!    big lakes -- index
   if (create_large_lakes) then
       start= 1 ;  count(1)= id ;  count(2)= jd 
       rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, lake_idx(3:idp2,3:jdp2,n))
   endif

   if (screen_small_wbd) then
!    waterbod data
       start= 1 ;  count(1)= id ;  count(2)= jd 
       rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, wbd(3:idp2,3:jdp2,n))
   endif


!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

if (screen_small_wbd .or. create_large_lakes) deallocate (wbd)
if (create_large_lakes) deallocate (lake_idx)
deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area, travel, lake_frac, ctn2)

   
contains

! ----------------------------------------------------------------------
subroutine create_halo (ntl, id, jd, itw, ite, its, itn, field)
! ----------------------------------------------------------------------

implicit none

integer :: n, i, j, ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4

integer, intent(in)     :: ntl, id, jd
integer, intent(in)     :: itw(:), ite(:), its(:), itn(:)

real, intent(inout)     :: field(:,:,:)

ip1= id + 1  ;  jp1= jd + 1
ip2= id + 2  ;  jp2= jd + 2
ip3= id + 3  ;  jp3= jd + 3
ip4= id + 4  ;  jp4= jd + 4


if (ntl == 1) then
    field(1,:,ntl)= field(ip1,:,ntl)
    field(2,:,ntl)= field(ip2,:,ntl)
    field(ip3,:,ntl)= field(3,:,ntl)
    field(ip4,:,ntl)= field(4,:,ntl)
else
    do n= 1,ntl
! define tiles to the west, east, south, and north
       if (mod(n,2) == 0) then
           do j= 3,jp2
              field(1,j,n)=     field(ip1,j,itw(n))  ! western edge
              field(2,j,n)=     field(ip2,j,itw(n))  ! western edge
           enddo
       
           do j= 3,jp2
              i= ip2-j+3
              field(ip3,j,n)=  field(i,3,ite(n))     ! eastern edge
              field(ip4,j,n)=  field(i,4,ite(n))     ! eastern edge
           enddo
          
           do i= 3,ip2
              j= jp2-i+3
              field(i,1,n)=     field(ip1,j,its(n))  ! southern edge
              field(i,2,n)=     field(ip2,j,its(n))  ! southern edge
           enddo

           do i= 3,ip2
              field(i,jp3,n)=  field(i,3,itn(n))     ! northern edge
              field(i,jp4,n)=  field(i,4,itn(n))     ! northern edge
           enddo
       else
           do j= 3,jp2
              i= ip2-j+3
              field(1,j,n)=     field(i,jp1,itw(n))  ! western edge
              field(2,j,n)=     field(i,jp2,itw(n))  ! western edge
           enddo
       
           do j= 3,jp2
              field(ip3,j,n)=  field(3,j,ite(n))     ! eastern edge
              field(ip4,j,n)=  field(4,j,ite(n))     ! eastern edge
           enddo
       
           do i= 3,ip2
              field(i,1,n)=     field(i,jp1,its(n))  ! southern edge
              field(i,2,n)=     field(i,jp2,its(n))  ! southern edge
           enddo
       
           do i= 3,ip2
              j= jp2-i+3
              field(i,jp3,n)=  field(3,j,itn(n))     ! northern edge
              field(i,jp4,n)=  field(4,j,itn(n))     ! northern edge
           enddo
       endif
    enddo

endif

end subroutine create_halo


end






