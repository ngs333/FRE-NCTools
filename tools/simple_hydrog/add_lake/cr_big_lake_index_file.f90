program cr_big_lake_index_file

implicit none

integer, parameter :: maxdims= 3
integer, parameter :: nlkmx= 1000000
integer, parameter :: nlake_def= 16, nlake_defp1= nlake_def+1

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, k, m, j1, i1, id, jd, id_fine, jd_fine
integer :: idp1_fine, jdp1_fine, ktr2, rcode, rcode2, varid, ncid
integer :: attnum, latid, lonid, ktr3, ktr4, londim, latdim
real :: pi, dtr, sum, mval_lmask, f1, lake_mask_val
character(len=8)   :: var
character(len=40)  :: var_units, lake_mask_var
character(len=200) :: fname, lmask_latlon_field

character(len=8)   :: vname_glcc= 'WaterBod'
   
integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (nlake_defp1)       :: ktr
integer, dimension (nlkmx,nlake_defp1) :: ilake, jlake

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


   character(len=11), dimension (nlake_defp1) :: lake_name = &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ', &
   'Caspian    ' /)
   
integer, allocatable, dimension (:,:) :: big_lake, ichk
real, allocatable, dimension (:)      :: lonf, latf, lonfb, latfb
real, allocatable, dimension (:,:)    :: arlatf, lmask

pi= 4.*atan(1.)
dtr= pi/180.

open (10, file= 'out.cr_big_lake_index_file', form= 'formatted')

read (5,'(a)') lmask_latlon_field
read (5,'(a)') lake_mask_var
read (5,*) lake_mask_val
close (5)

write (6,*) 'fort.5 data read'

! --------------------------------------------------------------
!  read lake-var field that has been fregridded from
!   c768 to 1/16-deg resolution
! --------------------------------------------------------------

rcode= NF_OPEN (trim(lmask_latlon_field), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open lat/lon lmask netcdf file"  
    write (6,*) trim(lmask_latlon_field)
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
    
allocate (lmask(id_fine,jd_fine))

lmask= mval_mdl
!write (6,*) 'read field: ', lake_mask_var
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
rcode= nf_get_vara_double (ncid, varid, start, count, lmask)

mval_lmask= 1.e+20
rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_lmask)
!write (6,*) 'mval= ', mval_lmask
   
var_units= ' '
rcode= nf_get_att_text (ncid, varid, "units", var_units)
!write (6,*) 'units= ', var_units

rcode= nf_close (ncid)
    
deallocate (arlatf)

allocate (big_lake(id_fine,jd_fine), ichk(id_fine,jd_fine))

! assign big-lake indices
ktr= 0 ;  ilake= 0 ;  jlake= 0 ;  frlake= 0.
big_lake= 0 ; ktr2= 0 ; ktr3= 0 ; ktr4= 0 ; ichk= 0
do j= 1,jd_fine
   do i= 1,id_fine
      if (lmask(i,j) >= lake_mask_val) then
          ktr2= ktr2 + 1
          do n= 1,nlake_defp1
             if (latf(j) > bounds(n) .and. latf(j) <= boundn(n) .and. &
                 lonf(i) > boundw(n) .and. lonf(i) <= bounde(n)) then
                   big_lake(i,j)= n
                   ktr3= ktr3 + 1
                   ktr(n)= ktr(n) + 1
                   ilake(ktr(n),n)= i
                   jlake(ktr(n),n)= j
                   frlake(ktr(n),n)= 1.
                   if (ichk(i,j) == 0) then
                      ichk(i,j)= n
                   else
                      ktr4= ktr4 + 1
                      write (6,*) j, i, ichk(i,j), n, lake_name(ichk(i,j)), lake_name(n)
                      ichk(i,j)= n
                   endif
             endif
          enddo
      endif
   enddo
enddo

write (6,*) "ktr2= ", ktr2, ", ktr3= ", ktr3, ", ktr4= ", ktr4


! sort lake grid cells by j
do l= 1,nlake_defp1
   do n= 1,ktr(l)
      do m= 1,ktr(l)-1
         if (jlake(m,l) > jlake(m+1,l)) then
             j1= jlake(m,l)
             jlake(m,l)= jlake(m+1,l)
             jlake(m+1,l)= j1
             
             i1= ilake(m,l)
             ilake(m,l)= ilake(m+1,l)
             ilake(m+1,l)= i1
             
             f1= frlake(n,l)
             frlake(n,l)= frlake(m+1,l)
             frlake(m+1,l)= f1
         endif
      enddo
   enddo
enddo

! then sort lake grid cells by i
do l= 1,nlake_defp1
   do n= 1,ktr(l)-1
      do m= 1,ktr(l)-1
         if (jlake(m,l) == jlake(m+1,l) .and. ilake(m,l) > ilake(m+1,l)) then
             i1= ilake(m,l)
             ilake(m,l)= ilake(m+1,l)
             ilake(m+1,l)= i1
                 
             f1= frlake(n,l)
             frlake(n,l)= frlake(m+1,l)
             frlake(m+1,l)= f1
         endif
      enddo
   enddo
enddo

open (11, file= 'lake_list_16', form= 'formatted')
do l= 1,nlake_def
   write (6,*) l, ', ktr= ', ktr(l), lake_name(l)
   do n= 1,ktr(l)
      if (n > 1) then
          write (11,'(3i5,f6.1,i5)') 1, jlake(n,l), ilake(n,l), &
                 frlake(n,l), l
      else
          write (11,'(3i5,f6.1,i5,3x,a)') 1, jlake(n,l), ilake(n,l), &
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
          write (11,'(3i5,f6.1,i5)') 1, jlake(n,l), ilake(n,l), &
                 frlake(n,l), l
      else
          write (11,'(3i5,f6.1,i5,3x,a)') 1, jlake(n,l), ilake(n,l), &
                 frlake(n,l), l, lake_name(l)
      endif
   enddo
enddo
close (11)
    
close (10)


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
    
deallocate (big_lake, latf, lonf, latfb, lonfb, lmask)
   
stop

end






