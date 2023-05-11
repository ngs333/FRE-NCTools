program mod_lake_indices

implicit none

integer, parameter :: nlake= 17

integer :: n_idx_lk, n_idx_ca, n_rtg_lk, np, l, ktr
integer :: n_idx_lk_new, n_idx_ca_new, ichk
real :: lfr

character(len=11), dimension(nlake) :: lake_name= &
(/ 'Michigan   ', 'Huron      ', 'Superior   ', 'Victoria   ', &
   'Tanganyika ', 'Baikal     ', 'Great Bear ', 'Malawi     ', &
   'Great Slave', 'Erie       ', 'Winnipeg   ', 'Ontario    ', &
   'Balkhash   ', 'Ladoga     ', 'Aral       ', 'Chad       ', &
   'Caspian    ' /)

integer, allocatable, dimension (:)      :: t_idxl, j_idxl, i_idxl, ilk
integer, allocatable, dimension (:)      :: t_idxc, j_idxc, i_idxc, ica
integer, allocatable, dimension (:)      :: t_idxl_new, j_idxl_new, i_idxl_new, ilk_new
integer, allocatable, dimension (:)      :: t_idxc_new, j_idxc_new, i_idxc_new, ica_new
integer, allocatable, dimension (:)      :: t_rtg, j_rtg, i_rtg

read (5,*) n_idx_lk
read (5,*) n_idx_ca
read (5,*) n_rtg_lk

write (6,*) "n_idx_lk= ", n_idx_lk
write (6,*) "n_idx_ca= ", n_idx_ca
write (6,*) "n_rtg_lk= ", n_rtg_lk

allocate (t_idxl(n_idx_lk), j_idxl(n_idx_lk), i_idxl(n_idx_lk), ilk(n_idx_lk))
allocate (t_idxl_new(n_idx_lk), j_idxl_new(n_idx_lk), i_idxl_new(n_idx_lk), ilk_new(n_idx_lk))

open (20, form= "formatted")
do np= 1,n_idx_lk
   read (20,*) t_idxl(np), j_idxl(np), i_idxl(np), lfr, ilk(np)
enddo
close (20)
write (6,*) "lake indices read"

allocate (t_idxc(n_idx_ca), j_idxc(n_idx_ca), i_idxc(n_idx_ca), ica(n_idx_ca))
allocate (t_idxc_new(n_idx_ca), j_idxc_new(n_idx_ca), i_idxc_new(n_idx_ca), ica_new(n_idx_ca))

open (21, form= "formatted")
do np= 1,n_idx_ca
   read (21,*) t_idxc(np), j_idxc(np), i_idxc(np), lfr, ica(np)
enddo
close (21)
write (6,*) "caspian indices read"

allocate (t_rtg(n_rtg_lk), j_rtg(n_rtg_lk), i_rtg(n_rtg_lk))

open (22, form= "formatted")
do np= 1,n_rtg_lk
   read (22,*) t_rtg(np), j_rtg(np), i_rtg(np)
enddo
close (22)
write (6,*) "lake routings read"

write (6,*) "check lakes:"
ktr= 0
do np= 1,n_idx_lk
   do l= 1,n_rtg_lk
      if (t_idxl(np) == t_rtg(l) .and. j_idxl(np) == j_rtg(l) .and. &
          i_idxl(np) == i_rtg(l)) then
          ktr= ktr+1
          t_idxl_new(ktr)= t_idxl(np)
          j_idxl_new(ktr)= j_idxl(np)
          i_idxl_new(ktr)= i_idxl(np)
          ilk_new(ktr)= ilk(np)
          go to 20
      endif
   enddo
   write (6,*) "remove point: ", t_idxl(np), j_idxl(np), i_idxl(np), ilk(np)
20 continue
enddo
n_idx_lk_new= ktr
write (6,*) "n_idx_lk_new= ", n_idx_lk_new

write (6,*) "check caspian:"
ktr= 0
do np= 1,n_idx_ca
   do l= 1,n_rtg_lk
      if (t_idxc(np) == t_rtg(l) .and. j_idxc(np) == j_rtg(l) .and. &
          i_idxc(np) == i_rtg(l)) then
          ktr= ktr+1
          t_idxc_new(ktr)= t_idxc(np)
          j_idxc_new(ktr)= j_idxc(np)
          i_idxc_new(ktr)= i_idxc(np)
          ica_new(ktr)= ica(np)
          go to 40
      endif
   enddo
   write (6,*) "remove point: ", t_idxc(np), j_idxc(np), i_idxc(np), ica(np)
40 continue
enddo
n_idx_ca_new= ktr
write (6,*) "n_idx_ca_new= ", n_idx_ca_new

open (25, file= "lake_list_16", form= "formatted")
ichk= 0
do np= 1,n_idx_lk_new
   if (ilk_new(np) == ichk) then
       write (25,'(3i5,f6.1,i5)')      t_idxl_new(np), j_idxl_new(np), i_idxl_new(np), 1.0, ilk_new(np)
   else
       write (25,'(3i5,f6.1,i5,3x,a)') t_idxl_new(np), j_idxl_new(np), i_idxl_new(np), 1.0, ilk_new(np), &
                                       lake_name(ilk_new(np))
       ichk= ilk_new(np)
   endif
enddo
close (25)

open (26, file= "lake_list_caspian", form= "formatted")
do np= 1,n_idx_ca_new
   if (np > 1) then
       write (26,'(3i5,f6.1,i5)')      t_idxc_new(np), j_idxc_new(np), i_idxc_new(np), 1.0, ica_new(np)
   else
       write (26,'(3i5,f6.1,i5,3x,a)') t_idxc_new(np), j_idxc_new(np), i_idxc_new(np), 1.0, ica_new(np), &
                                       lake_name(ica_new(np))
   endif
enddo
close (26)


deallocate (t_idxl, j_idxl, i_idxl, ilk)
deallocate (t_idxc, j_idxc, i_idxc, ica)
deallocate (t_idxl_new, j_idxl_new, i_idxl_new, ilk_new)
deallocate (t_idxc_new, j_idxc_new, i_idxc_new, ica_new)
deallocate (t_rtg, j_rtg, i_rtg)

stop
end


