
subroutine load_system ()
  use M_system
  integer iatom
  integer in1
  integer ispec
  logical zindata

  open (unit = 69, file = 'input.xyz', status = 'old')
  read (69, *) natoms
  close (unit = 69)

  call allocate_system()
  open (unit = 69, file = 'input.xyz', status = 'old')
  read (69, *) natoms
  read (69,*)
  do iatom = 1, natoms
   read (69,*) symbol(iatom),ratom(:,iatom)
   zindata = .false.
   do ispec = 1, nspecies
     if (trim(symbol(iatom)) .eq. symbolA(ispec)) then
       zindata = .true.
       imass(iatom) = ispec
     end if
   end do
  end do
  close (unit = 69)

  call scf_loop ()

  call getenergy ()

  ! call postscf () cuando queramos hacer DOS
  ! call getforces()
end subroutine


