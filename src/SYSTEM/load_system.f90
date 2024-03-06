
subroutine load_system ()
  use M_system
  use M_constants!solo para ini const
  use M_fdata, only: symbolA, nspecies
  implicit none
  integer iatom, ikpoint
  integer in1
  integer ispec
  logical zindata

  write(*,*) symbolA, nspecies
  open (unit = 69, file = 'input.xyz', status = 'old')
  read (69, *) natoms
  read (69,*)
  allocate (ratom (3, natoms))
  allocate (symbol (natoms))
  allocate (imass (natoms))
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

  !Latice vectors
  a1vec(1) = 100
  a1vec(2) = 0
  a1vec(3) = 0
  a2vec(1) = 0
  a2vec(2) = 100
  a2vec(3) = 0
  a3vec(1) = 0
  a3vec(2) = 0
  a3vec(3) = 100
  
  !kpts cargamos Gamma solo
  nkpoints = 1
  allocate (special_k(3, nkpoints))
  allocate (special_k_orig(3, nkpoints))
  allocate (scale_k(3, nkpoints))
  allocate (weight_k(nkpoints))
  allocate (weight_k_orig(nkpoints))
  special_k_orig(:,1) = 0
  weight_k_orig(1) = 1
  weight_k(1) = 1

  do ikpoint = 1, nkpoints
    special_k(:,ikpoint) = special_k_orig(:,ikpoint)
    weight_k(ikpoint) = weight_k_orig(ikpoint)
  end do
  
  call allocate_system()
  
  call scf_loop ()
  
  call getenergy ()

  ! call postscf () cuando queramos hacer DOS
  ! call getforces()
end subroutine


