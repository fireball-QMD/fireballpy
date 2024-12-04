! Set int options
subroutine set_options(dipole_method, charges_method, &
    & ismolecule, isgamma, total_charge, mixer_method, &
    & max_iter, mix_order, beta, w0, tol)
  use iso_c_binding
  use M_system, only: igamma, icluster, idipole, iqout, &
    & qstate, ialgmix, max_scf_iterations, idmix, w02, bmix, sigmatol
  implicit none
  integer(c_long), intent(in) :: dipole_method, charges_method, &
    & ismolecule, isgamma, total_charge, mixer_method, max_iter, mix_order
  real(c_double), intent(in) :: beta, w0, tol
  idipole = dipole_method
  iqout = charges_method
  icluster = ismolecule
  igamma = isgamma
  qstate = total_charge
  ialgmix = mixer_method
  max_scf_iterations = max_iter
  idmix = mix_order
  bmix = beta
  w02 = w0*w0
  sigmatol = tol
end subroutine set_options

! Set cell vectors
subroutine set_cell(a1, a2, a3)
  use iso_c_binding
  use M_system, only: a1vec, a2vec, a3vec
  implicit none
  real(c_double), dimension(3), intent(in) :: a1, a2, a3
  a1vec = a1
  a2vec = a2
  a3vec = a3
end subroutine set_cell

! Set coordinates
subroutine set_coords(naux, z, xyz)
  use iso_c_binding
  use M_system, only: natoms, ratom, symbol, imass
  use M_fdata, only : symbolA, nspecies, nzx
  implicit none
  integer(c_long), intent(in) :: naux
  integer(c_long), dimension(naux), intent(in) :: z
  real(c_double), dimension(3, naux), intent(in) :: xyz
  integer(c_long) :: iatom, ispec
  natoms = naux
  if (allocated(ratom)) deallocate(ratom)
  if (allocated(symbol)) deallocate(symbol)
  if (allocated(imass)) deallocate(imass)
  allocate (ratom (3, naux))
  allocate (symbol (naux))
  allocate (imass (naux))
  do iatom = 1, naux
    ratom(:, iatom) = xyz(:, iatom)
    do ispec = 1, nspecies
      if (z(iatom) .eq. nzx(ispec)) then
        imass(iatom) = ispec
        symbol(iatom) = symbolA(ispec)
        exit
      end if
    end do
  end do
end subroutine set_coords

! Provide input charges
subroutine set_initial_charges(natoms, nsh_max, qinput)
  use iso_c_binding
  use M_system, only : Qin
  implicit none
  integer(c_long), intent(in) :: nsh_max, natoms
  real(c_double), dimension(nsh_max, natoms), intent(in) :: qinput
  Qin = qinput
end subroutine set_initial_charges

! Faster set coordinates for simulations
subroutine update_coords(natoms, xyz)
  use iso_c_binding
  use M_system, only: ratom
  implicit none
  integer(c_long), intent(in) :: natoms
  real(c_double), dimension(3, natoms), intent(in) :: xyz
  ratom = xyz
end subroutine update_coords
 
! Set kpoints as needed by Fireball
subroutine set_kpoints(nkpts, kpts, weights)
  use iso_c_binding
  use M_system, only: nkpoints, special_k, weight_k
  implicit none
  integer(c_long), intent(in) :: nkpts
  real(c_double), dimension(3, nkpts), intent(in) :: kpts
  real(c_double), dimension(nkpts), intent(in) :: weights
  if (allocated(special_k)) deallocate(special_k)
  if (allocated(weight_k)) deallocate(weight_k)
  allocate (special_k(3, nkpts))
  allocate (weight_k(nkpts))
  nkpoints = nkpts
  weight_k = weights
  special_k = kpts
end subroutine set_kpoints

! Load FData in the module
subroutine loadfdata_from_path(fdata_path)
  use iso_c_binding
  use M_fdata, only: fdataLocation
  implicit none
  character(len=400), intent(in) :: fdata_path
  fdatalocation=trim(fdata_path)
  call load_fdata()
end subroutine loadfdata_from_path
 
! Allocate all needed arrays. Python gc handles deallocate after work
subroutine call_allocate_system()
  use iso_c_binding
  implicit none
  call allocate_system()
end subroutine call_allocate_system

! Get util system size information
subroutine get_sizes(nsh, norbs_new, norbs)
  use iso_c_binding
  use M_system, only: norbitals_new, norbitals
  use M_fdata, only: nsh_max
  integer(c_long), intent(out) :: nsh, norbs_new, norbs
  nsh = nsh_max
  norbs_new = norbitals_new
  norbs = norbitals
end subroutine get_sizes

! Return initial charges
subroutine get_initial_charges(natoms, nsh_max, qinitial)
  use iso_c_binding
  use M_fdata, only : Qneutral, nssh
  use M_system, only : imass
  implicit none
  integer(c_long), intent(in) :: nsh_max, natoms
  real(c_double), dimension(nsh_max, natoms), intent(inout) :: qinitial
  integer(c_long) :: iatom, issh, in1
  do iatom=1,natoms
    in1 = imass(iatom)
    do issh=1,nssh(in1)
      qinitial(issh, iatom) = Qneutral(issh, in1)
    end do
    do issh=nssh(in1)+1,nsh_max
      qinitial(issh, iatom) = 0.0d0
    end do
  end do
end subroutine get_initial_charges

! Compute the SCF loop
subroutine scf(natoms, nshell, nkpts, nbands, &
    & verbose, fix_charges, shell_charges, eigenvalues, &
    & converged, errno_out, energy, fermi_level, charges)
  use iso_c_binding
  use M_system, only: errno, scf_achieved, etot, efermi, eigen_k, Qin, imass, ifixcharge
  use M_fdata, only: Qneutral, nssh
  implicit none
  integer(c_long), intent(in) :: natoms, nshell, nkpts, nbands
  logical, intent(in) :: verbose, fix_charges
  real(c_double), dimension(nshell, natoms), intent(inout) :: shell_charges
  real(c_double), dimension(nbands, nkpts), intent(inout) :: eigenvalues
  logical, intent(out) :: converged
  integer(c_long), intent(out) :: errno_out
  real(c_double), intent(out) :: energy, fermi_level
  real(c_double), dimension(natoms), intent(out) :: charges
  integer(c_long) :: iatom, issh, in1

  errno = 0
  if(fix_charges) then
    ifixcharge = 1_c_long
  else
    ifixcharge = 0_c_long
  end if
  call scf_loop (verbose)
  converged = scf_achieved
  errno_out = errno
  energy = etot
  fermi_level = efermi
  do iatom = 1, natoms
    in1 = imass(iatom)
    charges(iatom) = 0.0d0
    do issh = 1, nssh(in1)
      charges(iatom) = charges(iatom) + Qneutral(issh,in1) - Qin(issh,iatom)
    end do
  end do
  shell_charges = Qin
  eigenvalues = eigen_k(1:nbands, :)
end subroutine scf

! Get the forces in each atom
subroutine calc_forces(natoms, forces, errno_out)
  use iso_c_binding
  use M_system, only : ftot, errno
  implicit none
  integer(c_long), intent(in) :: natoms
  real(c_double), dimension(3, natoms), intent(inout) :: forces
  integer(c_long), intent(out) :: errno_out
  errno = 0
  call getforces()
  errno_out = errno
  forces = ftot
end subroutine calc_forces

! Get hamiltonian and overlap matrix
subroutine get_hs(norbitals, natoms, orbitals, sdat, hdat)
  use iso_c_binding
  use M_system, only : imass
  use M_fdata, only : num_orb
  implicit none
  integer(c_long), intent(in) :: norbitals, natoms
  integer(c_long), dimension(natoms), intent(inout) :: orbitals
  real(c_double), dimension(norbitals, norbitals), intent(inout) :: sdat, hdat
  integer(c_long) :: iatom
  call geth(sdat, hdat)
  do iatom = 1, natoms
    orbitals(iatom) = num_orb(imass(iatom))
  end do
end subroutine get_hs
