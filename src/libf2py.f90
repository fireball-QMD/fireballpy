! Set int options
subroutine set_options(dipole_method, charges_method, fix_charges, &
    & ismolecule, isgamma, total_charge, mixer_method, &
    & max_iter, mix_order, beta, w0, tol)
  use iso_c_binding
  use M_system, only: igamma, icluster, idipole, iqout, ifixcharge, &
    & qstate, ialgmix, max_scf_iterations, idmix, w02, bmix, sigmatol
  implicit none
  integer(c_long), intent(in) :: dipole_method, charges_method, fix_charges, &
    & ismolecule, isgamma, total_charge, mixer_method, max_iter, mix_order
  real(c_double), intent(in) :: beta, w0, tol
  idipole = dipole_method
  iqout = charges_method
  ifixcharge = fix_charges
  icluster = ismolecule
  qstate = total_charge
  ialgmix = mixer_method
  max_scf_iterations = max_iter
  idmix = mix_order
  w02 = w0*w0
  bmix = beta
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
subroutine set_charges(natoms, nsh_max, qinput)
  use iso_c_binding
  use M_system, only : Qin
  implicit none
  integer(c_long), intent(in) :: nsh_max, natoms
  real(c_double), dimension(nsh_max, natoms), intent(in) :: qinput
  Qin = qinput
end subroutine set_charges

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

! Compute the SCF loop
subroutine scf(verbose, errno_out)
  use iso_c_binding
  use M_system, only: errno
  implicit none
  logical, intent(in) :: verbose
  integer(c_long), intent(out) :: errno_out
  errno = 0
  call scf_loop (verbose)
  errno_out = errno
end subroutine scf

! Execute Dassembles for forces
subroutine calc_forces(errno_out)
  use iso_c_binding
  use M_system, only: errno
  implicit none
  integer(c_long), intent(out) :: errno_out
  errno = 0
  call getforces()
  errno_out = errno
end subroutine calc_forces

! Get util system size information
subroutine get_sizes(nsh, norbs)
  use iso_c_binding
  use M_system, only: norbitals_new
  use M_fdata, only: nsh_max
  integer(c_long), intent(out) :: nsh, norbs
  nsh = nsh_max
  norbs = norbitals_new
end subroutine get_sizes

! Return energy, fermi level
subroutine get_energies(e, ef)
  use iso_c_binding
  use M_system, only: etot, efermi
  implicit none
  real(c_double), intent(out) :: e, ef
  e = etot
  ef = efermi
end subroutine get_energies

! Return eigenvalues
subroutine get_eigenvalues(norbitals_new, nkpoints, eig)
  use iso_c_binding
  use M_system, only: eigen_k
  implicit none
  integer(c_long), intent(in) :: norbitals_new, nkpoints
  real(c_double), dimension(norbitals_new, nkpoints), intent(inout) :: eig
  eig = eigen_k(1:norbitals_new, :)
end subroutine get_eigenvalues

! Return partial charges and shell charges
subroutine get_charges(natoms, nsh_max, qpartial, qshell)
  use iso_c_binding
  use M_fdata, only : Qneutral, nssh
  use M_system, only : Qin, imass
  implicit none
  integer(c_long), intent(in) :: nsh_max, natoms
  real(c_double), dimension(natoms), intent(inout) :: qpartial
  real(c_double), dimension(nsh_max, natoms), intent(inout) :: qshell
  integer(c_long) :: iatom, issh, in1
  do iatom=1,natoms
    in1 = imass(iatom)
    qpartial(iatom) = 0.0d0
    do issh=1,nssh(in1)
      qshell(issh, iatom) = Qin(issh, iatom)
      qpartial(iatom) = qpartial(iatom) + Qneutral(issh,in1) - Qin(issh,iatom)
    end do
    do issh=nssh(in1)+1,nsh_max
      qshell(issh, iatom) = 0.0d0
    end do
  end do
end subroutine get_charges

! Get the forces in each atom
subroutine get_forces(natoms, forces)
  use iso_c_binding
  use M_system, only : ftot
  implicit none
  integer(c_long), intent(in) :: natoms
  real(c_double), dimension(3, natoms), intent(inout) :: forces
  forces = ftot
end subroutine get_forces
