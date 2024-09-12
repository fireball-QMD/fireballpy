subroutine get_errno(errno_out)
  use iso_c_binding
  use M_system, only: errno
  implicit none
  integer(c_long), intent(out) :: errno_out
  errno_out = errno
end subroutine get_errno

subroutine set_igamma(aux)
  use iso_c_binding
  use M_system, only : igamma
  implicit none
  integer(c_long), intent(in):: aux
  igamma=aux
end subroutine set_igamma
 
subroutine get_igamma(igamma_out)
  use iso_c_binding
  use M_system, only : igamma
  implicit none
  integer(c_long), intent(out) :: igamma_out
  igamma_out=igamma
end subroutine get_igamma
 
subroutine set_icluster(aux)
  use iso_c_binding
  use M_system, only : icluster
  implicit none
  integer(c_long), intent(in):: aux
  icluster=aux
end subroutine set_icluster

subroutine get_icluster(icluster_out)
  use iso_c_binding
  use M_system, only : icluster
  implicit none
  integer(c_long), intent(out) :: icluster_out
  icluster_out=icluster
end subroutine get_icluster

subroutine set_ifixcharge(aux)
  use iso_c_binding
  use M_system, only : ifixcharge
  implicit none
  integer(c_long), intent(in):: aux
  ifixcharge=aux
end subroutine set_ifixcharge
 
subroutine set_idipole(aux)
  use iso_c_binding
  use M_system, only : idipole
  implicit none
  integer(c_long), intent(in):: aux
  idipole=aux
end subroutine set_idipole

subroutine set_iqout(aux)
  use iso_c_binding
  use M_system, only : iqout
  implicit none
  integer(c_long), intent(in):: aux
  iqout=aux
end subroutine set_iqout

subroutine get_idipole(idipole_out)
  use iso_c_binding
  use M_system, only : idipole
  implicit none
  integer(c_long), intent(out) :: idipole_out
  idipole_out=idipole
end subroutine get_idipole

subroutine get_iqout(iqout_out)
  use iso_c_binding
  use M_system, only : iqout
  implicit none
  integer(c_long), intent(out) :: iqout_out
  iqout_out=iqout
end subroutine get_iqout

subroutine get_nssh(iaux, nssh_out)
  use iso_c_binding
  use M_fdata, only : nssh
  use M_system, only : imass
  implicit none
  integer(c_long), intent(in):: iaux
  integer(c_long), intent(out) :: nssh_out
  nssh_out=nssh(imass(iaux))
end subroutine get_nssh

subroutine get_etot(etot_out)
  use iso_c_binding
  use M_system, only : etot
  implicit none
  real(c_double), intent(out) :: etot_out
  etot_out=etot
end subroutine get_etot

subroutine get_efermi(efermi_out)
  use iso_c_binding
  use M_system, only : efermi
  implicit none
  real(c_double), intent(out) :: efermi_out
  efermi_out=efermi
end subroutine get_efermi

subroutine get_norbitals_new(norbitals_new_out)
  use iso_c_binding
  use M_system, only : norbitals_new
  implicit none
  integer(c_long), intent(out) :: norbitals_new_out
  norbitals_new_out = norbitals_new
end subroutine get_norbitals_new

subroutine get_eigen(iaux,jaux,eigen_out)
  use iso_c_binding
  use M_system, only : eigen_k !(imu,ikpoint)
  implicit none
  integer(c_long), intent(in):: iaux
  integer(c_long), intent(in):: jaux
  real(c_double), intent(out) :: eigen_out
  eigen_out = eigen_k(iaux,jaux)
end subroutine get_eigen

subroutine get_atom_force(iaux,jaux,atom_force_out)
  use iso_c_binding
  use M_system, only : ftot
  implicit none
  integer(c_long), intent(in):: iaux
  integer(c_long), intent(in):: jaux
  real(c_double), intent(out) :: atom_force_out
  atom_force_out = ftot(jaux,iaux)
end subroutine get_atom_force

subroutine get_shell_atom_charge(iauxssh,iauxatom,shell_atom_charge_out)
  use iso_c_binding
  use M_system, only : Qin
  implicit none
  integer(c_long), intent(in):: iauxatom
  integer(c_long), intent(in):: iauxssh
  real(c_double), intent(out) :: shell_atom_charge_out
  shell_atom_charge_out = Qin(iauxssh,iauxatom)
end subroutine get_shell_atom_charge
 
subroutine get_partial_charge(iatom,nssh,partial_charge_out)
  use iso_c_binding
  use M_system, only : Qin, Q0_TOT
  implicit none
  integer(c_long) issh
  integer(c_long), intent(in):: iatom, nssh
  real(c_double), intent(out) :: partial_charge_out
  partial_charge_out = Q0_TOT(iatom)
  do  issh = 1, nssh
    partial_charge_out = partial_charge_out - Qin(issh,iatom)
  end do
end subroutine get_partial_charge

subroutine set_shell_atom_charge(natoms,nssh,qaux)
  use iso_c_binding
  use M_system, only : Qin
  implicit none
  integer(c_long), intent(in):: natoms, nssh
  real(c_double), dimension(natoms, nssh), intent(in):: qaux
  integer(c_long) :: iauxssh, iauxatom
  do iauxatom=1,natoms
    do iauxssh=1,nssh
      Qin(iauxssh, iauxatom) = qaux(iauxatom, iauxssh)
    end do
  end do
end subroutine set_shell_atom_charge

subroutine get_fdata_is_load(aux)
  use iso_c_binding
  use M_fdata, only : fdata_is_load
  implicit none
  integer(c_long), intent(out) :: aux
  aux = fdata_is_load
end subroutine get_fdata_is_load

subroutine loadfdata_from_path(fdatafile)
  use iso_c_binding
  use M_fdata, only: fdataLocation, fdata_is_load
  implicit none
  character(len=400),intent(in):: fdatafile
  fdatalocation=trim(fdatafile)
  fdata_is_load = 1
  call load_fdata()
end subroutine loadfdata_from_path
 
subroutine set_cell(lvs)
  use iso_c_binding
  use M_system, only: a1vec, a2vec, a3vec
  implicit none
  real(c_double), dimension(3,3), intent(in) :: lvs
  a1vec(1) = lvs(1,1)
  a1vec(2) = lvs(1,2)
  a1vec(3) = lvs(1,3)
  a2vec(1) = lvs(2,1)
  a2vec(2) = lvs(2,2)
  a2vec(3) = lvs(2,3)
  a3vec(1) = lvs(3,1)
  a3vec(2) = lvs(3,2)
  a3vec(3) = lvs(3,3)
end subroutine set_cell

subroutine set_coords(naux, z, xyz)
  use iso_c_binding
  use M_system, only: natoms, ratom, symbol, imass
  use M_fdata, only : symbolA, nspecies, nzx 
  implicit none
  integer(c_long), intent(in) :: naux
  integer(c_long), dimension(naux), intent(in) :: z
  real(c_double), dimension(naux,3), intent(in) :: xyz
  integer(c_long) :: iatom,ispec
  natoms=naux
  if (allocated(ratom)) deallocate(ratom)
  if (allocated(symbol)) deallocate(symbol)
  if (allocated(imass)) deallocate(imass)
  allocate (ratom (3, natoms))
  allocate (symbol (natoms))
  allocate (imass (natoms))
  do iatom = 1, natoms
    ratom(:, iatom) = xyz(iatom, :)
    do ispec = 1, nspecies
      if ( z(iatom) .eq. nzx(ispec)) then
        imass(iatom) = ispec
        symbol(iatom) = symbolA(ispec)
      end if
   end do
  end do
end subroutine set_coords
 
subroutine set_coords_xyz(naux,xyz)
  use iso_c_binding
  use M_system, only: natoms, ratom, ishiftO, shifter
  implicit none
  integer(c_long), intent(in) :: naux
  real(c_double), dimension(naux,3), intent(in) :: xyz
  integer(c_long) iatom
  real(c_double) :: distance

  do iatom = 1, natoms
    ratom(:, iatom) = xyz(iatom, :)
  end do

  ishiftO = 0
  do iatom = 1, natoms
    distance = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2
    distance = sqrt(distance)
    if (distance .lt. 1.0d-4) ishiftO = 1
  end do

  if (ishiftO .eq. 1) then
    do iatom = 1, natoms
      ratom(:,iatom) = ratom(:,iatom) + shifter
    end do
  end if
end subroutine set_coords_xyz

subroutine load_cell_100()
  use iso_c_binding
  use M_system, only: a1vec, a2vec, a3vec
  implicit none
  a1vec(1) = 100
  a1vec(2) = 0
  a1vec(3) = 0
  a2vec(1) = 0
  a2vec(2) = 100
  a2vec(3) = 0
  a3vec(1) = 0
  a3vec(2) = 0
  a3vec(3) = 100
end subroutine load_cell_100
 
subroutine set_kpoints(naux,kpts)
  use iso_c_binding
  use M_system, only: nkpoints, special_k, special_k_orig, scale_k, weight_k, weight_k_orig
  implicit none
  integer(c_long), intent(in) :: naux
  real(c_double), dimension(naux,3), intent(in) :: kpts
  integer(c_long) :: ikpoint
  real(c_double) :: sum_weight
  nkpoints=naux
  if (allocated(special_k)) deallocate(special_k)
  if (allocated(special_k_orig)) deallocate(special_k_orig)
  if (allocated(scale_k)) deallocate(scale_k)
  if (allocated(weight_k)) deallocate(weight_k)
  if (allocated(weight_k_orig)) deallocate(weight_k_orig)
  allocate (special_k(3, nkpoints))
  allocate (special_k_orig(3, nkpoints))
  allocate (scale_k(3, nkpoints))
  allocate (weight_k(nkpoints))
  allocate (weight_k_orig(nkpoints))
  sum_weight = 0.0d0
  do ikpoint = 1, nkpoints
    special_k_orig(:,ikpoint) = kpts(ikpoint, :)
    weight_k_orig(ikpoint) = 1.0d0/nkpoints
    sum_weight = sum_weight + weight_k_orig(ikpoint)
  end do
  do ikpoint = 1, nkpoints
    special_k(:,ikpoint) = special_k_orig(:,ikpoint)
    weight_k(ikpoint) = weight_k_orig(ikpoint)
  end do
end subroutine set_kpoints

subroutine call_allocate_system()
  use iso_c_binding
  implicit none
  call allocate_system()
end subroutine call_allocate_system

subroutine call_scf_loop(verbose)
  use iso_c_binding
  implicit none
  logical, intent(in) :: verbose
  call scf_loop (verbose)
end subroutine call_scf_loop

subroutine call_getforces()
  use iso_c_binding
  implicit none
  call getforces()
end subroutine call_getforces

subroutine set_mixer_params(iam, msi, imix, bbmix, ww02, stol, wis)
  use iso_c_binding
  use M_system, only: wi, ialgmix, max_scf_iterations, idmix, w02, bmix, sigmatol
  implicit none
  integer(c_long), intent(in) :: iam
  integer(c_long), intent(in) :: msi
  integer(c_long), intent(in) :: imix
  real(c_double), intent(in) :: ww02
  real(c_double), intent(in) :: bbmix
  real(c_double), intent(in) :: stol
  real(c_double), intent(in), dimension(msi) :: wis
  allocate(wi(msi))
  ialgmix = iam
  max_scf_iterations = msi
  idmix = imix
  w02 = ww02
  bmix = bbmix
  sigmatol = stol
  wi(:) = wis(:)
end subroutine set_mixer_params
