subroutine set_igamma(aux)
  use M_system, only : igamma
  implicit none
  integer, intent(in):: aux
  igamma=aux
end subroutine set_igamma
 
integer function get_igamma()
  use M_system, only : igamma
  get_igamma=igamma
  return
end function get_igamma
 
subroutine set_icluster(aux)
  use M_system, only : icluster
  implicit none
  integer, intent(in):: aux
  icluster=aux
end subroutine set_icluster

integer function get_icluster()
  use M_system, only : icluster
  get_icluster=icluster
  return
end function get_icluster

subroutine set_ifixcharge(aux)
  use M_system, only : ifixcharge
  implicit none
  integer, intent(in):: aux
  ifixcharge=aux
end subroutine set_ifixcharge
 
subroutine set_idipole(aux)
  use M_system, only : idipole
  implicit none
  integer, intent(in):: aux
  idipole=aux
end subroutine set_idipole

subroutine set_iqout(aux)
  use M_system, only : iqout
  implicit none
  integer, intent(in):: aux
  iqout=aux
end subroutine set_iqout

integer function get_idipole()
  use M_system, only : idipole
  get_idipole=idipole
  return
end function get_idipole

integer function get_iqout()
  use M_system, only : iqout
  get_iqout=iqout
  return
end function get_iqout


integer function get_nssh(iaux)
  use M_fdata, only : nssh
  use M_system, only : imass
  implicit none
  integer, intent(in):: iaux
  get_nssh=nssh(imass(iaux))
  return
end function get_nssh

real*8 function get_etot()
  use M_system, only : etot
  get_etot=etot
  return
end function get_etot

integer function get_norbitals_new()
  use M_system, only : norbitals_new
  get_norbitals_new = norbitals_new
  return
end function get_norbitals_new

real*8 function get_eigen(iaux,jaux)
  use M_system, only : eigen_k !(imu,ikpoint)
  implicit none
  integer, intent(in):: iaux
  integer, intent(in):: jaux
  get_eigen = eigen_k(iaux,jaux)
  return
end function get_eigen

real*8 function get_atom_force(iaux,jaux)
  use M_system, only : ftot
  implicit none
  integer, intent(in):: iaux
  integer, intent(in):: jaux
  get_atom_force = ftot(jaux,iaux)
  return
end function get_atom_force

real*8 function get_shell_atom_charge(iauxssh,iauxatom)
  use M_system, only : Qin
  implicit none
  integer, intent(in):: iauxatom
  integer, intent(in):: iauxssh
  get_shell_atom_charge = Qin(iauxssh,iauxatom)
end function get_shell_atom_charge

integer function get_fdata_is_load()
  use M_fdata, only : fdata_is_load
  get_fdata_is_load = fdata_is_load
end function get_fdata_is_load


subroutine loadfdata_from_path(fdatafile)
  use M_fdata
  implicit none
  character(len=400),intent(in):: fdatafile
  fdatalocation=trim(fdatafile)
  fdata_is_load = 1
  call load_fdata()
end subroutine loadfdata_from_path
 
subroutine set_cell(lvs)
  use M_system
  implicit none
  real*8, dimension(3,3), intent(in) :: lvs
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
  use M_system
  use M_fdata, only : symbolA, nspecies, nzx 
  implicit none
  integer, intent(in) :: naux
  integer, dimension(naux), intent(in) :: z
  real*8, dimension(naux,3), intent(in) :: xyz
  integer :: iatom,ispec
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
  use M_system
  implicit none
  integer, intent(in) :: naux
  real*8, dimension(naux,3), intent(in) :: xyz
  integer iatom
  real*8 distance

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
  use M_system
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
  use M_system
  implicit none
  integer, intent(in) :: naux
  real*8, dimension(naux,3), intent(in) :: kpts
  integer :: ikpoint
  real*8 :: sum_weight
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
  if (allocated(ioccupy_k)) deallocate( ioccupy_k)
  allocate(ioccupy_k (norbitals, nkpoints))
  if (allocated(foccupy)) deallocate(foccupy) 
  allocate(foccupy (norbitals, nkpoints))
  if (allocated(eigen_k)) deallocate(eigen_k) 
  allocate (eigen_k (norbitals, nkpoints))
  if (allocated(bbnkre_o)) deallocate(bbnkre_o)
  allocate (bbnkre_o(norbitals,norbitals,nkpoints))

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
  use M_system
  implicit none
  call allocate_system()
end subroutine call_allocate_system

 
subroutine call_scf_loop()
  use M_system
  implicit none
  call scf_loop ()
end subroutine call_scf_loop

subroutine call_getenergy()
  use M_system
  implicit none
  call getenergy ()
end subroutine call_getenergy

subroutine call_getforces()
  use M_system
  implicit none
  call getforces()
end subroutine call_getforces

