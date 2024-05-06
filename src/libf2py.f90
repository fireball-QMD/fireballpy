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

 

subroutine print_atoms_positions()
  use M_system
  implicit none
  integer iatom
  do iatom = 1, natoms
    write (*,'(3x,a2, 3(2x,f10.5))') symbol(iatom), ratom(:,iatom)
  end do
end subroutine print_atoms_positions


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
 
subroutine loadlvs_from_file(lvsfile)
  use M_system
  implicit none
  character(len=400),intent(in):: lvsfile
  open (unit = 72, file = trim(lvsfile), status = 'old')
  read (72,*) a1vec(:)
  read (72,*) a2vec(:)
  read (72,*) a3vec(:)
  close(72)
 end subroutine loadlvs_from_file

subroutine set_kpoints(naux,kpts)
  use M_system
  implicit none
  integer, intent(in) :: naux
  real*8, dimension(naux,3), intent(in) :: kpts
  integer :: ikpoint
  real*8 :: sum_weight
  nkpoints=naux
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

subroutine loadkpts_from_file(kptsfile)
  use M_system
  implicit none
  character(len=400),intent(in):: kptsfile
  integer :: ikpoint
  real*8 :: sum_weight
  open (unit = 54, file = kptsfile, status = 'old')
  read (54,*) nkpoints
  allocate (special_k(3, nkpoints))
  allocate (special_k_orig(3, nkpoints))
  allocate (scale_k(3, nkpoints))
  allocate (weight_k(nkpoints))
  allocate (weight_k_orig(nkpoints))
  sum_weight = 0.0d0
  do ikpoint = 1, nkpoints
    read (54,*) special_k_orig(:,ikpoint), weight_k_orig(ikpoint)
    sum_weight = sum_weight + weight_k_orig(ikpoint)
  end do
  close (unit = 54)
  do ikpoint = 1, nkpoints
    special_k(:,ikpoint) = special_k_orig(:,ikpoint)
    weight_k(ikpoint) = weight_k_orig(ikpoint)
  end do
end subroutine loadkpts_from_file

subroutine rescal_structure(rescal)
  use M_system
  implicit none
  real*8, intent(in)::rescal
  integer :: iatom,ikpoint
  a1vec(:)=a1vec(:)*rescal
  a2vec(:)=a2vec(:)*rescal
  a3vec(:)=a3vec(:)*rescal
  do iatom = 1, natoms
    ratom(:,iatom)=ratom(:,iatom)*rescal
  end do
  do ikpoint = 1, nkpoints
    special_k_orig(:,ikpoint)=special_k_orig(:,ikpoint)/rescal
  end do

  do ikpoint = 1, nkpoints
    special_k(:,ikpoint) = special_k_orig(:,ikpoint)
    weight_k(ikpoint) = weight_k_orig(ikpoint)
  end do
  write(*,'(3x,a12,a1,F6.3)') 'rescal         ','=',rescal
end subroutine rescal_structure
 
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

subroutine info_fdata()
  use M_fdata
  use M_system
  implicit none
  integer :: ispec
  print*, '- Fireballpy is a minimal version of the fireball program.'
  print*, '  itheory     = 1 !FIX DOGS'
  print*, '  itheory_xc  = 2 !FIX McWEDA'
  write(*,'(3x,a12,a1,i2)') 'icluster      ','=',icluster
  write(*,'(3x,a12,a1,i2)') 'iforce        ','=',iforce
  write(*,'(3x,a12,a1,i2)') 'idipole       ','=',idipole
  write(*,'(3x,a12,a1,i2)') 'iqout         ','=',iqout
  write(*,'(3x,a12,a1,i2)') 'igamma         ','=',igamma
  do ispec = 1, nspecies
    write (*,'(a,i2,a,a2,a,i2,a,i2)') '   spec = ',ispec,'; ele = ',symbolA(ispec),'; Z = ',nzx(ispec), '; nssh = ',nssh(ispec)  
  end do
end subroutine info_fdata
 
subroutine info_energy(out_energy)
  use M_system
  implicit none
  real*8, intent(out) :: out_energy
  write(*,'(3x,A,I4,A,F12.10,A,L1)') 'Kscf =',Kscf,'; sigma =',sigma,'; scf_achieved =',scf_achieved

  write (*,*) ' ---------- T H E  T O T A L  E N E R G Y ----------- '
  write (*,*) '  '
  write (*,'(2x,A, f15.6)') '           ebs = ',ebs
  write (*,'(2x,A, f15.6)') '     uii - uee = ',uiiuee
  write (*,'(2x,A, f15.6)') '     etotxc_1c = ',etotxc_1c
  write (*,'(2x,A, f15.6)') '     etotxc_1c = ',uxcdcc
  write (*,'(2x,A, f15.6)') '          ETOT = ',etot
  write (*,'(2x,A, f15.6)') '     Etot/atom = ',etotper
  write (*,'(2x,A, f15.6)') ' Atomic Energy = ',atomic_energy
  write (*,'(2x,A, f15.6)') '     CohesiveE = ',etot - atomic_energy
  write (*,'(2x,A, f15.6)') '   Fermi Level = ',efermi
  write (*,*) '  '
  write (*,'(2x,A, f15.6)')' Cohesive Energy per atom  = ', (etot - atomic_energy)/natoms
  write (*,*) ' ----------------------------------------------------- '

  out_energy = etot
end subroutine info_energy

subroutine info_forces(naux, out_forces)
  use M_system
  implicit none
  integer, intent(in) :: naux
  real*8, dimension(naux, 3), intent(out) :: out_forces
  integer :: iatom
  write (*,*) ' The grand total force (eV/A): '
  do iatom = 1, natoms
    out_forces(iatom, :) = ftot(:, iatom)
    write (*,'(2x,A,i4, A, 3e14.6)') ' iatom = ', iatom, ' ftot      = ',ftot(:,iatom)
  end do
end subroutine info_forces

