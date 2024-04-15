subroutine set_igamma(aux)
  use M_system, only : igamma
  implicit none
  integer, intent(in):: aux
  igamma=aux
end
 
subroutine set_icluster(aux)
  use M_system, only : icluster
  implicit none
  integer, intent(in):: aux
  icluster=aux
end

subroutine set_idipole(aux)
  use M_system, only : idipole
  implicit none
  integer, intent(in):: aux
  idipole=aux
end

subroutine set_iqout(aux)
  use M_system, only : iqout
  implicit none
  integer, intent(in):: aux
  iqout=aux
end


integer function get_nssh(iaux)
 use M_fdata, only : nssh
 use M_system, only : imass
 get_nssh=nssh(imass(iaux))
 return
end

real function get_etot()
  use M_system, only : etot
  get_etot=etot
  return
end


real function get_atom_force(iaux,jaux)
  use M_system, only : ftot
  implicit none
  integer, intent(in):: iaux
  integer, intent(in):: jaux
  get_atom_force = ftot(jaux,iaux)
  return
end function get_atom_force

real function get_shell_atom_charge(iauxssh,iauxatom)
  use M_system, only : Qin
  implicit none
  integer, intent(in):: iauxatom
  integer, intent(in):: iauxssh
  get_shell_atom_charge = Qin(iauxssh,iauxatom)
end function get_shell_atom_charge

 
subroutine loadfdata_from_path(fdatafile)
  use M_fdata
  implicit none
  character(len=400),intent(in):: fdatafile
  fdatalocation=trim(fdatafile)
  call load_fdata()
end

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
end


subroutine print_atoms_positions()
  use M_system
  implicit none
  integer iatom
  do iatom = 1, natoms
    write (*,'(3x,a2, 3(2x,f10.5))') symbol(iatom), ratom(:,iatom)
  end do
end 


subroutine loadlvs_100()
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
end
 
subroutine loadlvs_from_file(lvsfile)
  use M_system
  implicit none
  character(len=400),intent(in):: lvsfile
  open (unit = 72, file = trim(lvsfile), status = 'old')
  read (72,*) a1vec(:)
  read (72,*) a2vec(:)
  read (72,*) a3vec(:)
  close(72)
 end                       
 
subroutine loadkpts_gamma()
  use M_system
  implicit none
  nkpoints = 1
  allocate (special_k(3, nkpoints))
  allocate (special_k_orig(3, nkpoints))
  allocate (scale_k(3, nkpoints))
  allocate (weight_k(nkpoints))
  allocate (weight_k_orig(nkpoints))
  special_k_orig(:,1) = 0
  weight_k_orig(1) = 1
  weight_k(1) = 1
  special_k(:,1) = special_k_orig(:,1)
  weight_k(1) = weight_k_orig(1)
end
 
subroutine loadkpts_from_file(kptsfile)
  use M_system
  implicit none
  character(len=400),intent(in):: kptsfile
  integer :: ikpoint
  real :: sum_weight
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
end
 
 
subroutine rescal_structure(rescal)
  use M_system
  implicit none
  real,intent(in)::rescal
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
end
 
subroutine call_allocate_system()
  use M_system
  implicit none
  call allocate_system()
end

 
subroutine call_scf_loop()
  use M_system
  implicit none
  call scf_loop ()
end

subroutine call_getenergy()
  use M_system
  implicit none
  call getenergy ()
end

subroutine call_getforces()
  use M_system
  implicit none
  call getforces()
end


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
end
 
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
end

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
end 
 
