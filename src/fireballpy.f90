program pyreball

  use M_fdata
  use M_system
  implicit none
  integer :: ispec
  integer :: issh  
  integer :: iatom 
  integer :: nucz
  integer :: ikpoint
  integer :: in1
  logical :: zindata
  integer :: iostat
  character(len=400) :: archivo
     
  print*, '- Fireballpy is a minimal version of the fireball program.'
  print*, '  itheory     = 1 !FIX DOGS' 
  print*, '  itheory_xc  = 2 !FIX McWEDA'
  write(*,'(3x,a12,a1,i2)') 'icluster      ','=',icluster
  write(*,'(3x,a12,a1,i2)') 'iforce        ','=',iforce
  write(*,'(3x,a12,a1,i2)') 'idipole       ','=',idipole
  write(*,'(3x,a12,a1,i2)') 'iqout         ','=',iqout   
  print*,''

  fdatalocation='/home/dani/FB/git/create/coutput'
  fdatalocation='/home/dani/Fdata_HC-new/'
   
 
  archivo=trim(fdataLocation)//'/info.dat'
  inquire(file=archivo, exist=iostat)
  if (iostat .ne.0) then
    call load_fdata()
    print*,'- Read fdata from : ',trim(archivo)  
    do ispec = 1, nspecies
      write (*,'(a,i2,a,a2,a,i2,a,i2)') '   spec = ',ispec,'; ele = ',symbolA(ispec),'; Z = ',nzx(ispec), '; nssh = ',nssh(ispec)        
    end do
  else
    print*,'Problemas con :',trim(archivo)
    stop
  end if
  print*,''
  

  archivo='input.bas'
  inquire(file=archivo, exist=iostat)
  if (iostat .ne.0) then
    print*,'- Read atoms positions from :',trim(archivo)
    open (unit = 69, file = 'input.bas', status = 'old')
    read (69, *) natoms
    allocate (ratom (3, natoms))
    allocate (symbol (natoms))
    allocate (imass (natoms))
    !write(*,'(i5)') natoms
    !print*,''
    do iatom = 1, natoms
     read (69,*) nucz,ratom(:,iatom)
     zindata = .false.
     do ispec = 1, nspecies
     if ( nucz .eq. nzx(ispec)) then
       zindata = .true.
       imass(iatom) = ispec
       symbol(iatom)=symbolA(ispec)
       write (*,'(3x,a2, 3(2x,f10.5))') symbol(iatom), ratom(:,iatom)
     end if
     end do
    end do
    close (unit = 69)
   else
    print*,'Problemas con :',trim(archivo)
    stop
   end if
   print*,''


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

  write(*,'(A,F20.6,A,I4,A,F12.10,A,L1)') 'EBS = ',ebs,'; Kscf =',Kscf,'; sigma =',sigma,'; scf_achieved =',scf_achieved
 
  write(*,*)'========== CHARGES ====== '
  do iatom = 1, natoms
    in1 = imass(iatom)
     write (*,'(2x, 10f14.8)') (Qin(issh,iatom), issh = 1, nssh(in1))
  end do
 

  call getenergy ()
  write (*,*) ' ---------- T H E  T O T A L  E N E R G Y ----------- '
  write (*,*) '  '
  write (*,502) ebs
  write (*,503) uiiuee
  write (*,504) etotxc_1c
  write (*,505) uxcdcc
  write (*,507) etot
  write (*,508) etotper
  write (*,509) atomic_energy
  write (*,510) etot - atomic_energy
  write (*,512) efermi
  write (*,*) '  '
  write (*,511) (etot - atomic_energy)/natoms
  write (*,*) ' ----------------------------------------------------- '




  if (iforce .eq. 1) then
    call getforces()
    write (*,*) '====================== '
    write (*,*) ' The grand total force (eV/A): '
    do iatom = 1, natoms
      write (*,130)  iatom, ftot(:,iatom)
    end do
  end if

130     format (2x, ' iatom = ', i4, ' ftot      = ', 3e14.6)
100     format (2x, 70('='))
500     format (2x, ' Time step = ', i6, ' SCF step = ', i3)
501     format (2x, ' Time step = ', i6)
502     format (2x, '           ebs = ', f15.6)
503     format (2x, '     uii - uee = ', f15.6)
504     format (2x, '     etotxc_1c = ', f15.6)
505     format (2x, '        uxcdcc = ', f15.6)
507     format (2x, '          ETOT = ', f15.6)
508     format (2x, '     Etot/atom = ', f15.6)
509     format (2x, ' Atomic Energy = ', f15.6)
510     format (2x, '     CohesiveE = ', f15.6)
511     format (2x, ' Cohesive Energy per atom  = ', f15.6)
512     format (2x, '   Fermi Level = ', f15.6)

end program pyreball

