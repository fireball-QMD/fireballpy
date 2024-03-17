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

  call getenergy ()

  if (iforce .eq. 1) call getforces()

end program pyreball

