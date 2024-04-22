subroutine read_1c ()
  use M_fdata
  implicit none
  integer iline
  integer Nlines_vdip1c_max
  integer trash
  integer in1
  integer ins
  integer issh
  integer isorp
  integer itype
  integer jssh
  integer kssh
  integer kkssh
  integer numsh
  integer, dimension (nsh_max) :: imask
  integer ideriv
  integer iissh, jjssh

  real(8), dimension (nspecies) :: idshell

  logical skip_it

  character (len=70) message
  character (len = 200) extension
  character (len = 200) filename
  character (len = 200) root
  character(2) :: auxz  

  allocate(exc1c0 (nspecies,nsh_max,nsh_max))
  allocate(nuxc1c (nspecies,nsh_max,nsh_max))
  allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
  allocate(d2exc1c (nspecies,nsh_max,nsh_max))
  allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
  allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))

  exc1c0=0.0d0
  nuxc1c=0.0d0
  dexc1c=0.0d0
  d2exc1c=0.0d0
  dnuxc1c=0.0d0
  d2nuxc1c=0.0d0

  do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/xc1c_dqi.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
      read (36,100) message
    end do
    read (36,*) itype, numsh
    do issh = 1, numsh
      read (36,*) (exc1c0(in1,issh,jssh),jssh=1,numsh)
    end do
    read (36,*)
    do issh = 1, numsh
      read (36,*) (nuxc1c(in1,issh,jssh),jssh=1,numsh)
    end do
    close(36)
  end do !in1 .. nspecies

  do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/nuxc1crho.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
      read (36,100) message
    end do
    do kssh = 1, nssh(in1)
      read (36,*) itype, numsh, kkssh
      do issh = 1, numsh
        read (36,*) (dnuxc1c(in1,issh,jssh,kssh),jssh=1,numsh)
      end do
    end do 
    close(36)
  end do 

  do in1 = 1, nspecies
    write (auxz,'(i2.2)') nzx(in1)
    root = trim(fdataLocation)//'/exc1crho.'//auxz//'.dat'
    open (unit = 36, file = root, status = 'unknown')
    do iline = 1, 6
      read (36,100) message
    end do
    do kssh = 1, nssh(in1)
      read (36,*) itype, numsh, kkssh
      do issh = 1, numsh
        read (36,*) (dexc1c(in1,issh,jssh,kssh),jssh=1,numsh)
      end do
    end do
    close(36)
  end do 

100   format (a70)
end subroutine read_1c

