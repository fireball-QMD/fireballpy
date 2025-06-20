subroutine read_1c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nspecies, nzx, fdataLocation, onecfname, exc1c0, nuxc1c, &
    & dnuxc1c, dexc1c, d2exc1c, d2nuxc1c, nsh_max, nssh
  implicit none
  integer :: iline, in1, issh, itype, jssh, kssh, kkssh, numsh
  character (len=3) :: auxz
  character (len=1000) :: root

  if (allocated(exc1c0)) deallocate(exc1c0)
  if (allocated(nuxc1c)) deallocate(nuxc1c)
  if (allocated(dexc1c)) deallocate(dexc1c)
  if (allocated(d2exc1c)) deallocate(d2exc1c)
  if (allocated(dnuxc1c)) deallocate(dnuxc1c)
  if (allocated(d2nuxc1c)) deallocate(d2nuxc1c)
  allocate(exc1c0 (nspecies,nsh_max,nsh_max))
  allocate(nuxc1c (nspecies,nsh_max,nsh_max))
  allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
  allocate(d2exc1c (nspecies,nsh_max,nsh_max))
  allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
  allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))
  exc1c0 = 0.0d0
  nuxc1c = 0.0d0
  dexc1c = 0.0d0
  d2exc1c = 0.0d0
  dnuxc1c = 0.0d0
  d2nuxc1c = 0.0d0

  do in1 = 1, nspecies
    write (auxz,'(''.'',i2.2)') nzx(in1)
    root = trim(fdataLocation) // trim(onecfname(1)) // auxz // '.dat'
    open (unit = 36, file = root, status = 'old')
    do iline = 1, 6
      read (36,*)
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
  end do

  do in1 = 1, nspecies
    write (auxz,'(''.'',i2.2)') nzx(in1)
    root = trim(fdataLocation) // trim(onecfname(2)) // auxz // '.dat'
    open (unit = 36, file = root, status = 'old')
    do iline = 1, 6
      read (36,*)
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
    write (auxz,'(''.'',i2.2)') nzx(in1)
    root = trim(fdataLocation) // trim(onecfname(3)) // auxz // '.dat'
    open (unit = 36, file = root, status = 'old')
    do iline = 1, 6
      read (36,*)
    end do
    do kssh = 1, nssh(in1)
      read (36,*) itype, numsh, kkssh
      do issh = 1, numsh
        read (36,*) (dexc1c(in1,issh,jssh,kssh),jssh=1,numsh)
      end do
    end do
    close(36)
  end do
end subroutine read_1c
