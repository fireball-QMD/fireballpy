subroutine read_1c ()
  use M_fdata, only: nspecies, nzx, fdataLocation, onecfname, exc1c0, nuxc1c, &
    & dnuxc1c, dexc1c, nsh_max, nssh
  implicit none
  logical skip_it
  integer iline, Nlines_vdip1c_max, trash, in1, ins, issh, isorp, &
    & itype, jssh, kssh, kkssh, numsh, ideriv, iissh, jjssh
  integer, dimension (nsh_max) :: imask
  real*8, dimension (nspecies) :: idshell
  character (len=3) :: auxz
  character (len=1000) filename, root

  onecfname(1) = "xc1c_dqi "
  onecfname(2) = "nuxc1crho"
  onecfname(3) = "exc1crho "

  exc1c0 = 0.0d0
  nuxc1c = 0.0d0
  dexc1c = 0.0d0
  dnuxc1c = 0.0d0

  do in1 = 1, nspecies
    write (auxz,"(""."",i2.2)") nzx(in1)
    root = trim(fdataLocation)//trim(onecfname(1))
    filename = trim(root)//auxz//".dat"
    open (unit = 36, file = filename, status = "old")
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
    write (auxz,"(""."",i2.2)") nzx(in1)
    root = trim(fdataLocation)//trim(onecfname(2))
    filename = trim(root)//auxz//".dat"
    open (unit = 36, file = filename, status = "old")
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
    write (auxz,"(""."",i2.2)") nzx(in1)
    root = trim(fdataLocation)//trim(onecfname(3))
    filename = trim(root)//auxz//".dat"
    open (unit = 36, file = filename, status = "old")
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
