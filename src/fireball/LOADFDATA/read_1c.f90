subroutine read_1c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nspecies, nzx, fdataLocation, onecfname, exc_1c_0, vxc_1c_0, &
    &  gxc_1c, fxc_1c, nsh_max, nssh 
  implicit none
  integer :: iline, in1, issh, itype, jssh, kssh, kkssh, numsh
  character (len=3) :: auxz
  character (len=1000) :: root

  if (allocated(exc_1c_0)) deallocate(exc_1c_0)
  if (allocated(vxc_1c_0)) deallocate(vxc_1c_0)
  if (allocated(gxc_1c)) deallocate(gxc_1c)
  if (allocated(fxc_1c)) deallocate(fxc_1c)
  allocate(exc_1c_0 (nspecies,nsh_max))
  allocate(gxc_1c (nspecies,nsh_max,nsh_max,nsh_max))
  allocate(vxc_1c_0 (nspecies,nsh_max,nsh_max))
  allocate(fxc_1c (nspecies,nsh_max,nsh_max))
  exc_1c_0 = 0.0d0
  vxc_1c_0 = 0.0d0
  gxc_1c = 0.0d0
  fxc_1c = 0.0d0

  do in1 = 1, nspecies
    write (auxz,'(''.'',i2.2)') nzx(in1)
    ! leemos onecenter_xc.06.dat 
    root = trim(fdataLocation) // trim(onecfname(1)) // auxz // '.dat'
    open (unit = 36, file = root, status = 'old')
    read (36,*) numsh
    do issh = 1, numsh
      read (36,*) (vxc_1c_0(in1,issh,jssh),jssh=1,numsh)
    end do
    read (36,*)
    do kssh= 1, numsh
      do issh = 1, numsh
        read (36,*) (gxc_1c(in1,issh,jssh,kssh),jssh=1,numsh)
      end do
      read (36,*)
    end do
    read (36,*) (exc_1c_0(in1,jssh), jssh = 1, numsh)
    read (36,*)
    do kssh = 1, numsh
      read (36,*) (fxc_1c(in1,jssh,kssh),jssh=1,numsh)
      read (36,*)
    end do
    read (36,*)
    close (unit = 36)
  end do
end subroutine read_1c
