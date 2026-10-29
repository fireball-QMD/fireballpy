subroutine read_1c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use :: constants, only: FNAME_LEN
  use :: utils, only: utils_itos, utils_open
  use M_fdata, only: nspecies, nzx, fdataLocation, onecfname, exc_1c_0, vxc_1c_0, &
    &  gxc_1c, fxc_1c, nsh_max, nssh, Qref
  implicit none
  integer :: iline, io, in1, issh, itype, jssh, kssh, numsh
  character (len=FNAME_LEN) :: root
  real(double), allocatable :: temp(:,:)

  if (allocated(exc_1c_0)) deallocate(exc_1c_0)
  if (allocated(vxc_1c_0)) deallocate(vxc_1c_0)
  if (allocated(gxc_1c)) deallocate(gxc_1c)
  if (allocated(fxc_1c)) deallocate(fxc_1c)

  allocate(exc_1c_0 (nsh_max,nspecies))
  allocate(fxc_1c (nsh_max,nsh_max,nspecies))

  allocate(gxc_1c (nsh_max,nsh_max,nsh_max,nspecies))
  allocate(vxc_1c_0 (nsh_max,nsh_max,nspecies))

  exc_1c_0 = 0.0d0
  vxc_1c_0 = 0.0d0
  gxc_1c = 0.0d0
  fxc_1c = 0.0d0

  do in1 = 1, nspecies
    numsh = nssh(in1)
    root = trim(fdataLocation)//trim(onecfname(1))//"."//trim(utils_itos(nzx(in1), lpad=2))//".dat"
    io = utils_open(root, "r", skip=5)
    allocate(temp(numsh, numsh))
    do issh = 1, numsh
      read (io,*) (vxc_1c_0(jssh,issh,in1),jssh=1,numsh)
    end do
    do kssh= 1, numsh
      read (io,*)
      do issh = 1, numsh
        read (io,*) (gxc_1c(jssh,issh,kssh,in1),jssh=1,numsh)
      end do
    end do
    read (io,*)
    do issh = 1, numsh
      read (io,*) (temp(jssh,issh),jssh=1,numsh)
    end do
    do issh = 1, numsh
      exc_1c_0(issh, in1) = temp(issh,issh)
    end do
    do kssh= 1, numsh
      read (io,*)
      do issh = 1, numsh
        read (io,*) (temp(jssh,issh),jssh=1,numsh)
      end do
      do issh = 1, numsh
        fxc_1c(issh,kssh,in1) = temp(issh,issh)
      end do
    end do
    deallocate(temp)
    close (unit = io)
  end do
end subroutine read_1c
