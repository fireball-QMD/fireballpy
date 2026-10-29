subroutine read_2c (interaction)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use :: constants, only: FNAME_LEN
  use :: utils, only: utils_itos, utils_open
  use :: math, only: math_interp_new
  use M_fdata, only: fdataLocation, twocfname, nspecies, nssh, nsshPP, nzx, &
    & ind2c, index_max2c, index_maxS, index_maxPP, index_max2cDipX, index_max2cDipY, nsh_max, initype, splineint_2c, &
    & TWOCENTER_VNA_L, TWOCENTER_VNA_R, TWOCENTER_VNA_A, TWOCENTER_VXC_0, TWOCENTER_VXC_L, TWOCENTER_VXC_R, TWOCENTER_DEN_L, &
    & TWOCENTER_DEN_R, TWOCENTER_DEN_A, TWOCENTER_DENS_L, TWOCENTER_DENS_R, TWOCENTER_DENS_A, TWOCENTER_VNL, &
    & TWOCENTER_COULOMB, TWOCENTER_OVERLAPS, TWOCENTER_DIP_Y, TWOCENTER_DIP_X
  implicit none
  integer, intent (in) :: interaction
  integer :: in1, in2, iounit, isorp, issh, isub2c, itype, numz, num_nonzero, header, ipoint, integral
  character (len = FNAME_LEN) :: filename, root
  real(double), allocatable :: zarr(:), xintegral_2c(:,:)

  select case (interaction)
  case (TWOCENTER_VNL)
    header = 4
  case (TWOCENTER_VXC_0, TWOCENTER_VXC_L, TWOCENTER_VXC_R)
    header = 5
  case default
    header = 2
  end select
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      select case (interaction)
      case (TWOCENTER_VNA_L, TWOCENTER_VXC_L, TWOCENTER_DEN_L, TWOCENTER_DENS_L)
        isub2c = nssh(in1)
      case (TWOCENTER_VNA_R, TWOCENTER_VNA_A, TWOCENTER_VXC_R, TWOCENTER_DEN_R, &
          TWOCENTER_DEN_A, TWOCENTER_DENS_R, TWOCENTER_DENS_A)
        isub2c = nssh(in2)
      case default
        isub2c = 0
      end select
      do isorp = initype(interaction), isub2c
        root = trim(fdataLocation)//trim(twocfname(interaction))
        if (isub2c > 0) then
          root = trim(root)//"_"//trim(utils_itos(isorp, lpad=2))
        end if
        filename = trim(root)//"."//trim(utils_itos(nzx(in1), lpad=2))//"."//trim(utils_itos(nzx(in2), lpad=2))//".dat"
        itype = ind2c(interaction, isorp)
        select case (interaction)
        case (TWOCENTER_VNA_A, TWOCENTER_DEN_A)
          num_nonzero = index_max2c(in1, in1)
        case (TWOCENTER_VNL)
          num_nonzero = index_maxPP(in1, in2)
        case (TWOCENTER_DIP_X)
          num_nonzero = index_max2cDipX(in1, in2)
        case (TWOCENTER_DIP_Y)
          num_nonzero = index_max2cDipY(in1, in2)
        case (TWOCENTER_DENS_A)
          num_nonzero = index_maxS(in1, in1)
        case (TWOCENTER_DENS_L, TWOCENTER_DENS_R, TWOCENTER_OVERLAPS)
          num_nonzero = index_maxS(in1, in2)
        case (TWOCENTER_COULOMB)
          num_nonzero = nssh(in1)*nssh(in2)
        case default
          num_nonzero = index_max2c(in1, in2)
        end select
        if (num_nonzero == 0) cycle
        iounit = utils_open(filename, "r", skip=header)
        read (iounit,"(6x,i8)") numz
        read (iounit,*)
        if (allocated(zarr)) deallocate(zarr)
        if (allocated(xintegral_2c)) deallocate(xintegral_2c)
        allocate(zarr(numz), xintegral_2c(num_nonzero, numz))
        do ipoint = 1, numz
          read (iounit,*) zarr(ipoint), (xintegral_2c(integral,ipoint), integral = 1, num_nonzero)
        end do
        do integral = 1, num_nonzero
          splineint_2c(integral,itype,in1,in2) = math_interp_new(zarr, xintegral_2c(integral,:), rval=0.0d0)
        end do
        close(iounit)
      end do
    end do
  end do
end subroutine read_2c
