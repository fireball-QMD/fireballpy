module M_pairzero
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private
  public :: zero_pair_orb, zero_pair_sh

contains

  subroutine zero_pair_orb(a)
    use M_system, only: natoms, imass, neighn, neigh_j
    use M_fdata, only: num_orb
    real(dp), intent(inout) :: a(:,:,:,:,:)
    integer :: iatom, ineigh, in1, in2, inu
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
        in2 = imass(neigh_j(ineigh,iatom))
        do inu = 1, num_orb(in2)
          a(:, 1:num_orb(in1), inu, ineigh, iatom) = 0.0_dp
        end do
      end do
    end do
  end subroutine zero_pair_orb

  subroutine zero_pair_sh(a)
    use M_system, only: natoms, imass, neighn, neigh_j
    use M_fdata, only: nssh
    real(dp), intent(inout) :: a(:,:,:,:,:)
    integer :: iatom, ineigh, in1, in2, jssh
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
        in2 = imass(neigh_j(ineigh,iatom))
        do jssh = 1, nssh(in2)
          a(:, 1:nssh(in1), jssh, ineigh, iatom) = 0.0_dp
        end do
      end do
    end do
  end subroutine zero_pair_sh

end module M_pairzero
