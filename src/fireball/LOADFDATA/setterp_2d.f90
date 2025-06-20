subroutine setterp_2d ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: nspecies, icon3c, isorpmax, isorpmax_xc, hx_bcna, hy_bcna, hx_den3, hy_den3, &
    & x3cmax_bcna, y3cmax_bcna, x3cmax_den3, y3cmax_den3, numx3c_bcna, numy3c_bcna, numx3c_den3, numy3c_den3
  implicit none
  integer :: in1, in2, in3, index, isorp
  real(double) :: xmin, ymin

  xmin = 0.0d0
  ymin = 0.0d0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      do in3 = 1, nspecies
        index = icon3c(in1,in2,in3)
        do isorp = 0, isorpmax
          hx_bcna(isorp,index) = (x3cmax_bcna(isorp,index) - xmin)/real(numx3c_bcna(isorp,index) - 1)
          hy_bcna(isorp,index) = (y3cmax_bcna(isorp,index) - ymin)/real(numy3c_bcna(isorp,index) - 1)
        end do
          hx_den3(0,index) = 0.0d0
          hy_den3(0,index) = 0.0d0
        do isorp = 1, isorpmax_xc
          hx_den3(isorp,index) = (x3cmax_den3(isorp,index) - xmin)/real(numx3c_den3(isorp,index) - 1)
          hy_den3(isorp,index) = (y3cmax_den3(isorp,index) - ymin)/real(numy3c_den3(isorp,index) - 1)
        end do
      end do
    end do
  end do
end subroutine setterp_2d
