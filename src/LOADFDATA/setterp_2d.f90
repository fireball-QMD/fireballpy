subroutine setterp_2d ()
  use M_fdata    
  implicit none
  integer in1
  integer in2
  integer in3
  integer index
  integer isorp
  real(8) xmin
  real(8) ymin
  xmin = 0.0d0
  ymin = 0.0d0
  do in1 = 1, nspecies
    do in2 = 1, nspecies
      do in3 = 1, nspecies
        index = icon3c(in1,in2,in3)

        do isorp = 0, isorpmax
          hx_bcna(isorp,index) = (x3cmax_bcna(isorp,index) - xmin)/real(numx3c_bcna(isorp,index) - 1, 8)
          hy_bcna(isorp,index) = (y3cmax_bcna(isorp,index) - ymin)/real(numy3c_bcna(isorp,index) - 1, 8)
        end do

        do isorp = 1, isorpmax_xc
          hx_den3(isorp,index) = (x3cmax_den3(isorp,index) - xmin)/real(numx3c_den3(isorp,index) - 1, 8)
          hy_den3(isorp,index) = (y3cmax_den3(isorp,index) - ymin)/real(numy3c_den3(isorp,index) - 1, 8)
        end do

      end do   
    end do   
  end do  
end subroutine setterp_2d
