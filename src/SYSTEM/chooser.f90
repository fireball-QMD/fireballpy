        subroutine chooser (l, dmat, pmat, rmatrix)
        implicit none
        integer, intent(in) :: l
        real, intent(in), dimension(5, 5) :: dmat
        real, intent(in), dimension(3, 3) :: pmat
        real, intent(out), dimension(5, 5) :: rmatrix
        rmatrix = 0.0d0
        if (l .eq. 0) then
         rmatrix(1,1) = 1.0d0
        else if (l .eq. 1) then
         rmatrix(1:3,1:3) = pmat
        else if (l .eq. 2) then
         rmatrix = dmat
        end if
        return
        end

