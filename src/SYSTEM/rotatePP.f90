! rotatePP.f90
! This routine rotates a matrix from molecular to crystal coordinates.
!
! The variable eps is a 3x3 output of the subroutine epsilon.
! The variable dmat is a 5x5 matrix rotating d-orbitals.
! The variable pmat is a 3x3 matrix rotating p-orbitals.
!
! Here is the famous Ortega convention:
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                       1    2    3
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                  1     2   3     4      5
!
! ===========================================================================
subroutine rotatePP(in1,in2,eps,mmatrix,xmatrix)
 use M_system
 use M_fdata, only: nssh,lssh,nsshPP,lsshPP
 implicit none
 integer, intent(in) :: in1
 integer, intent(in) :: in2
 real, intent(in) :: eps (3, 3)
 real, intent(in) :: mmatrix (numorb_max, numorb_max)
 real, intent(out) :: xmatrix (numorb_max, numorb_max)
 
 integer issh
 integer jssh
 integer k1, k2
 integer n1, l1, m1
 integer n2, l2, m2
   real dmat (5, 5)
 real left (5, 5)
 real pmat (3, 3)
 real right (5, 5)
 
 ! IMPORTANT: for the pseudopotential-matrix-elements, the first orbital is a valence orbital and the second is a "pseudopotential" orbital
 call twister (eps, dmat, pmat)
 xmatrix=0.0d0
 n1 = 0
 do issh = 1, nssh(in1)
  l1 = lssh(issh,in1)
  call chooser (l1, dmat, pmat, left)
  n2 = 0
  do jssh = 1, nsshPP(in2)
   l2 = lsshPP(jssh,in2)
   call chooser (l2, dmat, pmat, right)
   do m2 = 1, 2*l2 + 1
    do m1 = 1, 2*l1 + 1
     do k2 = 1, 2*l2 + 1
      do k1 = 1, 2*l1 + 1
       xmatrix(n1+k1,n2+k2) = xmatrix(n1+k1,n2+k2) + left(k1,m1)*mmatrix(n1+m1,n2+m2)*right(k2,m2)
      end do
     end do
    end do
   end do
   n2 = n2 + 2*l2 + 1
  end do
  n1 = n1 + 2*l1 + 1
 end do
 return
end
 
