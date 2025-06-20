subroutine lowdin_charges()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: spin
  use M_system, only: igamma, natoms, degelec, imass, blowre, blowim, norbitals, nkpoints, ioccupy_k, foccupy, weight_k, Qout, &
    & QLowdin_TOT
  use M_fdata, only: nssh,lssh
  implicit none
  integer iatom
  integer ikpoint
  integer imu
  integer in1
  integer issh, mmu
  integer mqn
  integer iorbital
  real(double) aux1, aux2, aux3
  Qout = 0.0d0
  QLowdin_TOT = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    do ikpoint = 1, nkpoints
       aux1 = weight_k(ikpoint)*spin
       do iorbital = 1, norbitals
         if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
          aux2 = aux1*foccupy(iorbital,ikpoint)
          imu = 0
          do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
              imu = imu + 1
              mmu = imu + degelec(iatom)
              if (igamma .eq. 0) then
                aux3 = aux2*(blowre(mmu,iorbital,ikpoint)**2  + blowim(mmu,iorbital,ikpoint)**2)
              else
                aux3 = aux2*blowre(mmu,iorbital,ikpoint)**2
              end if
              Qout(issh,iatom) = Qout(issh,iatom) + aux3
              QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3
            end do
          end do
        end if
      end do !kpoints
    end do
  end do !atoms
end subroutine lowdin_charges
