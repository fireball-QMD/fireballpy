! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! gaunt.f90
! Program Description
! ===========================================================================
!       This routine calculates a radial integral needed in onecentervdip 
! ===========================================================================
! Code originally written by:
! Code written by:
! Daniel G. Trabada
! Departamento de físca teórica de la materia condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================

! Program Declaration
! ===========================================================================
        real(kind=dp) function radial_integral(itype, lqn, il1, il2, il3, il4, nspec_max, drr_rho, rcutoffa_max)
        use iso_fortran_env, only: dp => real64
        use :: wavefunctions, only: wf_atoms
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: il1
        integer, intent (in) :: il2
        integer, intent (in) :: il3
        integer, intent (in) :: il4
        integer, intent (in) :: lqn
        integer, intent (in) :: itype
        integer, intent (in) :: nspec_max
        real(kind=dp), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=dp), intent (in), dimension (nspec_max) :: rcutoffa_max

        integer :: irho, irhop
        integer nnrho
        real(kind=dp) drho
        real(kind=dp) rhomax
        real(kind=dp) rhomin

        real(kind=dp) :: sumr, r, sumrp, rp, psi1, psi2, psi3, psi4


        real(kind=dp), dimension (:), allocatable :: factor
        real(kind=dp), dimension (:), allocatable :: rpoint

!-------------------------------------------------------------------------------------

              rhomin = 0.0_dp
         rhomax = rcutoffa_max(itype)
         drho = drr_rho(itype)
         nnrho = int((rhomax - rhomin)/drho) + 1
         allocate (rpoint(nnrho))
         allocate (factor(nnrho))
          
         
        
         do irho = 1, nnrho
          rpoint(irho) = float(irho - 1)*drho
! Set up the Simpson rule factors:
          factor(irho) = 2.0_dp*drho/3.0_dp
          if (mod(irho,2) .eq. 0) factor(irho) = 4.0_dp*drho/3.0_dp
          if (irho .eq. 1 .or. irho .eq. nnrho) factor(irho) = drho/3.0_dp
         end do !irho
!-------------------------------------------------------------------------------------
       !START RADIAL INTEGRAL
       !-----------------------------------------------------------
! First integrate the even pieces and then the odd pieces
              sumr = 0.0_dp
              do irho = 1, nnrho
               r = rpoint(irho)
               if (r .lt. 1.0e-4_dp) r = 1.0e-4_dp
                psi1 = wf_atoms(itype)%get_psi(il1,r)
                psi2 = wf_atoms(itype)%get_psi(il2,r)
 
! ****************************************************************************
! Perform the radial integration over r'.
! Limits from 0 to r.
                sumrp = 0.0_dp
                do irhop = 1, nnrho
                 rp = rpoint(irhop)
                 if (rp .lt. 1.0e-4_dp) rp = 1.0e-4_dp
                 psi3 = wf_atoms(itype)%get_psi(il3,rp)
                 psi4 = wf_atoms(itype)%get_psi(il4,rp)
                 if (rp .le. r) then
                  sumrp = sumrp + factor(irhop)*psi3*psi4*rp**(lqn + 2)/r**(lqn + 1)   ! Limits from 0 to rcutoff.
                 else
                  sumrp = sumrp + factor(irhop)*psi3*psi4*r**lqn/rp**(lqn - 1)
                 end if
!                if (sumr .gt. 1.0e-4_dp ) then
!                    write(36,'(10F14.10)')  sumr, psi1,psi2,r,sumrp,psi3,psi3,rp
!                end if
                end do !rhop
! ****************************************************************************
                sumr = sumr + factor(irho)*sumrp*psi1*psi2*r**2  !*coefficient
               end do !irho
         !------------------------------------------------------------
         !END OF START RADIAL INTEGRAL
         !-----------------------------------------------------------
         !deallocate (rpoint)
         !deallocate (factor)
         radial_integral = sumr
        return
        end function
