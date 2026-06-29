! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
 
!
! fireball-qmd is a free (GPLv3) open project.

!
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

 
! x_1c.f
! Program Description
! ===========================================================================
!
!       This routine calculates the one-center integrals for the exact
! exchange interactions.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine x_1c (nsh_max, nspec, nspec_max, fraction, nsshxc,       &
     &                   lsshxc, drr_rho, rcutoffa_max, what, signature)
        use iso_fortran_env, only: dp => real64
        use :: wavefunctions, only: wf_atoms
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max
 
        integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
        integer, intent (in), dimension (nspec_max) :: nsshxc
 
        real(kind=dp), intent (in) :: fraction
 
        real(kind=dp), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=dp), intent (in), dimension (nspec_max) :: rcutoffa_max
 
        character (len=70), intent (in) :: signature
 
        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: lmax = 3
 
        real(kind=dp), parameter :: eq2 = 14.39975_dp
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer issh
        integer itype
        integer jssh
        integer lalpha
        integer lmu
        integer lqn
        integer malpha
        integer mmu
        integer mqn
        integer nnrho
 
        real(kind=dp) coefficient
        real(kind=dp) cg1
        real(kind=dp) cg2
        real(kind=dp) cg3
        real(kind=dp) cg4
        real(kind=dp) drho
        real(kind=dp) psi1
        real(kind=dp) psi2
        real(kind=dp) psi3
        real(kind=dp) psi4
        real(kind=dp) r
        real(kind=dp) rp
        real(kind=dp) rhomax
        real(kind=dp) rhomin
        real(kind=dp) sumr
        real(kind=dp) sumrp
 
        real(kind=dp), dimension (-lmax:lmax) :: answer
        real(kind=dp), dimension (:), allocatable :: factor
        real(kind=dp), dimension (:), allocatable :: rpoint
 
        real(kind=dp), external :: clebsch_gordon
 
! Procedure
! ===========================================================================
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/exchange_1c.dat', status = 'unknown')
 
! Set up the header for the output file.
        write (36,100)
        write (36,*) ' All one center matrix elements created by: '
        write (36,200) signature
 
        do itype = 1, nspec
         write (36,300) what(itype)
        end do
        write (36,100)
 
! Loop over each of the species, then loop over each shell of each species.
! For each shell calculate the exact exchange energy contribution. The total
! exchange-energy contribution will be the sum of each contribution of the
! shell multiplied by the net occupation number of that particular shell.
        do itype = 1, nspec
 
! Initialize the limits of integration for the radial integral.
! Set up the grid points for the rho integration.
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
         end do
 
! Loop over the orbitals. Only the mu = nu matrix elements survive in the
! one-center exchange interactions.
         do issh = 1, nsshxc(itype)
          lmu = lsshxc(itype,issh)
          do mmu = -lmu, lmu
 
! Loop over the orbitals again, this loop will provide the sum over all
! the states-alpha.  Each piece in the alpha sum is multiplied by the
! occupation number. The sum will be performed in the Fireball code.
! Just calculate each contribution, and write out the answer for each
! alpha piece.
           do jssh = 1, nsshxc(itype)
            lalpha = lsshxc(itype,jssh)
            do malpha = -lalpha, lalpha
             answer(malpha) = 0.0_dp
 
! Next perform a sum over all quantum numbers l (up to lmax) and all
! corresponding quantum numbers m.
             do lqn = 0, 2*lmax
              do mqn = -lqn, lqn
 
! Calculate the Clebsch-Gordon coefficient which results from the theta', phi'
! integration. If this Clebsch-Gordon coefficient is zero, then skip this
! lqn, and mqn in the sum.
               cg1 = clebsch_gordon (lalpha, 0, lqn, 0, lmu, 0)
               cg2 = clebsch_gordon (lalpha, malpha, lqn, mqn, lmu, mmu)
               cg3 = clebsch_gordon (lqn, 0, lalpha, 0, lmu, 0)
               cg4 = clebsch_gordon (lqn, mqn, lalpha, malpha, lmu, mmu)
               coefficient = cg1*cg2*cg3*cg4*(2.0_dp*real(lalpha) + 1)       &
     &                                      /(2.0_dp*real(lmu) + 1)
 
! ****************************************************************************
! Perform the radial integration. Only do this integral if the coefficient is
! non-zero.
               if (coefficient .gt. 1.0e-4_dp) then
 
! First integrate the even pieces and then the odd pieces
                sumr = 0.0_dp
                do irho = 1, nnrho
                 r = rpoint(irho)
                 if (r .lt. 1.0e-4_dp) r = 1.0e-4_dp
                 psi1 = wf_atoms(itype)%get_psi(issh,r)
                 psi2 = wf_atoms(itype)%get_psi(jssh,r)
 
! ****************************************************************************
! Perform the radial integration over r'.
! Limits from 0 to r.
                 sumrp = 0.0_dp
                 do irhop = 1, nnrho
                  rp = rpoint(irhop)
                  if (rp .lt. 1.0e-4_dp) rp = 1.0e-4_dp
                  psi3 = wf_atoms(itype)%get_psi(jssh,rp)
                  psi4 = wf_atoms(itype)%get_psi(issh,rp)
                  if (rp .le. r) then
                   sumrp =                                                  &
     &              sumrp + factor(irhop)*psi3*psi4*rp**(lqn + 2)/r**(lqn + 1)
 
! Limits from 0 to rcutoff.
                  else
                   sumrp = sumrp + factor(irhop)*psi3*psi4*r**lqn/rp**(lqn - 1)
                  end if
                 end do
! ****************************************************************************
 
                 sumr = sumr + coefficient*factor(irho)*sumrp*psi1*psi2*r**2
                end do
               end if
! ****************************************************************************
 
! End loop over lqn and mqn
              end do
             end do
             answer(malpha) = answer(malpha) + (eq2/2.0_dp)*fraction*sumr
 
! End the loop over malpha
            end do
 
! Write out the contribution from the alpha orbital at this point.
            write (36,400) (answer(malpha), malpha = -lalpha, lalpha)
 
! End the loop over jssh
           end do
 
! End loop over issh, mmu.
          end do
         end do
 
! End loop over the species.
         deallocate (rpoint)
         deallocate (factor)
        end do
 
        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/exchange_1c.dat '
        write (*,*) '  '
 
        close (unit = 36)
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (8d20.10)
 
        return
        end
