! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
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

!
! goverlap1c.f90
! Program Description
! ===========================================================================
!       This routine calculates the one-center integrals for the
! calculation of the  matrix elements
!  < mu | Grad nu >
! where Grad is the gradient wrt the position of the atom where mu and nu are
!
! we calculate two different integrals:
! 
!  F1 = Int fmu(r) d/dr (fnu(r)) r**2 dr 
!  F2 = Int fmu(r) fnu(r) r dr 
!
!fmu(r) and fnu(r) are the radial parts of the mu , nu orbitals
!
! ===========================================================================
! Code written from nuxc1crho.f90 by
! Jose Ortega Mateo
! Dept. Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine goverlap1c (nspec, nspec_max, nsh_max, wfmax_points, &
     &                            nsshxc, lsshxc, rcutoffa_max, what,   &
     &                            signature, drr_rho)
        use constants
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max
        integer, intent (in) :: wfmax_points

        integer, intent (in), dimension (nspec_max) :: nsshxc
	integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
 
 
        real(kind=long), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=long), intent (in), dimension (nspec_max) :: rcutoffa_max
 
        character (len=70), intent (in) :: signature

        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Output
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer in2
        integer irho
        integer issh
        integer jssh
        integer lssh
        integer nnrho
        integer nssh
        integer l1
        integer l2
 
        real(kind=long) dnuxc
        real(kind=long) dnuxcs
        real(kind=long) dq
        real(kind=long) drho
        real(kind=long) exc
        real(kind=long) dexc
        real(kind=long) factor
        real(kind=long) rcutoff
        real(kind=long) rho
        real(kind=long) rhomin
        real(kind=long) rhomax
        real(kind=long) rh
        real(kind=long) rhp
        real(kind=long) rhpp
        real(kind=long) vxc
 
        real(kind=long), dimension (:, :), allocatable :: answer1
        real(kind=long), dimension (:, :), allocatable :: answer2
        real(kind=long), dimension (:), allocatable :: rho1c
        real(kind=long), dimension (:), allocatable :: rhop1c
        real(kind=long), dimension (:), allocatable :: rhopp1c
        real(kind=long), dimension (:), allocatable :: xnocc_in
 
        real(kind=long), external :: psiofr
        real(kind=long), external :: dpsiofr

        character(80) :: fname
         
! Procedure
! ===========================================================================
! Create file name for F1

        fname = 'coutput/goverlapf1.dat'

! Open the file to store the onecenter data.
        write (*,*) ' open file ',fname
        open (unit = 36, file = fname , status = 'unknown')

! Set up the header for the output file.
        write (36,100)
        write (36,*) 'One center matrix elements for goverlap1c: F1 '
        write (36,*) ' created by: '
        write (36,200) signature

        do in2 = 1, nspec
         write (36,300) what(in2)
        end do
        write (36,100)
!
! Create file name for F2

        fname = 'coutput/goverlapf2.dat'

! Open the file to store the onecenter data.
        write (*,*) ' open file ',fname
        open (unit = 37, file = fname , status = 'unknown')

! Set up the header for the output file.
        write (37,100)
        write (37,*) 'One center matrix elements for goverlap1c: F2 '
        write (37,*) ' created by: '
        write (37,200) signature

        do in2 = 1, nspec
         write (37,300) what(in2)
        end do
        write (37,100)

 
 
! Loop over the species
        do in1 = 1, nspec

         nssh = nsshxc(in1)

 
         drho = drr_rho(in1)
         rcutoff = rcutoffa_max(in1)
 
! Loop over shell (density)
!         do lssh = 1,nsshxc(in1)
!
!          write (36,400)  in1, nssh, lssh
           write (36,800)  in1, nssh
           write (37,800)  in1, nssh
 
! ***************************************************************************
! First initialize the answer array
          allocate (answer1 (nssh, nssh))
          allocate (answer2 (nssh, nssh))
          answer1 = 0.0d0
          answer2 = 0.0d0
 
! Fix the endpoints and initialize the increments dz and drho.
          rhomin = 0.0d0
          rhomax = rcutoff
 
          nnrho = nint((rhomax - rhomin)/drho) + 1
 
! Here we loop over rho.
          do irho = 1, nnrho
           rho = rhomin + real(irho - 1, kind=long)*drho
 
           factor = 2.0d0*drho/3.0d0
           if (mod(irho, 2) .eq. 0) factor = 4.0d0*drho/3.0d0
           if (irho .eq. 1 .or. irho .eq. nnrho) factor = drho/3.0d0
 
 
! JOM-test
            write(*,*)'dpsiofr', rho , psiofr(1,2,rho), dpsiofr(1,2,rho)
           do issh = 1, nssh
            l1 = lsshxc(in1,issh)
            do jssh = 1, nssh
             l2 = lsshxc(in1,jssh)
              answer1(issh,jssh) = answer1(issh,jssh)                            &
     &         + factor*rho**2*psiofr(in1,issh,rho)*dpsiofr(in1,jssh,rho)  
              answer2(issh,jssh) = answer2(issh,jssh)                            &
     &         + factor*rho*psiofr(in1,issh,rho)*psiofr(in1,jssh,rho)  
            end do
           end do
          end do
          do issh = 1, nssh
           write (36,500) answer1(issh,1:nssh)
           write (37,500) answer2(issh,1:nssh)
          end do ! do issh
          deallocate (answer1)
          deallocate (answer2)
!         end do ! lssh 
        end do ! do in1 = 1, ispec
        write (36,*) '  '
        write (37,*) '  '
        write (*,*) '  '
        write (*,*) ' Writing output to:  coutput/goverlapf1.dat '
        write (*,*) ' Writing output to:  coutput/goverlapf2.dat '
        write (*,*) '  '
 
        close (unit = 36)
        close (unit = 37)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (2x, i3, 2x, i3, 2x, i3)
800     format (2x, i3, 2x, i3)
500     format (8d20.10)
600     format ('dExc=',5f12.6)
        return
        end