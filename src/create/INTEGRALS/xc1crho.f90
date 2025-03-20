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
! onecenternuxc.f90
! Program Description
! ===========================================================================
!       This routine calculates the one-center integrals for the exchange-
! correlation interactions of the extended Hubbard model.
!
!   int [n(ishell)n(jshell) Nuxc(natom)]
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine xc1crho (nspec, nspec_max, nsh_max, wfmax_points,   &
     &                            iexc, fraction, nsshxc, rcutoffa_max,      &
     &                            xnocc, dqorb, iderorb, what, signature,    &
     &                            drr_rho)
        use constants
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iexc
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max
        integer, intent (in) :: wfmax_points

        integer, intent (in), dimension (nspec_max) :: iderorb
        integer, intent (in), dimension (nspec_max) :: nsshxc
 
        real(kind=long), intent (in) :: fraction
 
        real(kind=long), intent (in), dimension (nspec_max) :: dqorb
        real(kind=long), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=long), intent (in), dimension (nspec_max) :: rcutoffa_max
        real(kind=long), intent (in), dimension (nsh_max, nspec_max) :: xnocc
 
        character (len=70), intent (in) :: signature

        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Output
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ideriv
        integer in1
        integer in2
        integer irho
        integer issh
        integer jssh
        integer lssh
        integer nnrho
        integer nssh
 
        real(kind=long) dnuxc
        real(kind=long) dnuxcs
        real(kind=long) dq
        real(kind=long) drho
        real(kind=long) exc
        real(kind=long) factor
        real(kind=long) rcutoff
        real(kind=long) rho
        real(kind=long) rhomin
        real(kind=long) rhomax
        real(kind=long) rh
        real(kind=long) rhp
        real(kind=long) rhpp
        real(kind=long) vxc
 
        real(kind=long), dimension (:, :), allocatable :: answer
        real(kind=long), dimension (:), allocatable :: rho1c
        real(kind=long), dimension (:), allocatable :: rhop1c
        real(kind=long), dimension (:), allocatable :: rhopp1c
        real(kind=long), dimension (:), allocatable :: xnocc_in
 
        real(kind=long), external :: psiofr

	character(2)  :: shell
        character(80) :: fname
         
! Procedure
! ===========================================================================
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/nuxc_onecenter.dat',                &
     &        status = 'unknown')
 
! Set up the header for the output file.
        write (36,100)
        write (36,*) ' All one center matrix elements '
        write (36,*) ' created by: '
        write (36,200) signature
 
        do in1 = 1, nspec
         write (36,300) what(in1)
        end do
        write (36,100)
 
! Loop over the different charge types (0, -1, or +1).
        ideriv = 0
 
! Loop over the species
        allocate (rho1c (wfmax_points))
        allocate (rhop1c (wfmax_points))
        allocate (rhopp1c (wfmax_points))
        do in1 = 1, nspec
         nssh = nsshxc(in1)
         write (36,400) in1, nssh
 
! Needed for charge corrections:
         dq = dqorb(in1)
         jssh = iderorb(in1)
 
         drho = drr_rho(in1)
         rcutoff = rcutoffa_max(in1)
         allocate (xnocc_in (nssh))
         xnocc_in(1:nssh) = xnocc(1:nssh,in1)
 
! Obtain the density and respective derivatives needed for evaluating the
! exchange-correlation interactions (LDA or GGA).
         call rho1c_store (in1, nsh_max, nssh, dq, jssh, drho, rcutoff,      &
     &                     xnocc_in, ideriv + 1, wfmax_points, rho1c, rhop1c,&
     &                     rhopp1c)

! Loop over shell (density)
         do lssh = 1,nsshxc(in1)

          write (*,*) 'lssh = ',lssh

! create filename
          write (shell,'(i2.2)') lssh
          fname = 'coutput/exc1crho.'//shell//'.dat'

! Open the file to store the onecenter data.
          write (*,*) ' open file ',fname
          open (unit = 36, file = fname , status = 'unknown')

! Set up the header for the output file.
          write (36,100)
          write (36,*) ' All one center matrix elements '
          write (36,*) ' created by: '
          write (36,200) signature

          do in2 = 1, nspec
           write (36,300) what(in2)
          end do
          write (36,100)
          write (36,400)  in1, nssh
 
! Integrals <i|exc(i)-mu(i)|i> and <i.nu|mu(i)|i.nu'>
! ***************************************************************************
! First initialize the answer array
          allocate (answer (nssh, nssh))
          answer = 0.0d0
 
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
 
! Compute the exchange correlation potential
           rh = rho1c(irho)*abohr**3
           rhp = rhop1c(irho)*abohr**4
           rhpp = rhopp1c(irho)*abohr**5
           call get_potxc1c (iexc, fraction, rho, rh, rhp, rhpp, exc, vxc,    &
     &                      dnuxc, dnuxcs)
 
! Convert to eV
           dnuxc = dnuxc*Hartree*abohr**3
 
           do issh = 1, nssh
            do jssh = 1, nssh
             answer(issh,jssh) = answer(issh,jssh)                            &
     &        + dnuxc*factor*rho**2*(psiofr(in1,lssh,rho)**2)                 &
     &         *(psiofr(in1,issh,rho)*psiofr(in1,jssh,rho)/(4.0d0*pi))
            end do
           end do
          end do
          do issh = 1, nssh
           write (36,500) answer(issh,1:nssh)
          end do
          deallocate (answer)
         end do
        end do
        deallocate (xnocc_in)
        write (36,*) '  '
        write (*,*) '  '
        write (*,*) ' Writing output to:  coutput/exc1crho.XX.dat '
        write (*,*) '  '
 
        close (unit = 36)

! Deallocate Arrays
! ===========================================================================
        deallocate (rho1c, rhop1c, rhopp1c)
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (2x, i3, 2x, i3)
500     format (8d20.10)
        return
        end