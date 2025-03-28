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


! onecenterxc.f90
! Program Description
! ===========================================================================
!       This routine calculates the one-center integrals for the exchange-
! correlation interactions.
!
!   int [n(i) (Exc - Muxc)] and <i.nu | mu[n(i)] | i.nu'>
!
!       The full charge transfer OLSXC version calculates exc and 
!  < imu|muxc|inu> separately. 
! ===========================================================================
! Original code from Juergen Fritsch
! 
! Code rewritten by:
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
        subroutine onecenterxc (nspec, nspec_max, nsh_max, wfmax_points,     &
     &                          iexc, fraction, nsshxc, lsshxc, rcutoffa_max,&
     &                          xnocc, dqorb, iderorb, what, signature,      &
     &                          drr_rho ,nzx )
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
        integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
 
        real(kind=long), intent (in) :: fraction
 
        real(kind=long), intent (in), dimension (nspec_max) :: dqorb
        real(kind=long), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=long), intent (in), dimension (nspec_max) :: rcutoffa_max
        real(kind=long), intent (in), dimension (nsh_max, nspec_max) :: xnocc
        integer, intent (in), dimension (nspec_max) :: nzx
        
        real(kind=long) tmp

        character (len=70) :: signature

        character (len=70), intent (in), dimension (nspec_max) :: what
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ideriv
        integer irho
        integer issh
        integer iissh
        integer in1
        integer jssh
        integer jjssh
        integer L1
        integer L2
        integer ndq
        integer nnrho
        integer nssh
 
        real(kind=long) dnuxc
        real(kind=long) dnuxcs
        real(kind=long) dq
        real(kind=long) drho
        real(kind=long) exc
        real(kind=long) dexcc
        real(kind=long) exc1c_0
        real(kind=long) factor
        real(kind=long) rcutoff
        real(kind=long) rho
        real(kind=long) rhomax
        real(kind=long) rhomin
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

        real(kind=long) qa,qb
        real(kind=long) dq1,dq2,f0,f1,f2
        real(kind=long) denom
        real(kind=long) dqi(nsh_max,3)
        real(kind=long) qmax(nsh_max)
        real(kind=long) eexc(nsh_max,nsh_max)
        real(kind=long) dexc(nsh_max)
        real(kind=long) d2exc(nsh_max,nsh_max)
        real(kind=long) vvxc(nsh_max,nsh_max)
        real(kind=long) dvxc(nsh_max,nsh_max,nsh_max)
        real(kind=long) d2vxc(nsh_max,nsh_max,nsh_max,nsh_max)
        integer imask (nsh_max+1)

        character(80) ::  root
        character(2) :: auxz
 
! Procedure
! ===========================================================================
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/xc_onecenter.dat', status = 'unknown')
 
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
        allocate (rho1c (wfmax_points))
        allocate (rhop1c (wfmax_points))
        allocate (rhopp1c (wfmax_points))
        do ideriv = 0, 2
         write (36,*) ' derivative ', ideriv
 
! Loop over the species
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
          call rho1c_store (in1, nsh_max, nssh, dq, jssh, drho, rcutoff,     &
     &                      xnocc_in, ideriv + 1, wfmax_points, rho1c,       &
     &                      rhop1c, rhopp1c)
 
! Integrals <i|exc(i)-mu(i)|i> and <i.nu|mu(i)|i.nu'>
! ***************************************************************************
! First initialize the answer array
          exc1c_0 = 0.0d0
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
! Convert to atomic units.
           rho = rho/abohr
           rh = rho1c(irho)*abohr**3
           rhp = rhop1c(irho)*abohr**4
           rhpp = rhopp1c(irho)*abohr**5
           call get_potxc1c (iexc, fraction, rho, rh, rhp, rhpp, exc, vxc,   &
     &                       dnuxc, dnuxcs, dexcc)
 
! Convert to eV units
           vxc = hartree*vxc
           exc = hartree*exc
           rho = rho*abohr
 
! Add to integral -- factor*rho*rho weight factor for the radial integral
           exc1c_0 = exc1c_0 + 4.0d0*pi*rho1c(irho)*(exc - vxc)*factor*rho**2
           do issh = 1, nssh
            L1 = lsshxc(in1,issh)
            do jssh = 1, nssh
             L2 = lsshxc(in1,jssh)
             if (L1 .eq. L2) then
              answer(issh,jssh) = answer(issh,jssh)                          &
     &         + psiofr(in1,issh,rho)*vxc*psiofr(in1,jssh,rho)*factor*rho**2
             end if
            end do
           end do
          end do
          write (36,500) exc1c_0
          do issh = 1, nssh
           write (36,501) answer(issh,1:nssh)
          end do
          deallocate (xnocc_in)
          deallocate (answer)
         end do
        end do
 
        write (36,*) '  '
        do in1 = 1, nspec
         write (36,600) iderorb(in1), dqorb(in1)
        end do
 
        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/xc_onecenter.dat '
        write (*,*) '  '
 
        close (unit = 36)

! ***************************************************************************
! OSLXC - charge transfer, <imu|V_xc(rho_0+rho')|inu>
! We only use charge correction on XC-energy now (2-center Vxc charge 
! correction must be included). The charge correction is evaluated but is not 
! used in FIREBALL. The information is read in FIREBALL, but not used.
! ***************************************************************************
! Open the file to store the onecenter data.
        open (unit = 36, file = 'coutput/xc1c_dqi.dat', status = 'unknown')
 
! Set up the header for the output file.
        write (36,100)
        write (36,*) ' All one center matrix elements '
        write (36,*) ' created by: '
        write (36,200) signature
 
        do in1 = 1, nspec
         write (36,300) what(in1)

         write (auxz,'(i2.2)') nzx(in1)
         root = 'coutput/xc1c_dqi.'//auxz//'.dat'
         open (unit = 360, file = root , status = 'unknown')
         write (360,100)
         write (360,*) ' All one center matrix elements '
         write (360,*) ' created by: '
         write (360,200) signature
         write (360,300) what(in1)
         write (360,100)
         close(360)

        end do
        write (36,100)

! set number of derivatives steps to take
        ndq = 3 

! Loop over the species
        do in1 = 1, nspec
         nssh = nsshxc(in1)

         write (36,400) in1, nssh
         write (auxz,'(i2.2)') nzx(in1)
         root = 'coutput/xc1c_dqi.'//auxz//'.dat'
         open (unit = 360, file = root , position='append', status = 'old')
         write (360,400)  in1, nssh

         allocate (xnocc_in (nssh))

! Needed for charge corrections:
         drho = drr_rho(in1)
         rcutoff = rcutoffa_max(in1)
! Set charges with respect to dqi
         do jssh = 1, nssh
          xnocc_in(jssh) = xnocc(jssh,in1)
         end do

! Obtain the density and respective derivatives needed for evaluating the
! exchange-correlation interactions (LDA or GGA).
! We have to avoid change xnocc_in !!
         call rho1c_store (in1, nsh_max, nssh, 0.0d0, 1, drho, rcutoff,    &
     &                       xnocc_in, 1, wfmax_points, rho1c, rhop1c, rhopp1c)
 
! Integrals <i|exc(i)|i> and <i.nu|mu(i)|i.nu'>
! ***************************************************************************
! First initialize the answer array
         eexc = 0.0d0
         vvxc = 0.0d0
 
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
! Convert to atomic units.
          rho = rho/abohr
          rh = rho1c(irho)*abohr**3
          rhp = rhop1c(irho)*abohr**4
          rhpp = rhopp1c(irho)*abohr**5
          call get_potxc1c (iexc, fraction, rho, rh, rhp, rhpp, exc, vxc,  &
   &                        dnuxc, dnuxcs, dexcc)
 
! Convert to eV units
          vxc = hartree*vxc
          exc = hartree*exc
          rho = rho*abohr
 
! Add to integral -- factor*rho*rho weight factor for the radial integral
          do issh = 1, nssh
           L1 = lsshxc(in1,issh)
            do jssh = 1, nssh
             L2 = lsshxc(in1,jssh)
             if (L1 .eq. L2) then
              vvxc(issh,jssh) =  vvxc(issh,jssh)                           &
     &          + psiofr(in1,issh,rho)*vxc*psiofr(in1,jssh,rho)*factor*rho**2
              eexc(issh,jssh) =  eexc(issh,jssh)                           &
     &          + psiofr(in1,issh,rho)*exc*psiofr(in1,jssh,rho)*factor*rho**2
             end if
            end do
          end do
         end do ! do irho = 1, nnrho

! Print exc terms
         do issh = 1, nssh
          write (36,501) (eexc(issh,jssh),jssh = 1, nssh)
          write (360,501) (eexc(issh,jssh),jssh = 1, nssh)
         end do

         write (36,*)
         write (360,*)
! Print vxc terms
         do issh = 1, nssh
          write (36,501) (vvxc(issh,jssh),jssh = 1, nssh)
          write (360,501) (vvxc(issh,jssh),jssh = 1, nssh)
         end do

         deallocate (xnocc_in)
        
        close(360)

        end do ! do in1 = 1, nspec

        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/xc1c_dqi.dat '
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
450     format (2x, i3, 2x, i3, 2x, i3)
470     format (2x, f7.3, 1x, f7.3, 1x, f7.3)
480     format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
500     format (8d20.10)
501     format (8d20.10)
600     format (1x, i3, 2x, f10.5)
 
        return
      end subroutine onecenterxc
