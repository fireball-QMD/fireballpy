! copyright info:
!
!                             @Copyright 1998
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
!                      Richard B. Evans
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

!
! fireball-qmd is a free (GPLv3) open project.

!
! Use, duplication, or disclosure of this software and its
! documentation by the Government is subject to restrictions
! as set forth in subdivision { (b) (3) (ii) } of the Rights
! in Technical Data and Computer Software clause at 52.227-7013.
!
! rho1c_store.f
! Program Description
! ====================================================================
!
! This routine calculates and stores the combined density of species
! 1 and 2 as a function of r, z and d.
!
! On input:
!            in1:     The species type of atom 1
!
!     By common blocks:
!            nsshxc:     Number of shells for each species
!
! On output: rho1c:      The density as a function of species type,
!                        r, z and d.  Output is placed in common block
!                        density located in wavefunctions.inc
!            rhop1c:     derivative with respect to rho
!            rhopp1c:    second derivative with respect to rho
!
! ====================================================================
! Code written by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ====================================================================
!
! Program Declaration
! ====================================================================
        subroutine rho1c_store (in1, nsh_max, nssh, dq, jssh, drho,
     1                          rcutoff, xnocc, ix, wfmax_points,
     2                          rho1c, rhop1c, rhopp1c)
        implicit none

! Input Declaration and Description
! ====================================================================
! Input
        integer ix
        integer in1
        integer jssh
        integer nsh_max
        integer nssh
        integer wfmax_points

        real*8 dq
        real*8 drho
        real*8 rcutoff

        real*8 xnocc (nsh_max)

! Output
        real*8 rho1c (wfmax_points)
        real*8 rhop1c (wfmax_points)
        real*8 rhopp1c (wfmax_points)

! Local Parameters and Data Declaration
! ====================================================================
        integer j0 (3)
        data j0 /0, -1, 1/

! Local Variable Declaration and Description
! ====================================================================
        integer irho
        integer issh
        integer nnrho

        real*8 dens
        real*8 psiofr
        real*8 rho
        real*8 rhomax
        real*8 rhomin
        real*8 xinv4pi

        external psiofr

! Procedure
! ====================================================================
! Initialize 1/4pi
        xinv4pi = 1.0d0/(4.0d0*3.141592653589793238462643d0)

! Fix the endpoints and initialize the increments dz and drho.
        rhomin = 0.0d0
        rhomax = rcutoff

        nnrho = nint((rhomax - rhomin)/drho) + 1

! Performs some checks to make sure that the dimensions on rho1c are
! alright.
        if (nnrho .gt. wfmax_points) then
         write (*,*) ' In rho1c_store.f, nnrho = ', nnrho
         write (*,*) ' The dimension nrho_points =', wfmax_points
         write (*,*) ' needs to be increased. '
         stop 'error in rho1c_store'
        end if

! Here we loop over rho computing the sum of the densities for species
! in1 at each value of rho.
        do irho = 1, nnrho
         rho = rhomin + dfloat(irho - 1)*drho

         dens = 0.0d0
         do issh = 1, nssh
          dens = dens + xnocc(issh)*psiofr(in1,issh,rho)**2
         end do

! Here the derivative with respect to the charge correction term
! is calculated.  The variable switch determines whether the correction is
! for the one, two or three center case.
         dens = dens + j0(ix)*dq*psiofr(in1,jssh,rho)**2
         rho1c(irho) = dens*xinv4pi
        end do

! Now calculate the derivatives - rhop1c and rhopp1c
        do irho = 2, nnrho - 1
         rho = rhomin + dfloat(irho - 1)*drho

! First derivatives:
         rhop1c(irho) = (rho1c(irho+1) - rho1c(irho-1))/(2.0d0*drho)

! Second derivatives:
         rhopp1c(irho) =
     1    (rho1c(irho+1) - 2.0d0*rho1c(irho) + rho1c(irho-1))/(drho**2)
        end do

! Endpoints
        rhop1c(1) = 2.0d0*rhop1c(2) - rhop1c(3)
        rhop1c(nnrho) = 2.0d0*rhop1c(nnrho-1) - rhop1c(nnrho-2)

        rhopp1c(1) = 2.0d0*rhopp1c(2) - rhopp1c(3)
        rhopp1c(nnrho) = 2.0d0*rhopp1c(nnrho-1) - rhopp1c(nnrho-2)

! Format Statements
! ====================================================================

        return
        end
