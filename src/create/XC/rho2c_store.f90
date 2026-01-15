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
! rho2c_store.f
! Program Description
! ====================================================================
!
! This routine calculates and stores the combined density of species
! 1 and 2 as a function of r, z and d.
!
! On input:
!            iexc:       The exchange-correlation option used
!            itype1:     The species type of atom 1
!            itype2:     The species type of atom 2
!            rcutoff1:   The radius cutoff of atom 1
!            rcutoff2:   The radius cutoff of atom 2
!
!     By common blocks:
!            nsshxc:     Number of shells for each species
!
! On output: rho2c:      The density as a function of species type,
!                        r, z and d.  Output is placed in common block
!                        density located in wavefunctions.inc
!            rhop2c:     derivative with respect to rho
!            rhopp2c:    second derivative with respect to rho
!            rhoz2c:     derivative with respect to z
!            rhozz2c:    second derivative with respect to rho
!            rhopz2c:    mixed derivative with respect to rho and z
!
! ====================================================================
! Code written by:
! Richard B. Evans
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
        subroutine rho2c_store (iexc, itype1, itype2, rcutoff1, &
     &                          rcutoff2, d, ix)
        use precision, only: wp
        implicit none

        include '../parameters.inc'
        include '../exchange.inc'
        include '../quadrature.inc'
        include '../wavefunctions.inc'

! Input Declaration and Description
! ====================================================================
! Input
        integer iexc
        integer ix
        integer itype1
        integer itype2

        real(kind=wp) d
        real(kind=wp) rcutoff1
        real(kind=wp) rcutoff2

! Local Parameters and Data Declaration
! ====================================================================
        integer j1(5)
        data j1 /1, 1, 1, 2, 2/

        integer j2(5)
        data j2 /0, -1, 1, -1, 1/

! Local Variable Declaration and Description
! ====================================================================
        integer irho
        integer issh
        integer iz
        integer j1at
        integer j1ch
        integer jssh

        integer in (2)

        real(kind=wp) dens
        real(kind=wp) dzraw
        real(kind=wp) psiofr
        real(kind=wp) r
        real(kind=wp) r1
        real(kind=wp) r2
        real(kind=wp) rho
        real(kind=wp) xinv4pi
        real(kind=wp) z1
        real(kind=wp) z2

        external psiofr

! Procedure
! ====================================================================
! Initialize 1/4pi
        xinv4pi = 1.0d0/(4.0d0*3.141592653589793238462643D0)

! First initialize in (index) array. This is for the charge correction
! calculated later.
        in(1) = itype1
        in(2) = itype2

        j1ch = j1(ix)
        j1at = in(j1ch)
        jssh = iderorb(j1at)

! Fix the endpoints and initialize the increments dz and drho.
        zmin = dmin1(-rcutoff1, d - rcutoff2)
        zmax = dmax1(rcutoff1, d + rcutoff2)

        rhomin = 0.0d0
        rhomax = dmax1(rcutoff1, rcutoff2)

        dzraw = dmin1(drr_rho(itype1),drr_rho(itype2))
        dz = dzraw*ixcgridfactor

! If the grid size is made too large or too small, then set to defaults.
        if (dz .gt. 0.05d0) dz = 0.05d0
        if (dz .lt. 0.002d0) dz = 0.002d0

        nnz = int((zmax - zmin)/dz) + 1

        drho = dz
        nnrho = int((rhomax - rhomin)/drho) + 1

! Performs some checks to make sure that the dimensions on rho2c are
! alright.
        if (nnrho .gt. nrho_points) then
         write (*,*) ' In rho2c_store.f, nnrho = ', nnrho
         write (*,*) ' The dimension nrho_points =', nrho_points
         write (*,*) ' needs to be increased. '
         stop 'error in rho2c_store'
        end if

        if (nnz .gt. nz_points) then
         write (*,*) ' In rho2c_store.f, nnz = ', nnz
         write (*,*) ' The dimension nz_points =', nz_points
         write (*,*) ' needs to be increased. '
         stop 'error in rho2c_store'
        end if

! **************************************************************************
! Here we loop over z and r computing the sum of the densities
! for species in1 and in2 at each value of d, z and r.
        do iz = 1, nnz
         z1 = zmin + (iz-1)*dz
         z2 = z1 - d

         do irho = 1, nnrho
          rho = rhomin + (irho-1)*drho
          r1 = sqrt(z1**2 + rho**2)
          r2 = sqrt(z2**2 + rho**2)

          dens = 0.0d0
          do issh = 1, nsshxc(itype1)
           dens = dens + xnocc(issh,itype1)*psiofr(itype1,issh,r1)**2
          end do
          do issh = 1, nsshxc(itype2)
           dens = dens + xnocc(issh,itype2)*psiofr(itype2,issh,r2)**2
          end do

! Here the derivative with respect to the charge correction term
! is calculated.
          if (j1ch .eq. 1) r = r1
          if (j1ch .eq. 2) r = r2

          dens = dens +j2(ix)*dqorb(j1at)*psiofr(j1at,jssh,r)**2
          rho2c(irho,iz) = dens*xinv4pi
         end do
        end do

! **************************************************************************
! Now calculate the derivatives
! Only calculate the derivatives if doing GGA exchange-correlation.
        if (iexc .eq. 4 .or. iexc .eq. 5 .or. iexc .eq. 6 &
     &      .or. iexc .eq. 9 .or. iexc .eq. 10) then

! **************************************************************************
! Calculate rhop2c and rhopp2c.
         do iz = 1, nnz
          do irho = 2, nnrho - 1

! First derivative:
           rhop2c(irho,iz) = &
     &      (rho2c(irho+1,iz) - rho2c(irho-1,iz))/(2.0d0*drho)

! Second derivative:
           rhopp2c(irho,iz) = &
     &      (rho2c(irho+1,iz) - 2.0d0*rho2c(irho,iz) &
     &       + rho2c(irho-1,iz))/(drho**2)
          end do

! Find endpoint values for the derivatives calculated above.
          rhop2c(1,iz) = 2.0d0*rhop2c(2,iz) - rhop2c(3,iz)
          rhop2c(nnrho,iz) = 2.0d0*rhop2c(nnrho-1,iz) &
     &                      - rhop2c(nnrho-2,iz)

          rhopp2c(1,iz) = 2.0d0*rhopp2c(2,iz) - rhopp2c(3,iz)
          rhopp2c(nnrho,iz) = 2.0d0*rhopp2c(nnrho-1,iz) &
     &                       - rhopp2c(nnrho-2,iz)
         end do

! **************************************************************************
! Calculate rhoz2c and rhozz2c.
         do irho = 1, nnrho
          do iz = 2, nnz - 1

! First derivative:
           rhoz2c(irho,iz) = &
     &      (rho2c(irho,iz+1) - rho2c(irho,iz-1))/(2.0d0*dz)

! Second derivative:
           rhozz2c(irho,iz) = &
     &      (rho2c(irho,iz+1) - 2.0d0*rho2c(irho,iz) &
     &       + rho2c(irho,iz-1))/(dz**2)
          end do

! Find enpoint values for the derivatives calculated above.
          rhoz2c(irho,1) = 2.0d0*rhoz2c(irho,2) - rhoz2c(irho,3)
          rhoz2c(irho,nnz) = 2.0d0*rhoz2c(irho,nnz-1) &
     &                      - rhoz2c(irho,nnz-2)

          rhozz2c(irho,1) = 2.0d0*rhozz2c(irho,2) - rhozz2c(irho,3)
          rhozz2c(irho,nnz) = 2.0d0*rhozz2c(irho,nnz-1) &
     &                       - rhozz2c(irho,nnz-2)
         end do

! **************************************************************************
! Now calculate the cross term derivatives - rhopz2c.
         do irho = 1, nnrho
          do iz = 2, nnz - 1
           rhopz2c(irho,iz) = &
     &      (rhop2c(irho,iz+1) - rhop2c(irho,iz-1))/(2.0d0*dz)
          end do

! Now calculate the derivatives at the endpoints.
          rhopz2c(irho,1) = 2.0d0*rhopz2c(irho,2) - rhopz2c(irho,3)
          rhopz2c(irho,nnz) = 2.0d0*rhopz2c(irho,nnz-1) &
     &                       - rhopz2c(irho,nnz-2)
         end do
        end if

! Format Statements
! ====================================================================

        return
        end