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
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
! twocenter_integral.f
! Program Description
! ===========================================================================
!      This code computes the actual integral of the general two-center
! matrix elements of the form <psi1|V(1)|psi2>.  Thus V(1) is located at
! the site of one of the orbitals.  The potential V(1) is something like Vxc
! for the exchange correlation potential, Vna for the neutral atom potential,
! or 1 for the overlap term.
!
! The integral is  performed in cylindrical coordinates over rho and z
! (the phi part having been done by hand and giving the ifactor's below).
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine twocenter_integral (interaction, isorp, ideriv, iexc,
     1                                 fraction, n1, l1, m1, n2, l2, m2,
     2                                 nz, nrho, d, itype1, itype2,
     3                                 rcutoff1, rcutoff2, sum, ispher)
        implicit none
        include '../parameters.inc'
        include '../exchange.inc'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer ideriv
        integer iexc            ! type of exchange-correlation approximation
        integer isorp
        integer itype1
        integer itype2
 
        integer interaction     ! interaction type
                                !   0 = oslxc density
                                !   1 = overlap
                                !   2 = vna ontop
                                !   3 = vna atom
                                !   4 = non-local
                                !   5 = xc ontop
                                !   6 = xc atom-atom
                                !   7 = xc correction
                                !   8 = z-dipole
                                !   9 = y-dipole
                                !  10 = x-dipole
                                !  11 = coulomb
                                !  12 = extended hubbard (n1*n2 dmuxc(n1+n2)/dn)
                                !  13 = extended hubbard spin-polarization 
 
        integer m1, m2          ! m quantum numbers
        integer l1, l2          ! s, p, d, or f
        integer n1, n2          ! which shell is being considered
        integer nrho            ! number of rho-points on grid
        integer nz              ! number of z-points on grid
        logical ispher          ! spherical approx.
 
        real*8 d
        real*8 fraction
        real*8 rcutoff1, rcutoff2
 
! Output
        real*8 sum
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer lmax
        parameter (lmax = 3)
 
        integer max_points
        parameter (max_points = 1001)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer isorp_new
        integer iwftype         ! interaction type
                                ! 0 => atom/ontop wavefunctions different sites
                                ! 1 => atom/atom wavefunctions same site
        integer ix
        integer iz
        integer nnrho
        integer nnz
 
! These are the factors which come from the coefficient in the Ylm
        real*8 clm (0:lmax, -lmax:lmax)
 
        real*8 density1
        real*8 density2
        real*8 density3
        real*8 dexc
        real*8 dpotxc12
        real*8 dpotxc12s
        real*8 drho
        real*8 dvxc
        real*8 dz
        real*8 factor
        real*8 phifactor
        real*8 pi
        real*8 psi1
        real*8 psi2
        real*8 psi3
        real*8 psiofr
        real*8 r1, r2
        real*8 rescaled_psi
        real*8 rho
        real*8 rhomax
        real*8 rhomin
        real*8 vofr
        real*8 vnnaofr
        real*8 vppofr
        real*8 vxc
        real*8 z1
        real*8 z2
        real*8 zmax
        real*8 zmin
 
        real*8 rhomult (max_points)
        real*8 zmult (max_points)
 
        external dexc
        external dpotxc12
        external dpotxc12s
        external dvxc
        external psiofr
        external rescaled_psi
        external vnnaofr
        external vppofr
        external vxc
 
! Procedure
! ===========================================================================
! The variable iwftype defines where the potential is located, and the type
! of integral we are doing; iwftype 0 => ontop, iwftype 1 => atom
        if (interaction .ne. 0 .and. interaction .ne. 3 .and. 
     1      interaction .ne. 4 .and. interaction .ne. 6) then
         iwftype = 0
        else if (interaction .eq. 3 .or. interaction .eq. 6) then
         iwftype = 1
        else if (interaction .eq. 4) then
         iwftype = 2
! oslxc: ontop_right and ontop_left 2c density 
        else if ((interaction .eq. 0) .and. (ideriv.eq.0)) then 
           iwftype = 1
        else if ((interaction .eq. 0) .and. (ideriv.ne.0)) then 
           iwftype = 0
! Mc-Weda: ontop_right and ontop_left 2c density charge transfer correction
        else if ((interaction .eq. 14) .and. (ideriv.eq.0)) then 
           iwftype = 1
        else if ((interaction .eq. 14) .and. (ideriv.ne.0)) then 
           iwftype = 0
        end if

 
! Initialize constants
        pi = 3.141592653589793238462643D0
 
! The Ylm coefficients
        clm(0,0)  = 1.0d0
 
        clm(1,-1) = sqrt(3.0d0/2.0d0)
        clm(1,0)  = sqrt(3.0d0)
        clm(1,1)  = sqrt(3.0d0/2.0d0)
 
        clm(2,-2) = sqrt(15.0d0/8.0d0)
        clm(2,-1) = sqrt(15.0d0/2.0d0)
        clm(2,0)  = sqrt(5.0d0/4.0d0)
        clm(2,1)  = sqrt(15.0d0/2.0d0)
        clm(2,2)  = sqrt(15.0d0/8.0d0)
 
        clm(3,-3) = sqrt(35.0d0/16.0d0)
        clm(3,-2) = sqrt(105.0d0/8.0d0)
        clm(3,-1) = sqrt(21.0d0/16.0d0)
        clm(3,0)  = sqrt(7.0d0/4.0d0)
        clm(3,1)  = sqrt(21.0d0/16.0d0)
        clm(3,2)  = sqrt(105.0d0/8.0d0)
        clm(3,3)  = sqrt(35.0d0/16.0d0)

        if (interaction .eq. 9) then

          clm(0,0)  = sqrt(1.0d0/2.0d0)

          clm(1,-1) = sqrt(3.0d0/2.0d0)
          clm(1,0)  = sqrt(3.0d0/2.0d0)
          clm(1,1)  = sqrt(3.0d0/2.0d0)

          clm(2,-2) = sqrt(15.0d0/32.0d0)
          clm(2,-1) = sqrt(15.0d0/2.0d0)
          clm(2,0)  = sqrt(5.0d0/8.0d0)
          clm(2,1)  = sqrt(15.0d0/2.0d0)
          clm(2,2)  = -sqrt(15.0d0/32.0d0)

        endif

        if (interaction .eq. 10) then

          clm(0,0)  = sqrt(1.0d0/2.0d0)

          clm(1,-1) = sqrt(3.0d0/2.0d0)
          clm(1,0)  = sqrt(3.0d0/2.0d0)
          clm(1,1)  = sqrt(3.0d0/2.0d0)

          clm(2,-2) = sqrt(15.0d0/32.0d0)
          clm(2,-1) = sqrt(15.0d0/2.0d0)
          clm(2,0)  = sqrt(5.0d0/8.0d0)
          clm(2,1)  = sqrt(15.0d0/2.0d0)
          clm(2,2)  = sqrt(15.0d0/32.0d0)

        endif
 
! Initialize the sum to zero
        sum = 0.0d0
 
! Set integration limits
        zmin = dmax1(-rcutoff1, d - rcutoff2)
        zmax = dmin1(rcutoff1, d + rcutoff2)
 
        rhomin = 0.0d0
        rhomax = dmin1(rcutoff1, rcutoff2)
 
        if (interaction .eq. 11 .or. 
     1      (interaction .eq. 3 .and. isorp .ne. 0)) then
         zmin = - rcutoff1
         zmax = rcutoff1
         rhomax = rcutoff1
        end if
 
! Strictly define what the density of the mesh should be.  Make the density of
! the number of points equivalent for all cases. Change the number of points
! to be integrated to be dependent upon the distance between the centers and
! this defined density.
        dz = (rcutoff1 + rcutoff2)/(dfloat(nz)*2.0d0)
        nnz = int((zmax - zmin)/dz)
        if (mod(nnz,2) .eq. 0) nnz = nnz + 1
 
        drho = max(rcutoff1,rcutoff2)/dfloat(nrho)
        nnrho = int((rhomax - rhomin)/drho)
        if (mod(nnrho,2) .eq. 0) nnrho = nnrho + 1
 
! Set up Simpson's rule factors. First for the z integration and then for
! the rho integration.
        zmult(1) = dz/3.0d0
        zmult(nnz) = dz/3.0d0
        do iz = 2, nnz - 1, 2
         zmult(iz) = 4.0d0*dz/3.0d0
        end do
        do iz = 3, nnz - 2, 2
         zmult(iz) = 2.0d0*dz/3.0d0
        end do
 
        rhomult(1) = drho/3.0d0
        rhomult(nnrho) = drho/3.0d0
        do irho = 2, nnrho - 1, 2
         rhomult(irho) = 4.0d0*drho/3.0d0
        end do
        do irho = 3, nnrho - 2, 2
         rhomult(irho) = 2.0d0*drho/3.0d0
        end do
 
! This factor comes from a result of multiplying the two Ylm's together.
! and after the factor of pi is multiplied out after the phi integration.
! For the coulomb case, we need spherically symmetric charge densities,
! so add up all m's squared
        if (interaction .eq. 7 .or. interaction .eq. 11
     1        .or. interaction .eq. 12 .or. interaction .eq. 13) then
         phifactor = 2.0d0*pi
        else
         phifactor = clm(l1,m1)*clm(l2,m2)/2.0d0
        end if
 
! Integration is over z (z-axis points from atom 1 to atom 2) and rho (rho is
! radial distance from z-axis).
        do iz = 1, nnz
         z1 = zmin + dfloat(iz-1)*dz
         z2 = z1 - d
         do irho = 1, nnrho
          rho = rhomin + dfloat(irho-1)*drho
          r1 = sqrt(z1**2 + rho**2)
          r2 = sqrt(z2**2 + rho**2)
 
          if (r1 .lt. rcutoff1) then
 
! Total integration factor
           factor = zmult(iz)*rhomult(irho)*phifactor
 
! Calculate psi for each of the two atoms. Again, the variables l1, l2 (= 0-3)
! means s, p, d, and f-state. The variables m1, m2 (= 0-3) implies sigma,
! pi, delta, or phi.
 
! Given the position r, first get radial part of psi evaluated at this r
! The wavefunction on the left ("bra") is multiplied by the potential
! which is located at the same site as this orbital.
! For most cases (interaction .ne. 3, 4, 6, 11) the two orbitals are located at
! at different sites (iwftype = 0, but for the atom-atom case they are located
! at the same site (iwftype = 1).
           if (iwftype .eq. 0) then
            psi1 = psiofr (itype1, n1, r1)      ! v (ontop) case
            psi2 = psiofr (itype2, n2, r2)
           else if (iwftype .eq. 1) then        ! v (atom) case
            psi1 = psiofr (itype1, n1, r1)
            psi2 = psiofr (itype1, n2, r1)
           else if (iwftype .eq. 2) then        ! vnlocal case
            psi1 = psiofr (itype1, n1, r1)
            psi2 = vppofr (itype2, n2, r2)
           end if

! jel-spher
           if(ispher) then 
             psi1 = sqrt( psi1**2.0d0 )
             psi2 = sqrt( psi2**2.0d0 )
           endif
! end jel-spher
 
           density1 = psi1**2/(4.0d0*pi)
           density2 = psi2**2/(4.0d0*pi)

! comment: should be added later (oslxc-harris option)
!           if (interaction .eq. 0) then
!            if (isorp .eq. 0) then
!             density3 = 0.0d0
!             do isorp_new = 1, nsshxc(itype2) 
!              psi3 = psiofr (itype2, isorp_new, r2)
!              density3 = 
!     1         density3 + xnocc(isorp_new,itype2)*psi3**2/(4.0d0*pi)
!             end do
!            else
!             psi3 = psiofr (itype2, isorp, r2)
!             density3 = psi3**2/(4.0d0*pi)
!            end if
!           end if
 
! *************************************************************************
! ontop => orbitals at two different sites
! atom-atom => orbitals at same site
! interaction = 0 => do oslxc density
! interaction = 1 => do overlap
! interaction = 2 => do neutral atom/ontop function => vnnaofr
! interaction = 3 => do neutral atom/atom
! interaction = 4 => do non-local
! interaction = 5 => do xc ontop function => vxc
! interaction = 6 => do xc atom function => dvxc
! interaction = 7 => do xc double counting correction function => dexc
! interaction = 8 => do z-dipole
! interaction = 9 => do y-dipole
! interaction = 10 => do x-dipole
! interaction = 11 => do coulomb
! interaction = 12 => do extended hubbard (n1*n2*dmuxc(n1+n2), dmuxc=dmuxc/dn).
! interaction = 13 => do extended hubbard spin dependent contribution
! interaction = 13 => do charge transfer correction
! *************************************************************************
           ix = ideriv + 1

           if (interaction .eq. 0 .and. ideriv .eq. 0) then
             psi3 = psiofr (itype2, isorp, r2)
             ! density psi**2/const
             vofr = psi3**2/(4.0d0*pi)
           else if (interaction .eq. 0 .and. ideriv .eq. 1) then
             psi3 = psiofr (itype1, isorp, r1)
             ! density psi**2/const
             vofr = psi3**2/(4.0d0*pi)
           else if (interaction .eq. 0 .and. ideriv .eq. 2) then
             psi3 = psiofr (itype2, isorp, r2)
             ! density psi**2/const
             vofr = psi3**2/(4.0d0*pi)
           else if (interaction .eq. 1 .or. interaction .eq. 4) then
            vofr = 1.0d0
           else if (interaction .eq. 2) then
            if (ideriv .eq. 1) vofr = vnnaofr (itype1, isorp, r1)
            if (ideriv .eq. 2) vofr = vnnaofr (itype2, isorp, r2)
           else if (interaction .eq. 3) then
            vofr = vnnaofr (itype2, isorp, r2)
           else if (interaction .eq. 5) then
            vofr = vxc (rho, z1, iexc, fraction)
           else if (interaction .eq. 6) then
            vofr = dvxc (itype1, itype2, rho, z1, r1, iexc, fraction,
     1                   ix)
           else if (interaction .eq. 7) then
            vofr = dexc (itype1, itype2, rho, z1, r1, r2, iexc,
     1                   fraction, ix)
           else if (interaction .eq. 8) then
            vofr = z1 - 0.5d0*d
           else if (interaction .eq. 9) then
            vofr = rho
           else if (interaction .eq. 10) then
            vofr = rho
           else if (interaction .eq. 11) then
            vofr = vnnaofr (itype2, n2, r2)

           else if (interaction .eq. 12) then
            vofr = dpotxc12 (rho, z1, iexc, fraction)
           else if (interaction .eq. 13) then
            vofr = dpotxc12s (rho, z1, iexc, fraction) 
! McWeda charge transfer correction
           else if (interaction .eq. 14 .and. ideriv .eq. 1) then
	     vofr = dpotxc12 (rho, z1, iexc, fraction)
             psi3 = psiofr (itype1, isorp, r1)
             ! density psi**2/const
             vofr = psi3**2/(4.0d0*pi)*vofr
           else if (interaction .eq. 14 .and. ideriv .eq. 2) then
	     vofr = dpotxc12 (rho, z1, iexc, fraction)
             psi3 = psiofr (itype2, isorp, r2)
             ! density psi**2/const
             vofr = psi3**2/(4.0d0*pi)*vofr
           end if
 
! Let's skip over all this stuff if we are NOT using wavefunctions.
! There may be other cases where we can skip all this stuff also.
           if (interaction .ne. 7 .and. interaction .ne. 11 .and.
     1         interaction .ne. 12. and .interaction .ne. 13) then
 
! *************************************************************************
! Add magic factors based on what type of orbital is involved in the integration
            psi1 = rescaled_psi (l1, m1, rho, r1, z1, psi1)
 
! If iwftype = 0, then the wavefunctions are at different sites, and r1 is
! different from r2.  If iwftype = 1, then r1 needs to be set equal to r2,
! just for this wavefunction part, since they are located at the same site.
            if (iwftype .eq. 1) then
             psi2 = rescaled_psi (l2, m2, rho, r1, z1, psi2)
            else
             psi2 = rescaled_psi (l2, m2, rho, r2, z2, psi2)
            end if
 
! For double counting correction to the exchange-correlation interaction,
! set the wavefunctions to 1.0d0. This term is not a matrix element, but
! rather a correction to the over counting which occured in the matrix
! elements evaluation.
           else if (interaction .eq. 7) then
            psi1 = 1.0d0
            psi2 = 1.0d0
           else if (interaction .eq. 11) then
            psi1 = density1
            psi2 = 1.0d0
           else if (interaction .eq. 12 .or. interaction .eq. 13) then
            psi1 = density1
            psi2 = density2
           end if
 
! This is the actual integral
           sum = sum + factor*psi1*vofr*psi2*rho

          end if
         end do
        end do
 
! Format Statements
! ===========================================================================
2222    format(I4,'  ', f10.7,'   ', f12.4)

        return
        end

