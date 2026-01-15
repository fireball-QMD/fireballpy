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
 
! readpsi.f
! Program Description
! ===========================================================================
!        This subroutine readpsi reads the values for the radial wavefunctions
! of the different (s,p,d,f) orbitals off the data file atomicsymbol.wfl
! which comes from rcatms.f. The normalization will also be checked here.
! Note: Psi is already in Angstrom units on the data file.
! The variable ispec is a flag which would indicate in the multi-species
! MD code which species to read.
 
! rc must be input to this subroutine in abohr units.
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
 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readpsi (ispec, issh, lqn, rcutoff, xnoccin, nzx, &
     &                      filein, iammaster)
        use precision, only: wp
        implicit none
 
        include '../parameters.inc'
        include '../wavefunctions.inc'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispec
        integer issh
        integer lqn
        integer nzx
 
        real(kind=wp) rcutoff
        real(kind=wp) xnoccin
 
        character(len=25) filein
 
        logical iammaster
 
! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=wp) abohr
        parameter (abohr = 0.529177249d0)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer inum
        integer ipoint
        integer iremainder
        integer lqnwf
        integer mesh
        integer nzxwf
 
        real(kind=wp) r
        real(kind=wp) rc
        real(kind=wp) rc_max
        real(kind=wp) rcutoffwf
        real(kind=wp) sum
        real(kind=wp) xnoccwf
 
        real(kind=wp) psitemp(wfmax_points)
 
        character(len=25) fileinwf
 
! Procedure
! ===========================================================================
! Open the input file
        open (unit = 88, file = filein, status = 'old')
 
        if (iammaster) then
         write (*,*) ' '
         write (*,*) '*-----------------------------------------------*'
         write (*,*) '|               Welcome to READPSI              |'
         write (*,*) '| Reading the radial wavefunction of your atom  |'
         write (*,*) '*-----------------------------------------------*'
         write (*,*) ' '
        end if ! end master
 
        read (88,90) fileinwf
        read (88,*) nzxwf
        read (88,*) mesh
        read (88,*) rcutoffwf, rc_max, xnoccwf
        read (88,*) lqnwf
        if (iammaster) write (*,91) fileinwf
 
! Perform some checks
        if (nzxwf .ne. nzx) then
         write (*,*) ' nzxwf = ', nzxwf, ' nzx = ', nzx
         write (*,*) ' The Z number the wavefunction file, for this '
         write (*,*) ' shell, does not match the cutoff radius '
         write (*,*) ' that you put into the create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readpsi'
        end if
 
        if (rcutoffwf .le. (rcutoff - 1.0d-2) .or. &
     &      rcutoffwf .ge. (rcutoff + 1.0d-2)) then
         write (*,*) ' rcutoffwf = ', rcutoffwf, ' rcutoff = ', rcutoff
         write (*,*) ' The cutoff radius in the wavefunction file, for '
         write (*,*) ' this shell, does not match the cutoff radius '
         write (*,*) ' that you put into your create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readpsi'
        end if
 
        if (xnoccwf .ne. xnoccin) then
         write (*,*) ' xnoccwf = ', xnoccwf, ' xnoccin = ', xnoccin
         write (*,*) ' The occupation number in the wavefunction file, '
         write (*,*) ' for this shell, does not match the occupation '
         write (*,*) ' number that you put into your create.input file.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readpsi'
        end if
 
        if (lqnwf .ne. lqn) then
         write (*,*) ' lqnwf = ', lqnwf, ' lqn = ', lqn
         write (*,*) ' The l quantum number in the wavefunction file, '
         write (*,*) ' for this shell, does not match the l quantum '
         write (*,*) ' number that you put into your create.input file.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readpsi'
        end if
 
        if(mesh .gt. wfmax_points) then
         write (*,*) ' Error error ***** in readpsi. '
         write (*,*) ' Dimension of wavefunction = ', wfmax_points
         write (*,*) ' We are asking for mesh = ', mesh
         write (*,*) ' Redimension wfmax_points. '
         stop 'error in readpsi'
        end if
 
! Set some things up
        rc = rcutoffwf*abohr
 
        if (iammaster) then
         write (*,*) '  '
         write (*,200) rcutoffwf, rc
        end if ! end master
        npoints(issh,ispec) = mesh
        rrc(issh,ispec) = rc
        if (issh .eq. 1) rrc_rho(ispec) = -1.0d0
        rrc_rho(ispec) = max(rrc_rho(ispec),rc)
 
        drr(issh,ispec) = rc/real(mesh - 1, kind=wp)
        if (issh .eq. 1) drr_rho(ispec) = 99.0d0
        drr_rho(ispec) = min(drr_rho(ispec),drr(issh,ispec))
 
        npoints_rho(ispec) = int(rrc_rho(ispec)/drr_rho(ispec)) + 1
        rrc_rho(ispec) = real(npoints_rho(ispec) - 1, kind=wp)*drr_rho(ispec)
 
! Shift the rrc_rho value just a little or else the exchange-correlation
! interactions for gradients go haywire at the endpoints where rho = 0.0d0
        rrc_rho(ispec) = rrc_rho(ispec) - 1.0d-4
 
        do ipoint = 1, npoints_rho(ispec)
         rr_rho(ipoint,ispec) = real(ipoint - 1, kind=wp)*drr_rho(ispec)
        end do
 
! Read in the points
        inum = idint(real(mesh, kind=wp)/4)
        iremainder = mesh - (inum*4)
        do ipoint = 1, mesh - iremainder, 4
         read (88,100) psitemp(ipoint), psitemp(ipoint+1), &
     &                 psitemp(ipoint+2), psitemp(ipoint+3)
        end do
 
        if (iremainder .eq. 1) then
         read (88,100) psitemp(mesh)
        else if (iremainder .eq. 2) then
         read (88,100) psitemp(mesh-1), psitemp(mesh)
        else if (iremainder .eq. 3) then
         read (88,100) psitemp(mesh-2), psitemp(mesh-1), psitemp(mesh)
        end if
 
        close (unit = 88)
 
! Write psitemp to the wavefunction psi
        do ipoint = 1, mesh
         psi(ipoint,issh,ispec) = psitemp(ipoint)
        end do
 
! Check normalization
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Checking normalization [NORM(l) should be 1]'
         write (*,*) '  '
        end if ! end master
 
        r = - drr(issh,ispec)
        do ipoint = 1, mesh
         r = r + drr(issh,ispec)
         rr(ipoint,issh,ispec) = r
        end do
 
        sum = 0.0d0
        do ipoint = 1, mesh
         if (ipoint .ne. 1 .or. ipoint .ne. mesh) then
          sum = sum + drr(issh,ispec)*rr(ipoint,issh,ispec)**2 &
     &                               *psi(ipoint,issh,ispec)**2
         else
          sum = sum + 0.5d0*drr(issh,ispec)*rr(ipoint,issh,ispec)**2 &
     &                                     *psi(ipoint,issh,ispec)**2
         end if
        end do

        if (superspline) call buildspline_1d (psi(1,issh,ispec),  &
     &                       psi_spline(1,1,issh,ispec), mesh, rcutoff)
 
        if (iammaster) then
         write (*,300) issh, sum
         write (*,*) '  '
         write (*,*) ' *--------------- END READPSI ------------------*'
         write (*,*) '  '
        end if ! end master
 
! Code to just make ftnchek happy, because rc_max is not ever used
        if (rc_max .eq. 0.0d0) return
 
! Format Statements
! ===========================================================================
90      format (a25)
91      format (2x,a25)
100     format (4d18.10)
200     format (' Rc read in =', f8.4, ' abohr = ', f8.4, ' Angstroms')
300     format (' NORM (shell = ', i1, ') = ', f16.12)
 
        return
        end