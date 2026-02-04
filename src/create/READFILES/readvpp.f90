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
 
 
! readvpp.f
! Program Description
! ===========================================================================
!        This subroutine readvpp reads the values for the pseudopotential
! of the different (s,p,d,f) orbitals from the data file atomicsymbol.pp
! which comes from a pseudopotential generator.
! Note: The potential is already in Angstrom units on the data file.
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
        subroutine readvpp (ispec, filein, nsshPP, lsshPP, iexc, &
     &                      fraction, iammaster)
        use precision, only: wp
        implicit none
 
        include '../parameters.inc'
        include '../pseudopotentials.inc'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer ispec
 
        character(len=25) filein
 
        logical iammaster
 
! Output
        integer iexc
 
        integer nsshPP (nspec_max)
        integer lsshPP (nspec_max, nsh_max)
 
        real(kind=wp) fraction
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iline
        integer ipoint
        integer issh
        integer mesh
 
        real(kind=wp) alpha
        real(kind=wp) rshort
        real(kind=wp) vshort

        real(kind=wp) Z_val

        real(kind=wp) rcutoff
 
! Procedure
! ===========================================================================
! Open the input file
        open (unit = 88, file = filein, status = 'old')
 
        if (iammaster) then
         ! write (*,*) '  '
         ! write (*,*) ' *---------------------------------------------* '
         ! write (*,*) ' |               Welcome to READVPP            | '
         ! write (*,*) ' | Reading (non)-neutral potential of the atom | '
         ! write (*,*) ' *---------------------------------------------* '
         ! write (*,*) '  '
         ! write (*,91) filein
        end if ! end master
 
! There are 14 message lines in each pseudopotential file
        do iline = 1, 14
         read (88,*)
        end do
 
! Read in which exchange-correlation approximation we are using.
        read (88,*) iexc
        fraction = 0.0d0
        if (iexc .eq. 12) then
         rewind(88)
         do iline = 1, 14
          read (88,*)
         end do
         read (88,*) iexc, fraction
        end if
 
        if (iammaster) then
        ! write (*,*) ' We are reading in iexc from the pseudopotential '
        ! write (*,*) ' file. This tells us which exchange-correlation '
        ! write (*,*) ' approximation we are using. '
        ! write (*,*) ' You have chosen iexc = ', iexc
        ! write (*,*) ' The options that available are: '
        ! write (*,*) ' 1  LDA Wigner'
        ! write (*,*) ' 2  LDA Hedin/Lundqvist'
        ! write (*,*) &
     ! &  ' 3  LDA Ceperley/Alder Perdew/Zunger (1980) *** default'
        ! write (*,*) ' 4  GGA Perdew/Wang (1991)'
        ! write (*,*) ' 5  GGA Becke (1988) X, Perdew (1986)'
        ! write (*,*) ' 6  GGA Perdew/Burke/Ernzerhof (1996)'
        ! write (*,*) ' 7  LDA Zhao/Parr'
        ! write (*,*) ' 8  LDA Ceperley/Alder Perdew/Wang (1991)'
        ! write (*,*) ' 9  GGA Becke (1988) X, Lee/Yang/Parr (1988)'
        ! write (*,*) ' 10 GGA Perdew/Wang (1991) X, Lee/Yang/Parr (1988)'
        ! write (*,*) ' 11 LDA exchange only'
        ! write (*,*) ' 12 GGA Becke (1988) X, Lee/Yang/Parr (1988), but '
        ! write (*,*) '    with mixing of exact exchange. '
        end if ! end master
 
! Read the number of shells
        read (88,*) nsshPP(ispec)
        read (88,*) (lsshPP(ispec,issh), issh = 1, nsshPP(ispec))
        if (iammaster) then
         ! write (*,*) '  '
         ! write (*,*) ' # of pseudopotential shells = ', nsshPP(ispec)
         ! write (*,*) ' The L values = ', &
     ! &    (lsshPP(ispec,issh), issh = 1, nsshPP(ispec))
        end if ! end master
 
! jel-PP
! Read Z_val
        read (88,*) Z_val

! Read in alpha for longranged local part => -Z*e**2*erf(alpha*r)/r
        read (88,*) alpha

! jel-PP
! Read cutoff radius of PP
        read (88,*) rcPP(ispec)
 
! Read in the short-range local part  - this is not needed for the crtor
        read (88,*) mesh
        do ipoint = 1, mesh
         read (88,*) rshort, vshort
        end do
 
! Read in the pseudopotential - this is not needed for the crtor
        do issh = 1, nsshPP(ispec)
         read (88,100) mesh
         do ipoint = 1, mesh
          read (88,*) rshort, vshort
         end do
        end do
 
! Now read in the points for the non-local part
        do issh = 1, nsshPP(ispec)
         read (88,200) mesh, cl_pp(issh,ispec)
         if (mesh .gt. max_points_pp) then
          ! write (*,*) ' max_points_pp = ', max_points_pp, &
     ! &                ' mesh = ', mesh
          ! write (*,*) ' The number of mesh points in your file is '
          ! write (*,*) ' greater than the dimensioned number of points.'
          ! write (*,*) ' Double check everything and rerun creator.'
          stop 'error in readvpp'
         end if
 
         npoints_pp(issh,ispec) = mesh
         do ipoint = 1, mesh
          read (88,*) rr_pp(ipoint,issh,ispec), vpp(ipoint,issh,ispec)
         end do
         rrc_pp(issh,ispec) = rr_pp(mesh,issh,ispec)
         drr_pp(issh,ispec) = rr_pp(2,issh,ispec) - rr_pp(1,issh,ispec)
        end do
 
        close (unit = 88)

        do issh = 1, nsshPP(ispec)
         rcutoff = (mesh-1) * drr_pp(issh,ispec)
         if (superspline) call buildspline_1d (vpp(1,issh,ispec),    &
     &                        vpp_spline(1,1,issh,ispec), mesh, rcutoff)
        end do

 
      if (iammaster) then
       ! write (*,*) '  '
       ! write (*,*) ' *---------------- END READVPP -------------------*'
       ! write (*,*) '  '
      end if ! end master
 
! Code to just make ftnchek happy, because these are not ever used
        if (alpha*rshort*vshort .eq. 0) return
 
 
! Format Statements
! ===========================================================================
91      format (2x, a25)
100     format (12x, i5)
200     format (12x, i5, 4x, f14.7)
 
        return
        end
