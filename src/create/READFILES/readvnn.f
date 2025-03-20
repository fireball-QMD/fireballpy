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
 
 
 
! readvnn.f
! Program Description
! ===========================================================================
!        This subroutine readvnn reads the values for the (non)-neutral
! potential of the different (s,p,d,f) orbitals from the data file
! atomicsymbol.(spdf)_na which comes from navs.f
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
        subroutine readvnn (ispec, issh, rcutoff, nzx, filein, etotatom,
     1                      iammaster)
        implicit none
 
        include '../parameters.inc'
        include '../vnonneutral.inc'
 
! Argument Declaration and Description
! ===========================================================================
        integer ispec
        integer issh
        integer nzx
 
        real*8 etotatom
        real*8 rcutoff
 
        character*25 filein
 
        logical iammaster
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer mesh
        integer nzxvnn
 
        real*8 rcutoffvnn
 
        character*25 fileinvnn
 
! Procedure
! ===========================================================================
! Open the input file
        open (unit = 88, file = filein, status = 'old')
 
        read (88,90) fileinvnn
        read (88,*) nzxvnn
        read (88,*) rcutoffvnn
        read (88,*) mesh
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) '*-----------------------------------------------*'
         write (*,*) '|               Welcome to READVNN              |'
         write (*,*) '|  Reading (non)-neutral potential of the atom  |'
         write (*,*) '*-----------------------------------------------*'
         write (*,*) '  '
         write (*,91) fileinvnn
        end if ! end master
 
! We now add the total energy of the atom to the vna data file
! (Only na, not charged parts)
        if (issh .eq. 0) then
         read (88,*) etotatom
         if (iammaster) write (*,*) ' Atomic total energy for species:',
     1                   ispec,' is',etotatom
        end if
 
! Perform some checks
        if (mesh .gt. max_points_na) then
         write (*,*) ' max_points_na = ', max_points_na,
     1               ' mesh = ', mesh
         write (*,*) ' The number of mesh points in your file is '
         write (*,*) ' greater than the dimensioned number of points.'
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
        end if
 
        if (nzxvnn .ne. nzx) then
         write (*,*) ' nzxvnn = ', nzxvnn, ' nzx = ', nzx
         write (*,*) ' The cutoff radius in the wavefunction file, '
         write (*,*) ' for this shell, does not match the cutoff '
         write (*,*) ' radius that you put into the create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
        end if
 
        if (abs(rcutoffvnn - rcutoff) .gt. 1.0d-5) then
         write (*,*) ' rcutoffvnn = ', rcutoffvnn,
     1               ' rcutoff = ', rcutoff
         write (*,*) ' The cutoff radius in the wavefunction file, for '
         write (*,*) ' this shell, does not match the cutoff radius '
         write (*,*) ' that you put into your create.input file. '
         write (*,*) ' Double check everything and rerun creator.'
         stop 'error in readvnn'
        end if
 
! Read in the points
        npoints_na(issh,ispec) = mesh
        do ipoint = 1, mesh
         read (88,100) rr_na(ipoint,issh,ispec), vnna(ipoint,issh,ispec)
        end do
        drr_na(issh,ispec) = rr_na(2,issh,ispec) - rr_na(1,issh,ispec)

        if (superspline) call buildspline_1d (vnna(1,issh,ispec),   
     1                     vnna_spline(1,1,issh,ispec), mesh, rcutoff) 
 
        close (unit = 88)
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' *---------------- END READVNN -----------------*'
         write (*,*) '  '
        end if ! end master
 
! Format Statements
! ===========================================================================
90      format (a25)
91      format (2x,a25)
100     format (2d24.16)
 
        return
        end
