! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! readpsi.f90
! Program Description
! ===========================================================================
!        This subroutine readpsi reads the values for the radial wavefunctions
! of the different (s,p,d,f) orbitals off the data file atomicsymbol.wfl
! which comes from rcatms.f. The normalization will also be checked here.
! Note: Psi is already in Angstrom units on the data file.

! rc must be input to this subroutine in abohr units.
! ===========================================================================
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
        subroutine readpsi (nssh, filein, mesh_max, npoints, rrc, drr,      &
     &                      xnocc, rr, psi)
        use constants
        use precision
        use begin_input, only: outpath
        implicit none

! Argument Declaration and Description
! ===========================================================================
!Input
        integer, intent(in) :: mesh_max
        integer, intent(in) :: nssh

        character(len=11), intent(in), dimension (nssh) :: filein

! Output
        integer, intent(out), dimension (nssh) :: npoints

        real(kind=long), intent(out), dimension (nssh) :: drr
        real(kind=long), intent(out), dimension (nssh, mesh_max) :: psi
        real(kind=long), intent(out), dimension (nssh, mesh_max) :: rr
        real(kind=long), intent(out), dimension (nssh) :: rrc
        real(kind=long), intent(out), dimension (nssh) :: xnocc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer inum
        integer ipoint
        integer iremainder
        integer issh
        integer mesh
        integer nzx

        real(kind=long) r
        real(kind=long) rc
        real(kind=long) rc_max
        real(kind=long) sum

        real(kind=long), dimension (:), allocatable :: psitemp

        character(len=15) fileinwf

! Allocate Arrays
! ===========================================================================
        allocate (psitemp(mesh_max))

! Procedure
! ===========================================================================
        !write (*,*) ' '
        !write (*,*) '*--------------------------------------------------*'
        !write (*,*) '|                 Welcome to READPSI               |'
        !write (*,*) '| Reading now the radial wavefunction of your atom |'
        !write (*,*) '*--------------------------------------------------*'
        !write (*,*) ' '

! Open the input file
        do issh = 1, nssh
         open (unit = 21, file = trim(outpath)//trim(filein(issh)), status = 'old')

         read (21,90) fileinwf
         !write (*,90) fileinwf
         read (21,*) nzx
         read (21,*) mesh
         read (21,*) rrc(issh), rc_max, xnocc(issh)
         read (21,*)

! Set some things up
         rc = rrc(issh)*abohr
         npoints(issh) = mesh

         !write (*,*) '  '
         !write (*,200) rrc(issh), rc
         drr(issh) = rc/real(mesh - 1)

! Read in the points
         inum = int(real(mesh)/4.0d0)
         iremainder = mesh - (inum*4)
         psitemp = 0.0d0
         do ipoint = 1, mesh - iremainder, 4
          read (21,100) psitemp(ipoint), psitemp(ipoint+1),                 &
     &                  psitemp(ipoint+2), psitemp(ipoint+3)
         end do

         if (iremainder .eq. 1) then
          read (21,100) psitemp(mesh)
         else if (iremainder .eq. 2) then
          read (21,100) psitemp(mesh-1), psitemp(mesh)
         else if (iremainder .eq. 3) then
          read (21,100) psitemp(mesh-2), psitemp(mesh-1), psitemp(mesh)
         end if

         close (unit = 21)

! Write psitemp to psi - the real wavefunction for all orbitals
         psi(issh,1:mesh) = psitemp(1:mesh)
         psi(issh,mesh+1:mesh_max) = 0.0d0

! Check normalization
         !write (*,*) '  '
         !write (*,*) ' Now checking normalization [NORM(l) should be 1]:'
         !write (*,*) '  '

         r = - drr(issh)
         do ipoint = 1, mesh
          r = r + drr(issh)
          rr(issh,ipoint) = r
         end do

         sum = 0.0d0
         do ipoint = 1, mesh
          if (ipoint .ne. 1 .or. ipoint .ne. mesh) then
           sum = sum + drr(issh)*rr(issh,ipoint)**2*psitemp(ipoint)**2
          else
           sum = sum + 0.5d0*drr(issh)*rr(issh,ipoint)**2*psitemp(ipoint)**2
          end if
         end do

         !write (*,300) issh, sum
         !write (*,*) '  '
         !write (*,*) '  '

        end do

        !write (*,*) '  '
        !write (*,*) ' *----------------- END READPSI -------------------*'
        !write (*,*) '  '

! Deallocate Arrays
! ===========================================================================
        deallocate (psitemp)

! Format Statements
! ===========================================================================
90      format (2x, a15)
100     format (4d18.10)
200     format (' + Rc read in =', f8.4, ' abohr = ', f8.4, ' Angstroms')
300     format (' NORM (shell = ', i1, ') = ', f16.12)

        return
        end
