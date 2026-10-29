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
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use constants
	use :: utils, only: utils_open  
        use begin_input, only: outpath
        implicit none

! Argument Declaration and Description
! ===========================================================================
!Input
        integer, intent(in) :: nssh
        character(len=9), intent(in), dimension (nssh) :: filein
        integer, intent(in) :: mesh_max
! Output
        integer, intent(out), dimension (nssh) :: npoints
        real(kind=dp), intent(out), dimension (nssh) :: drr
        real(kind=dp), intent(out), dimension (nssh, mesh_max) :: psi
        real(kind=dp), intent(out), dimension (nssh, mesh_max) :: rr
        real(kind=dp), intent(out), dimension (nssh) :: rrc
        real(kind=dp), intent(out), dimension (nssh) :: xnocc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer issh, io

! Procedure
! ===========================================================================

! Open the input file
       do issh = 1, nssh
	 io = utils_open(trim(outpath)//trim(filein(issh)), "r", skip=1)
	 read (io, "(6x,f14.6)") xnocc(issh)
	 read (io, "(6x,i8)") npoints(issh)
	 psi(issh, :) = 0.0_dp
	 do ipoint = 1, npoints(issh)
	   read (io, "(2ES16.8)") rr(issh, ipoint), psi(issh, ipoint)
	 end do
	 rrc(issh) = rr(issh, npoints(issh))*invabohr
	 drr(issh) = rr(issh, npoints(issh))/real(npoints(issh) - 1, kind=dp)
       end do

! Format Statements
! ===========================================================================
        end
