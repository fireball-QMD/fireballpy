! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
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

! get_wfmesh_size.f90
! Program Description
! ===========================================================================
!       This program read in the first few lines of all the input files
! to determine the maximum mesh size for the wavefuntion.

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
        subroutine get_wfmesh_size (nssh, filein, mesh)
        use begin_input, only: outpath
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: nssh

        character(len=11), intent(in), dimension (nssh) :: filein

! Output
        integer, intent(out) :: mesh

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer issh
        integer meshin

! Procedure
! ===========================================================================
        mesh = -1001

! Open all of the input files
        do issh = 1, nssh
         open (12, file = trim(outpath)//trim(filein(issh)), status = 'old')
         read (12,*)
         read (12,*)
         read (12,*) meshin
         mesh = max(mesh,meshin)
         close (12)
        end do

! Format Statements
! ===========================================================================

        return
        end
