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

! get_nlmesh_size.f90
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
        subroutine get_nlmesh_size (ion, mesh)
        use begin_input, only: ppfile, ppionfile, outpath
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: ion

! Output
        integer, intent(out) :: mesh

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ioption
        integer ipoint
        integer issh
        integer nssh
        integer meshin

        integer, dimension (:), allocatable :: lssh

        real(kind=long) alpha
        real(kind=long) r_nl
        real(kind=long) r_short
        real(kind=long) v_nl
        real(kind=long) v_short
! jel-PP
        real(kind=long) Z_val
        real(kind=long) rc_PP

! Procedure
! ===========================================================================
        mesh = -1001

        if (ion .eq. 0) then
          open (unit = 15, file = trim(outpath)//trim(ppfile), status = 'old')
        elseif (ion .eq. 1) then
          open (unit = 15, file = trim(outpath)//trim(ppionfile), status = 'old')
        else
          stop
        end if

! Read in everything up until the point where the
! non-local pseudopotential is located.
! There are several message lines in the pseudopotential file
        do ipoint = 1, 14
         read (15,*)
        end do

! Read the exchange-correlation option
        read (15,*) ioption

! Read the number of shells
        read (15,*) nssh
        allocate (lssh(nssh))
        read (15,*) (lssh(issh), issh = 1, nssh)
        deallocate (lssh)

! jel-PP
! Read in Z_val
        read (15,*) Z_val

! Read in alpha
        read (15,*) alpha
!jel-PP
! Read in rc_PP
        read (15,*) rc_PP

! Read in the short-range local part  - this is not needed for the crtor
        read (15,*) meshin
        do ipoint = 1, meshin
         read (15,*) r_short, v_short
        end do

        do issh = 1, nssh
         read (15,100) meshin
         do ipoint = 1, meshin
          read (15,*) r_nl, v_nl
         end do
         mesh = max(mesh,meshin)
        end do

        close (unit = 15)

! Format Statements
! ===========================================================================
100     format (12x, i5)

        return
        end
