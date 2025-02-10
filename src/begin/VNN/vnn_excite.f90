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

! vnn_excite.f90
! Program Description
! ===========================================================================
!       This routine computes the local (non)-neutral atom potential vna(r)
! or vnn(s,p,d, or f) of an atom and puts it in a data file; e.g. 014.na,
! 014_s.na, etc.  This is a routine just for the excited states.
!
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
        subroutine vnn_excite
        use begin_input
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================
! grid sizes
        integer, parameter :: nrr = 399
        integer, parameter :: nrrp = 401

! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer ir
        integer irp
        integer issh
        integer itenth
        integer l
        integer mesh

        integer, dimension (:), allocatable :: npoints

        real(kind=long), external :: derf0
        real(kind=long) dr
        real(kind=long) drp1, drp2
        real(kind=long) fact
        real(kind=long), external :: psiofr
        real(kind=long) r
        real(kind=long) r21
        real(kind=long) r22
        real(kind=long) rc_max
        real(kind=long) rp1, rp2

        real(kind=long), dimension (:), allocatable :: drr
        real(kind=long), dimension (:,:), allocatable :: psi
        real(kind=long), dimension (:,:), allocatable :: rr
        real(kind=long), dimension (:), allocatable :: rrc
        real(kind=long), dimension (:), allocatable :: term1
        real(kind=long), dimension (:), allocatable :: term2
        real(kind=long), dimension (nssh) :: vrho
        real(kind=long), dimension (:,:), allocatable :: vvrho
        real(kind=long), dimension (:), allocatable :: xnocc
        real(kind=long), dimension (nrr) :: xrr

! Allocate Arrays
! ===========================================================================
        allocate (drr(nssh))
        allocate (npoints(nssh))
        allocate (rrc(nssh))
        allocate (term1(nssh))
        allocate (term2(nssh))
        allocate (vvrho(nssh,nrr))
        allocate (xnocc(nssh))

! Procedure
! ===========================================================================
! Initialize some constants
        !write (*,*) '  '
        !write (*,*) '  '
        !write (*,*) '  *------------------------------------------------*'
        !write (*,*) '  |            Starting subroutine vnn             |'
        !write (*,*) '  |               for excited states               |'
        !write (*,*) '  *------------------------------------------------*'
        !write (*,*) '  '
        !write (*,*) '  Atomic number = ', nznuc
        !write (*,*) '  '

        !write (*,*) ' We read the wavefunctions from: '
        do issh = 1, nssh
         !write (*,100) filename_ewf(issh)
        end do

! Read in the wavefunctions
! First get the maximum mesh size from all of the wavefuntion files
        call get_ewfmesh_size (nssh, filename_ewf, mesh)
        allocate (rr(nssh,mesh))
        allocate (psi(nssh,mesh))

        call readpsi_excite (nssh, filename_ewf, mesh, npoints, rrc, drr,    &
     &                       xnocc, rr, psi)

        rc_max = -1.0d0
        do issh = 1, nssh
         rc_max = max(rc_max,rrc(issh))
        end do
        rc_max = rc_max*0.529177249d0


! *****************************************************************************
! first set up vna matrix
        !write (*,*) '  '
        !write (*,*) ' Setting up V_na matrix for excited states. '
        !write (*,*) '   '
        !write (*,*) ' Write out the potential at 10 points:'
        !write (*,*) ' ====================================='
        dr = rc_max/real(nrr - 1)
        r = - dr

! integration over r
        do ir = 1, nrr
         r = r + dr
         term1 = 0.0d0
         term2 = 0.0d0

         drp1 = r/real(nrrp - 1)
         drp2 = (rc_max - r)/real(nrrp - 1)
         rp1 = 0.0d0
         rp2 = r

! integration over r'
         do irp = 1, nrrp
          fact = 4.0d0
          if (mod(irp,2) .ne. 0) fact = 2.0d0
          if (irp .eq. 1 .or. irp .eq. nrrp) fact = 1.0d0

          do issh = 1, nssh
           r21 = psiofr(issh, rp1, nssh, mesh, npoints, drr, rrc, rr, psi)**2
           r22 = psiofr(issh, rp2, nssh, mesh, npoints, drr, rrc, rr, psi)**2
           term1(issh) = term1(issh) + fact*rp1*rp1*r21
           term2(issh) = term2(issh) + fact*rp2*r22
          end do

          rp1 = rp1 + drp1
          rp2 = rp2 + drp2
         end do

         do issh = 1, nssh
          term1(issh) = (drp1/3.0d0)*term1(issh)
          term2(issh) = (drp2/3.0d0)*term2(issh)
          vrho(issh) = term2(issh)
          if (r .gt. 1.0d-5) vrho(issh) = vrho(issh) + term1(issh)/r
         end do

         vvrho(:,ir) = vrho
         xrr(ir) = r

! Write out 10 values
         itenth = int(float(nrr)/10.0d0)
         if (mod(ir,itenth) .eq. 1) then
          !write (*,*) '  '
          !write (*,200) ir, xrr(ir)
          do issh = 1, nssh
           !write (*,201) issh - 1, vvrho(issh,ir)
          end do
         end if
        end do

! Write out the information to the files
        !write (*,*) ' '
        !write (*,*) ' Finished setting up V_na matrix. '
        !write (*,*) ' Now write values on the following data files: '
        do issh = 1, nssh
         !write (*,100) filename_ena(issh)
         open (unit = 66, file = trim(outpath)//trim(filename_ena(issh)),  status = 'unknown')
         write (66,101) filename_ena(issh)
         write (66,*) nznuc
         write (66,*) rcutoff(issh)
         write (66,*) nrr

         do ipoint = 1, nrr
          write (66,300) xrr(ipoint), vvrho(issh,ipoint)
         end do
         close (unit = 66)

! For plotting purposes
      !   l = lam(issh)
      !   if (l .eq. 0) open (unit = 47, file = 'NAe_splot', status = 'unknown')
      !   if (l .eq. 1) open (unit = 47, file = 'NAe_pplot', status = 'unknown')
      !   if (l .eq. 2) open (unit = 47, file = 'NAe_dplot', status = 'unknown')
      !   if (l .eq. 3) open (unit = 47, file = 'NAe_fplot', status = 'unknown')
      !   do ipoint = 1, nrr
      !    write (47,301) xrr(ipoint), vvrho(issh,ipoint)
      !   end do
      !   close (unit = 47)
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (drr)
        deallocate (npoints)
        deallocate (psi)
        deallocate (rr)
        deallocate (rrc)
        deallocate (term1)
        deallocate (term2)
        deallocate (vvrho)
        deallocate (xnocc)

! Format Statements
! ===========================================================================
100     format (2x, a12)
101     format (a12)
200     format (2x, ' ir = ', i4, 2x, ' r = ', f9.5)
201     format (22x, ' Vvrho(issh = ', i1, ') = ', d16.5)
300     format (2d24.16)
301     format (2x, f6.3, 2x, f12.6)

        return
        end
