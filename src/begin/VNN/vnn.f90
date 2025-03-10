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
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! vnn.f90
! Program Description
! ===========================================================================
!       This routine computes the local (non)-neutral atom potential vna(r)
! or vnn(s,p,d, or f) of an atom and puts it in a data file; e.g. 014.na,
! 014_s.na, etc.
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
        subroutine vnn ()
        use begin_input
        use constants
        use pp_storage
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
        real(kind=long), parameter :: etotatom = 0.0d0

! Local Parameters and Data Declaration
! ===========================================================================
! grid sizes
        integer, parameter :: nrr = 399
        integer, parameter :: nrrp = 401
        integer, parameter :: mesh_rcatms = 1001
! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer ir
        integer irp
        integer issh
        integer itenth
        integer l
        integer mesh

        integer, dimension (:), allocatable :: npoints_vnn

        real(kind=long), external :: derf0
        real(kind=long) dr
        real(kind=long) drp1, drp2
        real(kind=long) fact
        real(kind=long), external :: psiofr
        real(kind=long) r_vnn
        real(kind=long) r21
        real(kind=long) r22
        real(kind=long) rc_max_vnn
        real(kind=long) rp1, rp2
        real(kind=long) vcore
        real(kind=long), external :: vshort
!dani.JOM
       real(kind=long), dimension (mesh_rcatms) :: r
       real(kind=long), dimension (mesh_rcatms) :: vc
       real(kind=long), dimension (:,:), allocatable :: vl
       integer ioptionpp
       real(kind=long) exmix
      
!        real(kind=long) xnz

        real(kind=long), dimension (:), allocatable :: drr
        real(kind=long), dimension (:,:), allocatable :: psi_vnn
        real(kind=long), dimension (:,:), allocatable :: rr
        real(kind=long), dimension (:), allocatable :: rrc
        real(kind=long), dimension (:), allocatable :: term1
        real(kind=long), dimension (:), allocatable :: term2
        real(kind=long), dimension (nrr) :: vneut
        real(kind=long), dimension (nssh) :: vrho
        real(kind=long), dimension (:,:), allocatable :: vvrho
        real(kind=long), dimension (:), allocatable :: xnocc
        real(kind=long), dimension (nrr) :: xrr

! Allocate Arrays
! ===========================================================================
        allocate (drr(nssh))
        allocate (npoints_vnn(nssh))
        allocate (rrc(nssh))
        allocate (term1(nssh))
        allocate (term2(nssh))
        allocate (vvrho(nssh,nrr))
        allocate (xnocc(nssh))
        allocate (vl(nssh,mesh_rcatms))

! Procedure
! ===========================================================================
! Initialize some constants
        !write (*,*) '  '
        !write (*,*) '  '
        !write (*,*) '  *------------------------------------------------*'
        !write (*,*) '  |            Starting subroutine vnn             |'
        !write (*,*) '  *------------------------------------------------*'
        !write (*,*) '  '
        !write (*,*) '  Atomic number = ', nznuc
        !write (*,*) '  '
!-------dani.JOM 
 call pp (0, mesh_rcatms, r, vc, vl, ioptionpp, exmix)
!------------------	


	
        !write (*,*) ' This program will read in data files containing the '
        !write (*,*) ' wavefunctions. These files were created by the rcatms.f '
        !write (*,*) ' subroutine. The file names should be according to our '
        !write (*,*) ' convention - '
        !write (*,*) ' 1. Begin with the chemical number (Z) (e.g. 14 for Si)'
        !write (*,*) ' which should be followed by the value for rc.'
        !write (*,*) ' 2. The extension should be .wf0, .wf1, .wf2, depending '
        !write (*,*) ' on whether the file is s, p or d. '
        !write (*,*) '  '
        !write (*,*) ' An example is 014_500.wf2 which is a wavefunction file '
        !write (*,*) ' for the d-orbital of Silicon and the wavefunction has '
        !write (*,*) ' a maximum cutoff of 5.00 Bohr. '
        !write (*,*) '  '

        !write (*,*) ' We read the wavefunctions from: '
        do issh = 1, nssh
         !write (*,100) filename_wf(issh)
        end do

        !write (*,*) '  '
        !write (*,*) ' Furthermore, the program will write out a data file'
        !write (*,*) ' containing the (non)-neutral atom (local) potential.'
        !write (*,*) ' Its name should according to our convention - '
        !write (*,*) ' 1. Begin with the chemical number (Z) (e.g. 14 for Si)'
        !write (*,*) ' which should be followed by the value for rc.'
        !write (*,*) ' 2. Signify to which shell the potential belongs. '
        !write (*,*) ' 3. The extension should be _s.na, _p.na, _d.na, or _f.na '
        !write (*,*) '  '
        !write (*,*) ' An example is 014_500_s.na which is the s-shell '
        !write (*,*) ' non-neutral atom potential file for Silicon with a '
        !write (*,*) ' maximum cutoff of 5.00 Bohr.'
        !write (*,*) '  '
        !write (*,*) ' We write the non-neutral atom potential to: '
        do issh = 0, nssh
         !write (*,100) filename_na(issh)
        end do
        !write (*,*) '  '

! Read in the wavefunctions
! First get the maximum mesh size from all of the wavefuntion files
        call get_wfmesh_size (nssh, filename_wf, mesh)
        allocate (rr(nssh,mesh))
        allocate (psi_vnn(nssh,mesh))
        call readpsi (nssh, filename_wf, mesh, npoints_vnn, rrc,        &
     &                drr, xnocc, rr, psi_vnn)

        rc_max_vnn = -1.0d0
        do issh = 1, nssh
         rc_max_vnn = max(rc_max_vnn,rrc(issh))
        end do
        rc_max_vnn = rc_max_vnn*0.529177249d0

!        xnz = 0.0d0
!        do issh = 1, nssh
!         xnz = xnz + xnocc(issh)
!        end do


! *****************************************************************************
! first set up vna matrix
        !write (*,*) '  '
        !write (*,*) ' Setting up V_na matrix. '
        !write (*,*) '   '
        !write (*,*) ' Write out the potential at 10 points:'
        !write (*,*) ' ====================================='
        dr = rc_max_vnn/real(nrr - 1)
        r_vnn = - dr

! integration over r
        do ir = 1, nrr
         r_vnn = r_vnn + dr
         term1 = 0.0d0
         term2 = 0.0d0

         drp1 = r_vnn/real(nrrp - 1)
         drp2 = (rc_max_vnn - r_vnn)/real(nrrp - 1)
         rp1 = 0.0d0
         rp2 = r_vnn

! integration over r'
         do irp = 1, nrrp
          fact = 4.0d0
          if (mod(irp,2) .ne. 0) fact = 2.0d0
          if (irp .eq. 1 .or. irp .eq. nrrp) fact = 1.0d0

          do issh = 1, nssh
           r21 = psiofr(issh, rp1, nssh, mesh, npoints_vnn, drr,        &
     &                  rrc, rr, psi_vnn)**2
           r22 = psiofr(issh, rp2, nssh, mesh, npoints_vnn, drr,        &
     &                  rrc, rr, psi_vnn)**2
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
          if (r_vnn .gt. 1.0d-5) vrho(issh) = vrho(issh) + term1(issh)/r_vnn
         end do

! Evaluate the core potential
! First the erf-piece local potential. This is the long-range piece.
         if (r_vnn .ne. 0.0d0) then
!          vcore = - (xnz/r_vnn)*derf0(alpha*r_vnn)
          vcore = - (Zval/r_vnn)*derf0(alpha*r_vnn)
         else
!          vcore = - xnz*2.0d0*alpha/sqrt(pi)
          vcore = - Zval*2.0d0*alpha/sqrt(pi)
         end if
         vneut(ir) = vcore + vshort(r_vnn)/eq2

! V_na is sum of all the pieces
         do issh = 1, nssh
          vneut(ir) = vneut(ir) + vrho(issh)*xnocc(issh)
         end do

         vvrho(:,ir) = vrho
         xrr(ir) = r_vnn

! Write out 10 values
         itenth = int(float(nrr)/10.0d0)
         if (mod(ir,itenth) .eq. 1) then
          !write (*,*) '  '
          !write (*,200) ir, xrr(ir), vneut(ir)
          do issh = 1, nssh
           !write (*,201) issh - 1, vvrho(issh,ir)
          end do
         end if
        end do

! Write out the information to the files
        !write (*,*) ' '
        !write (*,*) ' Finished setting up V_na matrix. '
        !write (*,*) ' Now write values on the following data files: '

        do issh = 0, nssh
         !write (*,100) filename_na(issh)
         open (unit = 66, file = trim(outpath)//trim(filename_na(issh)),  status = 'unknown')
         write (66,101) filename_na(issh)
         write (66,*) nznuc

         if (issh .eq. 0) write (66,*) rc_max_vnn/0.529177249d0
         if (issh .ne. 0) write (66,*) rcutoff(issh)
         write (66,*) nrr
         if (issh .eq. 0) write (66,*) etotatom

         do ipoint = 1, nrr
          if (issh .eq. 0) write (66,300) xrr(ipoint), vneut(ipoint)
          if (issh .ne. 0) write (66,300) xrr(ipoint), vvrho(issh,ipoint)
         end do
         close (unit = 66)


! For plotting purposes
!         l = -99
!dani.         if (issh .eq. 0) open (unit = 48, file = 'NA_plot', status = 'unknown')
!dani.         if (issh .ge. 1) l = lam(issh)
!dani.         if (l .eq. 0) open (unit = 48, file = 'NA_splot', status = 'unknown')
!dani.         if (l .eq. 1) open (unit = 48, file = 'NA_pplot', status = 'unknown')
!dani.         if (l .eq. 2) open (unit = 48, file = 'NA_dplot', status = 'unknown')
!dani.         if (l .eq. 3) open (unit = 48, file = 'NA_fplot', status = 'unknown')
!dani.         do ipoint = 1, nrr
!dani.          if (issh .eq. 0) write (48,301) xrr(ipoint), vneut(ipoint)
!dani.          if (issh .ne. 0) write (48,301) xrr(ipoint), vvrho(issh,ipoint)
!dani.         end do
!dani.         close (unit = 48)
        end do

! If the option to calculate excited states is chosen then call this routine
!dani.JOM 
         if (nexcite .eq. 1 .or. nexcite .eq. 2) then
         call vnn_excite
         end if


! Deallocate Arrays
! ===========================================================================
        deallocate (drr)
        deallocate (npoints_vnn)
        deallocate (psi_vnn)
        deallocate (rr)
        deallocate (rrc)
        deallocate (term1)
        deallocate (term2)
        deallocate (vvrho)
        deallocate (xnocc)
        deallocate (vl)

! Format Statements
! ===========================================================================
100     format (2x, a11)
101     format (a11)
200     format (2x, ' ir = ', i4, 2x, ' r = ', f9.5, 2x, ' V_na = ', d16.5)
201     format (22x, ' Vvrho(issh = ', i1, ') = ', d16.5)
300     format (2d24.16)
301     format (2x, f6.3, 2x, f12.6)

        return
        end
