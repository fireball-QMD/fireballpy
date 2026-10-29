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
        use, intrinsic :: iso_fortran_env, only: dp => real64
	use :: utils, only: utils_open, utils_itos
	implicit none

! Argument Declaration and Description
! ===========================================================================
        real(kind=dp), parameter :: etotatom = 0.0_dp

! Local Parameters and Data Declaration
! ===========================================================================
! grid sizes
        integer, parameter :: nrr = 2049
        integer, parameter :: nrrp = 2047
        integer, parameter :: mesh_rcatms = 2049
! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer ir
        integer irp
        integer issh
        integer itenth
        integer l, io
        integer mesh

        integer, dimension (:), allocatable :: npoints_vnn

        real(kind=dp), external :: derf0
        real(kind=dp) dr
        real(kind=dp) drp1, drp2
        real(kind=dp) fact
        real(kind=dp), external :: psiofr
        real(kind=dp) r_vnn
        real(kind=dp) r21
        real(kind=dp) r22
        real(kind=dp) rc_max_vnn
        real(kind=dp) rp1, rp2
        real(kind=dp) vcore
        real(kind=dp), external :: vshort
!dani.JOM
       real(kind=dp), dimension (mesh_rcatms) :: r
       real(kind=dp), dimension (mesh_rcatms) :: vc
       real(kind=dp), dimension (:,:), allocatable :: vl
       integer ioptionpp
       real(kind=dp) exmix
      
!        real(kind=dp) xnz

        real(kind=dp), dimension (:), allocatable :: drr
        real(kind=dp), dimension (:,:), allocatable :: psi_vnn
        real(kind=dp), dimension (:,:), allocatable :: rr
        real(kind=dp), dimension (:), allocatable :: rrc
        real(kind=dp), dimension (:), allocatable :: term1
        real(kind=dp), dimension (:), allocatable :: term2
        real(kind=dp), dimension (nrr) :: vneut
        real(kind=dp), dimension (nssh) :: vrho
        real(kind=dp), dimension (:,:), allocatable :: vvrho
        real(kind=dp), dimension (:), allocatable :: xnocc
        real(kind=dp), dimension (nrr) :: xrr
	character(len=9), dimension(:), allocatable :: filename_wf

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
	allocate (filename_wf(nssh))
	  
! Procedure
! ===========================================================================
! Initialize some constants
        mesh = mesh_rcatms

!-------dani.JOM 
        call pp (0, mesh_rcatms, r, vc, vl, ioptionpp, exmix)
!------------------	

! Read in the wavefunctions
	do issh = 1, nssh
	  l = lam(issh)
	  filename_wf(issh) = trim(utils_itos(nznuc, lpad=3))//"_"//lnames(l)//".wf "
	end do
        allocate (rr(nssh,mesh_rcatms))
        allocate (psi_vnn(nssh,mesh_rcatms))
        call readpsi (nssh, filename_wf, mesh, npoints_vnn, rrc,        &
     &                drr, xnocc, rr, psi_vnn)

	rc_max_vnn = -1.0_dp
        do issh = 1, nssh
         rc_max_vnn = max(rc_max_vnn,rrc(issh))
        end do
        rc_max_vnn = rc_max_vnn*abohr


! *****************************************************************************
! first set up vna matrix
        dr = rc_max_vnn/real(nrr - 1)
        r_vnn = - dr

! integration over r
        do ir = 1, nrr
         r_vnn = r_vnn + dr
         term1 = 0.0_dp
         term2 = 0.0_dp

         drp1 = r_vnn/real(nrrp - 1)
         drp2 = (rc_max_vnn - r_vnn)/real(nrrp - 1)
         rp1 = 0.0_dp
         rp2 = r_vnn

! integration over r'
         do irp = 1, nrrp
          fact = 4.0_dp
          if (mod(irp,2) .ne. 0) fact = 2.0_dp
          if (irp .eq. 1 .or. irp .eq. nrrp) fact = 1.0_dp

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
          term1(issh) = (drp1/3.0_dp)*term1(issh)
          term2(issh) = (drp2/3.0_dp)*term2(issh)
          vrho(issh) = term2(issh)
          if (r_vnn .gt. 1.0e-5_dp) vrho(issh) = vrho(issh) + term1(issh)/r_vnn
         end do

! Evaluate the core potential
! First the erf-piece local potential. This is the dp-range piece.
         if (r_vnn .ne. 0.0_dp) then
!          vcore = - (xnz/r_vnn)*derf0(alpha*r_vnn)
          vcore = - (Zval/r_vnn)*derf0(alpha*r_vnn)
         else
!          vcore = - xnz*2.0_dp*alpha/sqrt(pi)
          vcore = - Zval*2.0_dp*alpha/sqrt(pi)
         end if
         vneut(ir) = vcore + vshort(r_vnn)/eq2

! V_na is sum of all the pieces
         do issh = 1, nssh
          vneut(ir) = vneut(ir) + vrho(issh)*xnocc(issh)
         end do

         vvrho(:,ir) = vrho
         xrr(ir) = r_vnn
        end do

! Write out the information to the files

        do issh = 0, nssh
	  if (issh == 0) then
	    io = utils_open(trim(outpath)//trim(utils_itos(nznuc, lpad=3))//".na0", "w")
	    write (io, "('# Z = ',i3)") nznuc
	  else
	    l = lam(issh)
	    io = utils_open(trim(outpath)//trim(utils_itos(nznuc, lpad=3))//"_"//lnames(l)//".na", "w")
	    write (io, "('# Z, l = ',i3,2x,i2)") nznuc, l
	  end if
	  write (io, "('# n = ',i8)") nrr
	  if (issh == 0) then
	    do ipoint = 1, nrr
	      write (io, "(2ES16.8)") xrr(ipoint), vneut(ipoint)
	    end do
	  else
	    do ipoint = 1, nrr
	      write (io, "(2ES16.8)") xrr(ipoint), vvrho(issh,ipoint)
	    end do
	  end if
	  close(io)
        end do

! If the option to calculate excited states is chosen then call this routine
!dani.JOM 
         if (nexcite .eq. 1 .or. nexcite .eq. 2) call vnn_excite()


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
	deallocate (filename_wf)
        end
