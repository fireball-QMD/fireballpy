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

! rcatms_excite.f90
! Program Description
! ===========================================================================
!       Once the ground state wavefunctions are determined, the excited
! states are calculated (provided iexcite = 1).  The excited states are
! calculated from the ground state potential which was obtained from the
! iteratively solving the Schroedinger equation.
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
        subroutine rcatms_excite (mesh, mesh_psirc, rc_max, rcutoff_psirc,  &
     &                            r, vc, vl, vee, uxc)
        use begin_input
        use constants
        use, intrinsic :: iso_fortran_env, only: dp => real64
	use :: utils, only: utils_itos, utils_open
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mesh

        integer, intent (in), dimension (nssh) :: mesh_psirc

        real(kind=dp), intent (in) :: rc_max

        real(kind=dp), intent (in), dimension (mesh) :: r
        real(kind=dp), intent (in), dimension (nssh) :: rcutoff_psirc
        real(kind=dp), intent (in), dimension (mesh) :: vc
        real(kind=dp), intent (in), dimension (mesh) :: vee
        real(kind=dp), intent (in), dimension (nssh, mesh) :: vl
        real(kind=dp), intent (in), dimension (mesh) :: uxc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer inum
        integer ipoint, jpoint
        integer iremainder
        integer issh
        integer jssh
        integer l, io
        integer mesh_new

        real(kind=dp) eout
        real(kind=dp) slope

        real(kind=dp), dimension (:,:), allocatable :: psi
        real(kind=dp), dimension (:), allocatable :: psi0
        real(kind=dp), dimension (:), allocatable :: v
        real(kind=dp), dimension (:), allocatable :: xx
        real(kind=dp), dimension (:), allocatable :: yy
! optimized w.f.
        real(kind=dp), dimension (:,:), allocatable :: vXe

! Allocate Arrays
! ===========================================================================
        allocate (psi(nssh,mesh))
        allocate (psi0(mesh))
        allocate (v(mesh))
        allocate (vXe(nssh,mesh))

! Procedure
! ===========================================================================
! create additional X-potential for excited orbitals to optimize 
! pseudoatomic w.f.
        vXe = 0.0_dp
        if (ioptim .eq. 1) then 
           do issh = 1,nssh
              jssh = issh + nssh
              do ipoint = 1, mesh
! (r/r0)**6.0
!                  vXe(issh,ipoint) = (r(ipoint)/r0(issh))**2.0_dp
                 if (r(ipoint) .gt. r0(jssh)) then 
! siesta formulae, PRB 64, 235111 (2001)
                    vXe(issh,ipoint) =                                       &
     &                 v0(jssh)*exp( -1.0_dp*(rcutoff(issh) - r0(jssh))       &
     &                 / (r(ipoint)-r0(jssh)) ) /(rcutoff(issh)+0.01-r(ipoint))
                 else
                    vXe(issh,ipoint) = 0.0_dp
                 endif
              end do ! enddo ipoint
           enddo ! enddo  issh
        endif ! if(ioptim)

! Loop over shell
        do issh = 1, nssh
         l = lam(issh)

! Get the final l-dependent potential
!         v = 2.0_dp*(vc + vl(issh,:)) + vee + uxc
         v = 2.0_dp*(vc + vl(issh,:)) + vee + uxc + vXe(issh,:)

! Solve Schrodinger eqn for potential v
         call psirc (mesh_psirc(issh), 1, l, rcutoff_psirc(issh), v, eout, psi0)

! Record the wavefunction answer - pnl0 to an l-dependent wavefunction
         psi(issh,1:mesh_psirc(issh)) = psi0(1:mesh_psirc(issh))
         psi(issh,mesh_psirc(issh):mesh) = 0.0_dp
        end do

! Write wavefunction onto outputfile (named FILENAME)
! First convert to angstrom units :    r   ---> r    /  (0.529177249)
!                                     r(r) ---> r(r) / (0.529177249)**1.5
! r(r) = u(r) / r
! Loop over the different states
        allocate (xx(mesh))
        allocate (yy(mesh))

        do issh = 1, nssh
         mesh_new = mesh_psirc(issh)
         l = lam(issh)

! Find (approximately) the value at r=0
         xx(2) = r(2)*abohr
         xx(3) = r(3)*abohr

         yy(2) = psi(issh,2)/(abohr1_5*r(2))
         yy(3) = psi(issh,3)/(abohr1_5*r(3))

         slope = (yy(3) - yy(2))/(xx(3) - xx(2))
         xx(1) = 0.0_dp
         yy(1) = - slope*xx(2) + yy(2)

         xx(4:mesh_new-1) = r(4:mesh_new-1)*abohr
         yy(4:mesh_new-1) = psi(issh,4:mesh_new-1)/(abohr1_5*r(4:mesh_new-1))

         xx(mesh_new) = r(mesh_new)*abohr
         yy(mesh_new) = 0.0_dp

	 if (sav(issh)) then
	   io = utils_open(trim(outpath)//trim(utils_itos(nznuc, lpad=3))//"_"//lnames(l)//".ewf", "w")
	   write (io, "('# Z, l = ',i3,2x,i2)") nznuc, l
	   write (io, "('# q = ',f14.6)") 0.0_dp
	   write (io, "('# n = ',i8)") mesh_new
	   do jpoint = 1, mesh_new
	     write (io, "(2ES16.8)") xx(jpoint), yy(jpoint)
	   end do
	   close(io)
	 end if

        end do
        deallocate (xx)
        deallocate (yy)
        deallocate (vXe)

       ! dani.JOM
       ! call vnn_excite

! Deallocate Arrays
! ===========================================================================
        return
        end
