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
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mesh

        integer, intent (in), dimension (nssh) :: mesh_psirc

        real(kind=long), intent (in) :: rc_max

        real(kind=long), intent (in), dimension (mesh) :: r
        real(kind=long), intent (in), dimension (nssh) :: rcutoff_psirc
        real(kind=long), intent (in), dimension (mesh) :: vc
        real(kind=long), intent (in), dimension (mesh) :: vee
        real(kind=long), intent (in), dimension (nssh, mesh) :: vl
        real(kind=long), intent (in), dimension (mesh) :: uxc

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer inum
        integer ipoint, jpoint
        integer iremainder
        integer issh
        integer jssh
        integer l
        integer mesh_new

        real(kind=long) eout
        real(kind=long) slope

        real(kind=long), dimension (:,:), allocatable :: psi
        real(kind=long), dimension (:), allocatable :: psi0
        real(kind=long), dimension (:), allocatable :: v
        real(kind=long), dimension (:), allocatable :: xx
        real(kind=long), dimension (:), allocatable :: yy
! optimized w.f.
        real(kind=long), dimension (:,:), allocatable :: vXe

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
        vXe = 0.0d0
        if (ioptim .eq. 1) then 
           do issh = 1,nssh
              jssh = issh + nssh
              do ipoint = 1, mesh
! (r/r0)**6.0
!                  vXe(issh,ipoint) = (r(ipoint)/r0(issh))**2.0d0
                 if (r(ipoint) .gt. r0(jssh)) then 
! siesta formulae, PRB 64, 235111 (2001)
                    vXe(issh,ipoint) =                                       &
     &                 v0(jssh)*exp( -1.0d0*(rcutoff(issh) - r0(jssh))       &
     &                 / (r(ipoint)-r0(jssh)) ) /(rcutoff(issh)+0.01-r(ipoint))
                 else
                    vXe(issh,ipoint) = 0.0d0
                 endif
              end do ! enddo ipoint
           enddo ! enddo  issh
! write vX        
           do ipoint = 1, mesh
!              write (100,*) r(ipoint),(vXe(issh,ipoint),issh=1,nssh)
           enddo
        endif ! if(ioptim)

! Loop over shell
        do issh = 1, nssh
         l = lam(issh)

! Get the final l-dependent potential
!         v = 2.0d0*(vc + vl(issh,:)) + vee + uxc
         v = 2.0d0*(vc + vl(issh,:)) + vee + uxc + vXe(issh,:)

! Solve Schrodinger eqn for potential v
         call psirc (mesh_psirc(issh), 1, l, rcutoff_psirc(issh), v, eout, psi0)

! Record the wavefunction answer - pnl0 to an l-dependent wavefunction
         psi(issh,1:mesh_psirc(issh)) = psi0(1:mesh_psirc(issh))
         psi(issh,mesh_psirc(issh):mesh) = 0.0d0
        end do

! Write out the excited state wavefunction.
        !write (*,*) ' This program writes out the following data files: '
        !write (*,*) '  '
        !write (*,*) '  '
        !write (*,*) ' s.ewf [p.ewf...]  u_s(r) [u_p...] r = 100 points '
        !write (*,*) '  '

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
         if(sav(issh)) open (unit = 14, file = trim(outpath)//trim(filename_ewf(issh)), status = 'unknown')
         if(sav(issh)) write(14,101) filename_ewf(issh)
         if(sav(issh)) write(14,102) nznuc, atomname
         if(sav(issh)) write(14,103) mesh_new
         if(sav(issh)) write(14,104) rcutoff_psirc(issh), rc_max, 0.0d0
!                               ^
!                               |---------- zero occupation
         if(sav(issh)) write(14,105) lam(issh)

! Find (approximately) the value at r=0
         xx(2) = r(2)*abohr
         xx(3) = r(3)*abohr

         yy(2) = psi(issh,2)/(abohr15*r(2))
         yy(3) = psi(issh,3)/(abohr15*r(3))

         slope = (yy(3) - yy(2))/(xx(3) - xx(2))
         xx(1) = 0.0d0
         yy(1) = - slope*xx(2) + yy(2)

         xx(4:mesh_new-1) = r(4:mesh_new-1)*abohr
         yy(4:mesh_new-1) = psi(issh,4:mesh_new-1)/(abohr15*r(4:mesh_new-1))

         xx(mesh_new) = r(mesh_new)*abohr
         yy(mesh_new) = 0.0d0

         inum = idint(real(mesh_new)/4.d0)
         iremainder = mesh_new - (inum*4)
         do jpoint = 1, mesh_new - iremainder,4
          if(sav(issh)) write(14,400) yy(jpoint), yy(jpoint+1), yy(jpoint+2), yy(jpoint+3)
         end do

         if (iremainder .eq. 1) then
          if(sav(issh)) write(14,400) yy(mesh_new)
         else if (iremainder .eq. 2) then
          if(sav(issh)) write(14,400) yy(mesh_new-1), yy(mesh_new)
         else if (iremainder .eq. 3) then
          if(sav(issh)) write(14,400) yy(mesh_new-2), yy(mesh_new-1), yy(mesh_new)
         end if
         if(sav(issh)) close (unit = 14)

! Write out the wavefunctions for plotting purposes.
!dani.         if (l .eq. 0) open (unit = 17, file = 'sstate1', status = 'unknown')
!dani.         if (l .eq. 1) open (unit = 17, file = 'pstate1', status = 'unknown')
!dani.         if (l .eq. 2) open (unit = 17, file = 'dstate1', status = 'unknown')
!dani.         if (l .eq. 3) open (unit = 17, file = 'fstate1', status = 'unknown')
!dani.         do ipoint = 1, mesh_new
!dani.          write (17,401) r(ipoint), yy(ipoint)
!dani.         end do
!dani.         close (unit = 17)
        end do
        deallocate (xx)
        deallocate (yy)
        deallocate (vXe)

       ! dani.JOM
       ! call vnn_excite

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
101     format (2x, a12)
102     format (2x, i5, 2x, a12)
103     format (2x, i8)
104     format (2x, 2(2x,f7.4), 2x, f7.2)
105     format (2x, i3)
400     format (4d18.10)
401     format (2x, f6.3, 2x, f12.6)

        return
        end
