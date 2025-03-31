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
! last change Daniel G. Trabada, http://fireball.ftmc.uam.es/moodle/course/view.php?id=4
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! pp.f90
! Program Description
! ===========================================================================
! Calculate pseudopotential in r space.  The potential output is in Hartree
! units (program uses atomic units).

! Output: vc(i), vl(l,i), where l = angular momentum, and i = 1, mesh
! vc = core potential
! vl = non-local potentials
! thus the total potential for s states is vc(i)+vl(0,i),
!      the total potential for p states is vc(i)+vl(1,i),
!      the total potential for d states is vc(i)+vl(2,i), etc.
! total potential for other than considered states is just vc.
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
        subroutine pp (ion, mesh, r, vc, vl, ioptionpp, exmix)
        use begin_input
        use constants
        use pp_storage
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ion
        integer, intent (in) :: mesh
!       integer, intent (in) :: nssh

!       real(kind=long), intent (in) :: xnz

        real(kind=long), intent (in), dimension (mesh) :: r

! Output
        integer, intent (out) :: ioptionpp

        real(kind=long), intent (out) :: exmix

        real(kind=long), intent (out), dimension (mesh) :: vc
        real(kind=long), intent (out), dimension (nssh, mesh) :: vl

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: norder = 3

! Local Variable Declaration and Description
! ===========================================================================
        integer ileft
        integer imid
        integer ipoint
        integer iprod
        integer iright
        integer issh
        integer isum
        integer nsshPP
        integer nlmesh
        integer ppmesh
        integer :: imu = 1
        integer jssh
        
        integer, dimension (:), allocatable :: lsshPP
        integer, dimension (:), allocatable :: mapPP

        real(kind=long), external :: derf0
        real(kind=long) prod
        real(kind=long) rc_PP
        real(kind=long) rpoint
        real(kind=long) vcore
        real(kind=long), external :: vshort


! Allocate Arrays
! ===========================================================================
!        allocate (lsshPP(nssh))
        allocate (mapPP(nssh)) 

! Procedure
! ===========================================================================
        !write (*,*) '  '
        !write (*,*) ' Welcome to Fireballs2005 pseudopotential. '

! ***************************************************************************
! Find out the mesh size for the non-local piece of the pseudopotential.
        if (ion .eq. 0) then
          call get_nlmesh_size (ion, nlmesh)
          open (unit = 88, file = trim(outpath)//trim(ppfile), status = 'old')
        elseif (ion .eq. 1) then
          call get_nlmesh_size (ion, nlmesh)
          open (unit = 88, file = trim(outpath)//trim(ppionfile), status = 'old')
        else
          close (unit = 88)
          stop
        end if

! Read in the psuedopotential
        !write (*,*) ' open file ', ppfilein
        !open (unit = 88, file = trim(outpath)//trim(ppfilein), status = 'old')

! There are 14 message lines in the pseudopotential file
        do ipoint = 1, 7
         read (88,*)
        end do
        read (88,100) charge
        do ipoint = 9, 14
         read (88,*)
        end do

! Read which exchange-correlation approximation was used.
        read (88,*) ioptionpp, exmix

! Read the number of shells
! nssh is the number of L dependent pseudopotentials.
! If we have S and P, then nssh = 2.
        read (88,*) nsshPP
! If we must allocate lsshPP(nsshPP) read in the pp file not lsshPP(nssh) 
        allocate (lsshPP(nsshPP)) 
        !write (*,*) ' nssh read in ', nsshPP
        !write (*,*) ' nssh = ', nssh
        if (nsshPP .le. 0 .or. nsshPP .gt. 4)                                &
     &   stop ' Check nssh in pseudopotential file '
        read (88,*) (lsshPP(issh), issh = 1, nsshPP)
        !write (*,*) 'lsshPP = ', lsshPP(1:nsshPP)

! Read Zval
        read (88,*) Zval
        if (nzval_pp .gt. 0.0d0) Zval = nzval_pp
        write (*,*) ' Zval : ', Zval

! Read in alpha for longranged local part => -Z*e**2*erf(alpha*r)/r
        read (88,*) alpha
        !write (*,*) ' alpha read in: ', alpha

! Read Rcutoff of PP 
        read (88,*) rc_PP
        !write (*,*) ' rc_PP read in: ', rc_PP

! Read in the short-range local part  - this is not needed for the crtor
        read (88,*) ppmesh
        if (allocated(r_short) .and. ppmesh .gt. npoints_short)              &
         deallocate (r_short, v_short)
        npoints_short = ppmesh
        if (.not. allocated(r_short)) allocate (r_short(ppmesh))
        if (.not. allocated(v_short)) allocate (v_short(ppmesh))

        do ipoint = 1, ppmesh
         read (88,*) r_short(ipoint), v_short(ipoint)
        end do
        rrc_short = r_short(ppmesh)
        drr_short = rrc_short/real(ppmesh - 1)


! ***************************************************************************
! Now read in the non-local pseudopotential piece.
! First allocate the necessary arrays.
        allocate (drr_nl(nsshPP))
        allocate (npoints_nl(nsshPP))
        allocate (r_nl(nsshPP,nlmesh))
        allocate (rrc_nl(nsshPP))
        allocate (v_nl(nsshPP,nlmesh))

        do issh = 1, nsshPP

! We put the pseuodpotentials in as S, P, D etc.
         read (88,200) npoints_nl(issh)

! Allocate the non-local arrays
         do ipoint = 1, npoints_nl(issh)
          read (88,*) r_nl(issh,ipoint), v_nl(issh,ipoint)
         end do
         rrc_nl(issh) = r_nl(issh,npoints_nl(issh))
         drr_nl(issh) = rrc_nl(issh)/real(npoints_nl(issh) - 1)
        end do

        close (unit = 88)
        !write (*,*) ' Finished reading data from pseduopotential file.'


! ***************************************************************************
! Now compute the potential in real space.
! Units notice!!! The pseudopotential we read in is in eV, Angstrom units.
! This atomic program uses Rydberg-abohr units. Kevin is changing this
! program, but for NOW, we must convert our eV-Angstrom units to Rdy-abohr
! units.
        do ipoint = 1, mesh
         rpoint = r(ipoint)

! First the erf-piece local potential. This is the long-range piece.
         if (rpoint .gt. 1.0d-4) then
          vcore = - (Zval/rpoint)*derf0(alpha*abohr*rpoint)

         else
          vcore = - Zval*2.0d0*alpha*abohr/sqrt(pi)
         end if
         rpoint = rpoint*abohr
         vc(ipoint) = vcore + vshort(rpoint)/hartree
        end do

! reset non-local part 
        vl = 0.0d0
! Now the non-local angular momentum dependent potential.
        do issh = 1, nssh
        mapPP(issh)=0
! map non-local angular momentum into vl array index
         do jssh = 1,nsshPP
          if (lsshPP(jssh) .eq. lam(issh)) then
           imu = jssh
           mapPP(issh)=1
         !  !write (*,*) 'imu = ',imu,lam(issh) 
          endif
         end do ! do jssh
          if(mapPP(issh).eq.1) then
           !write (*,*) 'imu = ', imu,'map in issh =',issh
          else
           !write (*,*) 'imu = ', imu,'not map in issh =',issh
          end if
! Change to abohr

         if (mapPP(issh).eq.1) then

          do ipoint = 1, mesh
           rpoint = r(ipoint)*abohr
           vl(issh,ipoint) = 0.0d0

! The short stuff is zero at long range.
           if (rpoint .lt. rrc_nl(imu)) then
            imid = int(rpoint/drr_nl(imu)) + 1

! Find starting and ending points for the interpolation
            if (mod(norder + 1, 2) .eq. 0) then
             ileft = imid - ((norder-1)/2)
             iright = imid + ((norder+1)/2)
            else
             ileft = imid - (norder/2)
             iright = imid + (norder/2)
            end if

            if (ileft .lt. 1) then
             ileft = 1
             iright = norder + 1
            else if (iright .gt. npoints_nl(imu)) then
             ileft = npoints_nl(imu) - (norder+1)
             iright = npoints_nl(imu)
            end if

! Now interpolate
            do isum = ileft, iright
             prod = 1.0d0
             do iprod = ileft, iright
              if (iprod .ne. isum) then
               prod = prod*(rpoint - r_nl(imu,iprod))/                  &
     &                    (r_nl(imu,isum) - r_nl(imu,iprod))
              end if
             end do
             vl(issh,ipoint) = vl(issh,ipoint) + v_nl(imu,isum)*prod

            end do ! do isum
           end if ! if (rpoint)
          end do ! do ipoint
         end if ! if mapPP
        end do ! do issh
! change units
        vl = vl/hartree

        !write (*,*) '  '
        !write (*,*) ' Done creating pseudopotetnial. '
        !write (*,*) '  '

! Deallocate Arrays
! ===========================================================================
        deallocate (drr_nl)
        deallocate (lsshPP)
        deallocate (npoints_nl)
        deallocate (r_nl)
        deallocate (rrc_nl)
        deallocate (v_nl)

! Format Statements
! ===========================================================================
100     format (31x, f6.2)
200     format (12x, i5)

        return
        end
