! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! readtheory.f90
! Program Description
! ===========================================================================
!       This reads the information from the theory.input and the
! switch.input files.
!
! The two center interactions are defined as follows:
! ontop => orbitals at two different sites
! atom-atom => orbitals at same site
! interaction = 0 => do density_OSL
! interaction = 1 => do overlap
! interaction = 2 => do neutral atom/ontop function => vnnaofr
! interaction = 3 => do neutral atom/atom
! interaction = 4 => do non-local
! interaction = 5 => do xc ontop function => vxc
! interaction = 6 => do xc atom function => dvxc
! interaction = 7 => do xc double counting correction function => dexc
! interaction = 8 => do z-dipole
! interaction = 9 => do y-dipole
! interaction = 10 => do x-dipole
! interaction = 11 => do coulomb
! interaction = 12 => do extended hubbard (n1*n2*dmuxc(n1+n2), dmuxc=dmuxc/dn).
! ===========================================================================
! Code written by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readtheory (iammaster, ibcna, ibcxc, ikinetic, iswitch,   &
     &                         imuxc1c, inuxc1c, inuxc2c, isnuxc1c, isnuxc2c,&
     &                         V_intra_dip, itest, idogs, iharris, ihubbard, ispin,    &
     &                         ioomethod, ixc_opt, igauss, ngauss)
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        logical, intent (in) :: iammaster
 
! Output
! switch.input variables
        integer, intent (out) :: ibcna
        integer, intent (out) :: ibcxc
        integer, intent (out) :: ikinetic
        integer, intent (out) :: imuxc1c
        integer, intent (out) :: inuxc1c
        integer, intent (out) :: inuxc2c
        integer, intent (out) :: isnuxc1c
        integer, intent (out) :: isnuxc2c
        integer, intent (out) :: V_intra_dip
 
! two-center interactions
        integer, intent (out), dimension (0:11) :: iswitch 

! theory.input variables
        integer, intent (out) :: idogs
        integer, intent (out) :: igauss
        integer, intent (out) :: iharris
        integer, intent (out) :: ihubbard
        integer, intent (out) :: ioomethod
! jel-oslxc
        integer, intent (out) :: ixc_opt
        integer, intent (out) :: ispin
        integer, intent (out) :: itest
        integer, intent (out) :: ngauss

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer interaction
 
! Procedure
! ===========================================================================
! We now read in a theory.input file. This determines certain defaults
! for the switches which is dependent upon the the level of theory that
! is chosen.
 
! Initialize the switches. The default is to compute nothing.
        imuxc1c = 0
        ikinetic = 0
        iswitch = 0
        ibcna = 0
        ibcxc = 0
        inuxc1c = 0
        inuxc2c = 0
        isnuxc1c = 0
        isnuxc2c = 0
 
! Now we read in a theory.input file.
!        open (unit = 45, file = 'theory.input', status = 'old')
!        read (45,*) itest
!        read (45,*) iharris
!        read (45,*) idogs 
!        read (45,*) ihubbard
!        read (45,*) ispin
!! jel-oslxc
!        read (45,*) ixc_opt
!        read (45,*) ioomethod
!        read (45,*) igauss, ngauss
         itest = 1
         iharris = 0
         idogs = 0
         ihubbard = 0
         ispin = 0
         ixc_opt = 0
         ioomethod = 0
         igauss = 0
!        close (unit = 45)
        if (iammaster) then
!         write (*,*) '  '
!         write (*,*) ' Read in level of theory from theory.input. '
!         write (*,*) ' Note that the lowest level of theory is the '
!         write (*,*) ' default!! '
!         write (*,*) ' itest = 1 => read switch.input'
!         write (*,*) ' iharris = 1 => do all harris integrals. '
!         write (*,*) ' idogs = 1 => do additional dogs integrals. '
!         write (*,*) ' ihubbard = 1 => do additional hubbard integrals.'
!         write (*,*) ' ispin = 1 => do additional spin-polarization integrals.'
!         write (*,*) ' ixc_opt = 1 => do alternative option for '
!         write (*,*) '                Horsfield approach, based on SNXC'
!         write (*,*) ' ioomethod = 1 => do additional integrals for the'
!         write (*,*) '                  orbital occupancy method.'
!         write (*,*) ' igauss = 1 => compute NA and XC using gaussians '
!         write (*,*) ' ngauss is  accuracy level (1-4) of the gaussian fits. '
!         write (*,*) ' Latter options are invoked only if itest .ne. 1!'
!         write (*,*) '  '
!         write (*,*) ' Read in theory.input'
!         write (*,*) ' (itest, iharris, idogs, ihubbard) '
!         write (*,*) ' Example: 1 0 0 0 means read switch.input. '
!         write (*,*) ' Example: 0 1 0 0 means do harris integrals. '
!         write (*,*) ' Example: 0 1 1 0 means all dogs integrals. '
!         write (*,*) ' Example: 0 0 1 0 means no harris but '
!         write (*,*) '                  do additional dogs integrals. '
!         write (*,*) ' Example: 0 1 0 1 means all harris '
!         write (*,*) '                  + hubbard integrals. '
!         write (*,*) ' Example: 0 0 0 1 means ONLY hubbard integrals. '
!         write (*,*) ' Example: 0 0 0 0 means ONLY exchange integrals. '
!         write (*,*) '  '
!         write (*,*) ' (1 => compute, 0 => don''t compute)'
!         write (*,*) '  '
!         write (*,*) '     Level of Theory: '
!         write (*,*) ' itest     = ', itest
!         write (*,*) ' iharris   = ', iharris
!         write (*,*) ' idogs     = ', idogs
!         write (*,*) ' ihubbard  = ', ihubbard
!         write (*,*) ' ispin     = ', ispin
! jel-oslxc
!         write (*,*) ' ixc_opt   = ', ixc_opt
!         write (*,*) ' ioomethod = ', ioomethod
!         write (*,*) ' igauss, Gaussian Accuracy (Default=3) = ', igauss,ngauss
        end if ! end master
 
! Doing harris
        if (iharris .eq. 1) then
         imuxc1c = 1
         ikinetic = 1
         do interaction = 0, 7
          iswitch(interaction) = 1
         end do
         iswitch(11) = 1
         ibcna = 1
         ibcxc = 1
        end if
 
! Doing dogs
        if (idogs .eq. 1) then
         iswitch(0) = 1
         iswitch(2) = 1
         iswitch(3) = 1
         iswitch(5) = 1
         iswitch(6) = 1
         iswitch(7) = 1
         iswitch(8) = 1
         ibcna = 1
         ibcxc = 1
        end if
 
! Doing extended hubbard
        if (ihubbard .eq. 1) then
         inuxc1c = 1
         inuxc2c = 1
        end if
 
! Doing spin-polarization
        if (ispin .eq. 1) then
         isnuxc1c = 1
         isnuxc2c = 1
        end if

! Doing oslxc method
        if (ixc_opt .eq. 1) then
         iswitch(0) = 1
         iswitch(5) = 1
         iswitch(6) = 1
         iswitch(7) = 1
         ibcxc = 1 
        end if

! Doing orbital occupancy method
        if (ioomethod .eq. 1) then
         iswitch(5) = 1
         iswitch(6) = 1
         ibcxc = 0
        end if

! Doing gaussians 
        if (igauss .eq. 1) then
         ibcna = 0
         ibcxc = 0
        end if

! Only doing testing. Read switch.input file.
        if (itest .eq. 1) then
         if (iammaster) then
!          write (*,*) '  '
!          write (*,*) ' Switches read in from switch.input '
         end if ! end master
         !open (unit = 45, file = 'switch.input', status = 'old')
         !read (45,*) iswitch(0)
         !read (45,*) imuxc1c
         !read (45,*) ikinetic
         !do interaction = 1, 10
         ! read (45,*) iswitch(interaction)
         !end do
         !read (45,*) iswitch(11)
         !read (45,*) ibcna
         !read (45,*) ibcxc
         !read (45,*) inuxc1c
         !read (45,*) inuxc2c
         !read (45,*) isnuxc1c
         !read (45,*) isnuxc2c
         !read (45,*) V_intra_dip
         !close (unit = 45)
         iswitch(0) = 1
         imuxc1c = 1
         ikinetic = 1
         do interaction = 1, 10
          iswitch(interaction) = 1
         end do
         iswitch(11) = 1
         ibcna = 1
         ibcxc = 1
         inuxc1c = 1
         inuxc2c = 1
         isnuxc1c = 0
         isnuxc2c = 0
         V_intra_dip = 0
        end if  
 
        if (iammaster) then
!         write (*,*) '  '
!         write (*,*) ' (1 => compute, 0 => don''t compute)'
!         write (*,*) '  '
!         write (*,*) '        one-center integrals '
!         write (*,*) ' imuxc1c   = ', imuxc1c
!         write (*,*) '        two-center integrals '
!         write (*,*) ' ikinetic  = ', ikinetic
!         write (*,*) '        gaussian 3C-NA 3C-XC '
!         write (*,*) ' igauss    = ', igauss
!         write (*,*) '  '
!         write (*,*) ' interaction = 0 => do average density '
!         write (*,*) ' interaction = 1 => do overlap '
!         write (*,*) ' interaction = 2 => do neutral atom/ontop '
!         write (*,*) ' interaction = 3 => do neutral atom/atom '
!         write (*,*) ' interaction = 4 => do non-local '
!         write (*,*) ' interaction = 5 => do xc ontop function '
!         write (*,*) ' interaction = 6 => do xc atom function '
!         write (*,*) ' interaction = 7 => do xc double counting '
!         write (*,*) ' interaction = 8 => do z-dipole '
!         write (*,*) ' interaction = 9 => do y-dipole '
!         write (*,*) ' interaction = 10 => do x-dipole '
!         write (*,*) ' interaction = 11 => do coulomb '
!         write (*,*) '  '
         do interaction = 0, 11
!          write (*,*) ' interaction = ', interaction,                       &
!     &                ' iswitch = ', iswitch(interaction)
         end do
 
!         write (*,*) '  '
!         write (*,*) '      three-center integrals '
!         write (*,*) ' ibcxc     = ', ibcxc
!         write (*,*) ' ibcna     = ', ibcna

!         write (*,*) '  '
!         write (*,*) '  extended-Hubbard integrals '
!         write (*,*) ' inuxc1c   = ', inuxc1c
!         write (*,*) ' inuxc2c   = ', inuxc2c

!         write (*,*) '  '
!         write (*,*) '  extended-Hubbard (spin) integrals '
!         write (*,*) ' isnuxc1c   = ', isnuxc1c
!         write (*,*) ' isnuxc2c   = ', isnuxc2c
        end if ! end master
 
! Format Statements
! ===========================================================================
        return
      end subroutine readtheory
