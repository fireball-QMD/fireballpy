! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
!
! fireball-qmd is a free (GPLv3) open project.

!
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

 
! readcreate.f
! Program Description
! ===========================================================================
!       This reads the information from the create.input file
!       This also reads the *.input files for each atom listed
!
! Program Declaration
! ===========================================================================
        subroutine readcreate (nspec, iammaster, iammpi, atom, what, &
     &                         nssh, lssh, nzx, rcutoff, rcutoffa, &
     &                         rcutoffa_max, xmass, ppfile, napot, &
     &                         wavefxn)
        use precision, only: wp
        implicit none
 
        include '../parameters.inc'
        include '../exchange.inc'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer nspec
 
        logical iammaster
        logical iammpi
 
! Output
        integer lssh (nspec_max, nsh_max)      ! l quantum number
        integer nssh (nspec_max)               ! number of shells
        integer nzx (nspec_max)
 
        real(kind=wp) rcutoff (nspec_max, nsh_max) ! cutoff radius in bohr
        real(kind=wp) rcutoffa (nspec_max, nsh_max)! cutoff radius in angstroms
        real(kind=wp) rcutoffa_max (nspec_max)     ! cutoff radius in angstroms
        real(kind=wp) xmass (nspec_max)
 
        character(len=2) atom (nspec_max)
        character(len=25) napot (nspec_max, 0:nsh_max)
        character(len=25) ppfile (nspec_max)
        character(len=25) wavefxn (nspec_max, nsh_max)
        character(len=70) what (nspec_max)
 
        integer nzx_max
 
! Local Parameters and Data Declaration
! ===========================================================================
        character(len=2) periodic (103)
        data periodic  / 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', &
     &       'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
     &       'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', &
     &       'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', &
     &       'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', &
     &       'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', &
     &       'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', &
     &       'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', &
     &       'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', &
     &       'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', &
     &       'Es', 'Fm', 'Md', 'No', 'Lw' /
 
        real(kind=wp) abohr
        parameter (abohr = 0.529177249d0)
 
! Local Variable Declaration and Description
! ===========================================================================
        integer icontinue
        integer issh
        integer ispec
 
        real(kind=wp) add
 
        character(len=2) atomcheck
 
        character(len=15) inputfile (nspec_max)
 
        logical read_input
 
! Procedure
! ===========================================================================
! We now read in a create.input file. This determines the number of atoms
! and the types of atoms.
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Carefully edit create.input for your situation.'
         write (*,*) '  '
         write (*,*) ' We now read create.input! '
         write (*,*) '  '
        end if ! end master
 
        open (unit = 44, file = 'create.input', status = 'old')
        read (44,*) nspec
        if (nspec .gt. nspec_max) then
         write (*,*) '  '
         write (*,*) ' nspec > nspec_max '
         write (*,*) ' Redimension creator! '
         stop 'error in readcreate'
        end if
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Use the same convention for numbering as '
         write (*,*) ' the dynamics: Z1 < Z2 < Z3 < ... < Z(nspec) '
         write (*,*) ' where Z is atomic number for each of the atoms. '
         write (*,*) '  '
        end if ! end master
 
        nzx_max = 0
        do ispec = 1, nspec
         read (44,105) inputfile(ispec)
         inquire (file = inputfile(ispec), exist = read_input)
         if (read_input) then
          open (unit = 45, file = inputfile(ispec), status = 'old')
         else
          write (*,*) ' The following input file does not exist! '
          write (*,104) inputfile(ispec)
          stop 'error in readcreate'
         end if
 
         read (45,101) atom(ispec)
         read (45,*) nzx(ispec)
 
         nzx_max = max(nzx(ispec),nzx_max)
         if (nzx(ispec) .lt. nzx_max) then
          write (*,*) ' ispec = ', ispec
          write (*,*) ' Z(ispec) .lt. Z(ispec-1) '
          write (*,*) ' nzx(ispec) = ', nzx(ispec)
          write (*,*) ' nzx(ispec-1) = ', nzx(ispec-1)
          stop ' Must stop bad order.  Z1 < Z2 < Z3 ... violated'
         end if
 
! Check whether you put in the correct nz for that atom.
! Go through the periodic table and check.
         atomcheck = periodic(nzx(ispec))
         if (iammaster) then
          write (*,200) ispec, nzx(ispec), atom(ispec), atomcheck
         end if ! end master
         if (atom(ispec) .ne. atomcheck) stop ' wrong nz(nuc) for atom!!'
 
         read (45,*) xmass(ispec)
 
! Read filename for the pseudopotential
         read (45,102) ppfile(ispec)
 
! Read filename for the neutral atom potential
         read (45,102) napot(ispec,0)
 
! Now read stuff in related to the orbital information
         read (45,*) nssh(ispec)
         nsshxc(ispec) = nssh(ispec)
 
! Loop over the number of shells:
! Read the l quantum number (0, 1, 2, and 3 => s, p, d, and f), the occupation
! number, cutoff radius (in bohr), wavefunction, and neutral atom potential
! for each shell.
         do issh = 1, nssh(ispec)
          read (45,*) lssh(ispec,issh)
          lsshxc(ispec,issh) = lssh(ispec,issh)
          read (45,*) xnocc(issh,ispec)
          read (45,*) rcutoff(ispec,issh)
          read (45,102) wavefxn(ispec,issh)
          read (45,102) napot(ispec,issh)
         end do
 
! Check that xns + xnp .lt. nzx
         add = 0.0d0
         do issh = 1, nssh(ispec)
          add = add + xnocc(issh,ispec)
         end do
         if (add .gt. nzx(ispec)) then
          write (*,*) ' xnocc = ', &
     &                 (xnocc(issh,ispec), issh = 1, nssh(ispec))
          write (*,*) ' Nuclear Z = ', nzx(ispec)
          write (*,*) ' Huh? Does this make sense? '
          write (*,*) ' Sorry I must stop. How can the number of '
          write (*,*) ' electrons be larger than nuclear Z? '
          stop 'error in readcreate'
         end if
 
! Read information for the xc interactions - iderorb is the shell where
! the charge will be changed in the +-Q derivative stuff, and dqorb is
! the charge amount changed.
! HAO modified at April 8, 2005
!        read (45,*) iderorb(ispec)
!         read (45,*) dqorb(ispec)
! jel-dq
!         read (45,*) (dqint(issh,ispec),issh = 1, nssh(ispec))
! end jel-dq
          iderorb(ispec) = nssh(ispec)
            dqorb(ispec) = 0.5d0
            if (nssh(ispec) .eq. 1) dqorb(ispec) = 0.25d0
            do issh = 1, nssh(ispec)
              dqint(issh,ispec) = dqorb(ispec)/nssh(ispec)
            end do             
! End of HAO
        end do
        close (unit = 44)
        close (unit = 45)
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Done reading in create.input'
         write (*,*) ' Please check it and if it looks OK, proceed. '
         write (*,*) '  '
         write (*,*) ' Number of species = ', nspec
         write (*,*) ' '
 
! Check for continuation.
!         write (*,*) ' Insert 0 to continue if the above is OK. '
!         if (.not. iammpi) read (*,*) icontinue
        end if ! end master
!        if (iammpi) icontinue = 0
        icontinue = 0
        if (icontinue .ne. 0) stop ' OK I''ll stop. '
 
! Set up the useful Angstrom array rcutoffa, and what.
        do ispec = 1, nspec
         do issh = 1, nssh(ispec)
          rcutoffa(ispec,issh) = rcutoff(ispec,issh)*abohr
         end do
         write (what(ispec),400) atom(ispec), nzx(ispec), &
     &    (rcutoffa(ispec,issh), issh = 1, nssh(ispec))
        end do
 
! First find the largest rc.
        do ispec = 1, nspec
         rcutoffa_max(ispec) = -1.0d0
         do issh = 1, nssh(ispec)
          rcutoffa_max(ispec) = &
     &     max(rcutoffa_max(ispec),rcutoffa(ispec,issh))
         end do
         if (iammaster) then
          write (*,*) '  '
          write (*,*) ' For ispec = ', ispec
          write (*,*) ' Largest rcutoffa_max (Angstrom)= ',  &
     &  rcutoffa_max(ispec)
          write (*,*) '  '
         end if ! end master
        end do
 
! Format Statements
! ===========================================================================
101     format (a2)
102     format (a25)
104     format (2x, a25)
105     format (a15)
200     format (' Species = ',i2,' Nuclear Z = ',i3,' atomname = ',a2, &
     &          ' In periodic table = ', a2)
400     format (2x, a2, '(',i2,')', ': Rc''s (A) = ', 8f6.3)
 
        return
        end
