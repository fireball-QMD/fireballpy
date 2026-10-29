! copyright info:
!
!                             @Copyright 1998
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
!                      Kirk VanOpdorp
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

 
! kinetic.f90
! Program Description
! ===========================================================================
! Calculates kinetic energy matrix elements
!
! Calls psiofr and iofile2c
!
! Specifically, this routine will create kinetic energy matrix elements
! for two atomic species designated by in1,in2, which gets passed in
! the call list. appropriately named output files are  created
!
! ===========================================================================
! Original code by Wolfgang Windl with new additions by Gary Adams.
 
! Code rewritten by:
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
        subroutine kinetic (itype1, itype2, atom1, atom2, what1, what2,      &
     &                      nzx1, nzx2, rcutoff1, rcutoff2, nssh, lssh,      &
     &                      index_max, nleft, lleft, mleft, nright, lright,  &
     &                      mright, signature, iammaster)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use dimensions
        use quadrature
	use :: constants, only: pi, fname_len
	use :: utils, only: utils_open, utils_itos, utils_log, UTILS_LOG_INFO, UTILS_LOG_ERROR
	use :: indices, only: indices_twocenter_set, INDICES_TWOCENTER, INDICES_NAMES_LEN
        use :: wavefunctions, only: wf_atoms
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: index_max
        integer, intent (in) :: itype1
        integer, intent (in) :: itype2
        integer, intent (in) :: nzx1
        integer, intent (in) :: nzx2
 
! information relevant to the storage of the non-zero matrix elements
! for two and three center matrix elements.
!
! terms 1 .. index_max2c(in1,in2): all two-center matrix elements
!
! terms index_max2c(in1,in2) + 1 .. index_max3c(in1,in2):
!       all additional terms appearing in the three center case
!
        integer, intent (in), dimension (inter_max) :: nleft
        integer, intent (in), dimension (inter_max) :: lleft
        integer, intent (in), dimension (inter_max) :: mleft
        integer, intent (in), dimension (inter_max) :: nright
        integer, intent (in), dimension (inter_max) :: lright
        integer, intent (in), dimension (inter_max) :: mright
 
! l quantum number
        integer, intent (in), dimension (nspec_max, nsh_max) :: lssh

! number of shells
        integer, intent (in), dimension (nspec_max) :: nssh 
 
        character (len=2), intent (in) :: atom1, atom2
        character (len=70), intent (in) :: signature
        character (len=70), intent (in) :: what1, what2
 
        logical iammaster
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer kmax
        parameter (kmax = 6)        ! maximum number of surviving bessels
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer id
        integer imu
        integer index
        integer iounit
        integer iq
        integer ir
        integer isplit
        integer issh
        integer ktype
        integer l1, l2
        integer m1, m2
        integer n1, n2
        integer nqq
        integer nqtop
        integer nrr
        integer nsplit
 
        integer inj (2)
 
        real(kind=dp) ang_integral
        real(kind=dp) d
        real(kind=dp) dmax
        real(kind=dp) dq
        real(kind=dp) dr
        real(kind=dp) eV
        real(kind=dp) factor
        real(kind=dp), external :: jl2
        real(kind=dp) q
        real(kind=dp) qmax
        real(kind=dp) r
        real(kind=dp) rcutoff1
        real(kind=dp) rcutoff2
        real(kind=dp) rrmax
        real(kind=dp) rmin
        real(kind=dp) xnqtop
        real(kind=dp) xtra
 
        real(kind=dp) angular (0:kmax, inter_max)
        real(kind=dp) answer (0:kmax)
 
! The 5 in the esplit is to split q into 5 ranges for comparison purposes.
! This is similarly done for xnormq and qsplit below.
        real(kind=dp) esplit (5)
        real(kind=dp) qsplit (5)
 
        real(kind=dp) rmax (2)
        real(kind=dp) rjj (nsh_max, 2, 2*nqke + 1)
        real(kind=dp) sumj (nsh_max, 2)
        real(kind=dp) tkinetic (inter_max, nddke)
        real(kind=dp) xnormq (2, 0:nsh_max, 5)
 
        character (len = 1) ang
        character (len = fname_len) filename
	character (len = INDICES_NAMES_LEN), allocatable :: names(:)
	integer, allocatable :: ls1(:), ls2(:), s1(:), s2(:), nl1(:), nl2(:), nm1(:), nm2(:)
	integer nindex
 
        logical skip
 
! Procedure
! ===========================================================================
 
! ang is the Angstrom symbol (A)
        ang = char(197)
 
	! Open necessary file
	filename = "kinetic."//utils_itos(nzx1, lpad=2)//"."//utils_itos(nzx2, lpad=2)//".dat"
	! Con 'exists' (que no existe como especificador) skip quedaba sin
	! definir; y con la negacion la logica salia al reves: se saltaba el
	! calculo cuando el fichero NO estaba, asi que no se escribia nunca.
	inquire (file=trim(filename), exist=skip)
        if (skip) return
	iounit = utils_open(filename, "w")
 
        factor = 7.62_dp/2.0_dp
        qmax = sqrt(ecutke/factor)
 
! **************************************************************************
! This section finds the numbers and limits necessary for integration.
! Note: we use different Rc's for different atoms now;
! It is better if the atoms are really different (e.g. O & Na)
        rmax(1) = rcutoff1
        rmax(2) = rcutoff2
        nqq = 2*nqke + 1
        nrr = 2*nrke + 1
        dq = qmax/real(nqq - 1)
 
! This section caculates and stores R0(q), R1(q), R2(q), R3(q).
! RL(q)=int {jL(qr) * RL(r) * r**2} dr
! Remember that psi_nlm (r) = Y_lm(Omega_r) R_nl(r), and
!               psi_nlm (q) = 4 pi (-i)^l Y_lm(Omega_q) R_nl(q)
 
        q = - dq
        rmin = 0.0_dp
        inj(1) = itype1
        inj(2) = itype2
 
 
! ***************************************************************************
!
! E S T A B L I S H   A N G U L A R   I N T E G R A T I O N   F A C T O R S
! ***************************************************************************
        do index = 1, index_max
         l1 = lleft(index)
         m1 = mleft(index)
         l2 = lright(index)
         m2 = mright(index)
 
! The variable angular(k) is the integral of Ylm*Pk[cos(theta)]*Yl'm'.
! We only consider terms up to k = 6. This is for the following reason
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! For these values of l1 and l2 k is bounded from above by 6.
! See my notes or check this out for yourself if you don't believe me!
! kmax is dimensioned in parameters above.
         call Pintegral (l1, m1, l2, m2, kmax, answer)
 
         angular(0,index) = answer(0)
         angular(1:kmax,index) = 0.0_dp
         do imu = 1, kmax
          if (mod(l1-l2-imu,2) .eq. 0) then
           angular(imu,index) = answer(imu)*(2*imu+1)*((-1)**((l1-l2-imu)/2))
          end if
         end do
        end do
 
! ***************************************************************************
!
! F O U R I E R   T R A N S F O R M   O F   R A D I A L   C O M P O N E N T
! ***************************************************************************
        rjj = 0.0_dp
        do iq = 1, nqq
         q = q + dq
 
! Use Simpson's rule
! Now do integral over r for this fixed q value. The r integral goes from 0 to
! rcutoff(i). Loop over both atoms:
         do iatom = 1, 2
          dr = (rmax(iatom) - rmin)/real(nrr - 1)
          r = - dr
          do ir = 1, nrr
           r = r + dr
           factor = 4.0_dp/3.0_dp
           if (mod(ir,2) .eq. 1) factor = 2.0_dp/3.0_dp
           if (ir .eq. 1 .or. ir .eq. nrr) factor = 1.0_dp/3.0_dp
           do issh = 1, nssh(inj(iatom))
            rjj(issh,iatom,iq) = rjj(issh,iatom,iq)                          &
     &       + factor*jl2(lssh(inj(iatom),issh),q*r)                         &
     &               *wf_atoms(inj(iatom))%get_psi(issh,r)*r*r*dr
           end do
          end do
         end do
        end do
 
! Test the normalization
!    - we first normalize in q space.
!    - in q space we get a factor (4*pi)/((2pi)**3)=2/pi
!    - we then normalize in r space
        xtra = 2.0_dp/pi
        nsplit = 5
        do isplit = 1, nsplit
         xnqtop = real(nqq)*sqrt(real(isplit)/real(nsplit))
         nqtop = int(xnqtop)
         sumj = 0.0_dp
         q = - dq
         do iq = 1, nqtop
          q = q + dq
          factor = (4.0_dp/3.0_dp)
          if (mod(iq,2) .eq. 1) factor = (2.0_dp/3.0_dp)
          if (iq .eq. 1 .or. iq .eq. nqq) factor = (1.0_dp/3.0_dp)
          do iatom = 1, 2
           do issh = 1, nssh(inj(iatom))
            sumj(issh,iatom) = sumj(issh,iatom)                              &
     &       + factor*q*q*dq*rjj(issh,iatom,iq)**2*xtra
           end do
          end do
         end do
 
         qsplit(isplit) = q
         esplit(isplit) = ecutke*((q/qmax)**2)
         do iatom = 1, 2
          do issh = 1, nssh(inj(iatom))
           xnormq(iatom,issh,isplit) = sumj(issh,iatom)
          end do
         end do
        end do
 
! ***************************************************************************
! E N D   F O U R I E R   T R A N S F O R M
!
! ***************************************************************************
 
! ***************************************************************************
!
! C A L C U L A T E   K I N E T I C   E N E R G Y
! ***************************************************************************
        rrmax = rcutoff1 + rcutoff2
        dr = rrmax/real(nddke - 1)
        d = - dr
 
! eV= ((h/2pi)**2)/m/pi (in eV times angstroms squared)
        eV = 2.4255_dp
 
! Now calculate the kinetic energy for all distances.
! ***************************************************************************
        do id = 1, nddke
         d = d + dr
 
! Initialize
         tkinetic(1:index_max,id) = 0.0_dp
 
! Now the integral over q loop.
         q = - dq
         do iq = 1, nqq
          q = q + dq
 
! Simpson rule for integration.
          factor = (4.0_dp/3.0_dp)*eV
          if (mod(iq,2) .eq. 1) factor = (2.0_dp/3.0_dp)*eV
          if (iq .eq. 1 .or. iq .eq. nqq) factor = (1.0_dp/3.0_dp)*eV
 
          do index = 1, index_max
           n1 = nleft(index)
           n2 = nright(index)
 
! We only consider terms up to k = 6. This is for the following reason.
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! The maximum value of k is 6!
           ang_integral = 0.0_dp
           do ktype = 0, 6
            ang_integral = ang_integral + angular(ktype,index)*jl2(ktype,q*d)
           end do
           tkinetic(index,id) = tkinetic(index,id)                           &
     &      + factor*rjj(n1,1,iq)*rjj(n2,2,iq)*ang_integral*q**4*dq
          end do
         end do
        end do
 
! ***************************************************************************
!
! W R I T E   T H E   O U T P U T
! ***************************************************************************
! CRP: use new format :)
      call utils_log(UTILS_LOG_INFO, "  - "//trim(filename))
      dmax = rcutoff1 + rcutoff2
      write (iounit, "('# Z1, Z2 = ',i3,', ',i3)") nzx1, nzx2
      write (iounit, "('# Rc1, Rc2 = ',f14.6,', ',f14.6)") rcutoff1, rcutoff2
      if (index_max == 0) then
        write (iounit, "('# n = ',i8)") 0
        write (iounit, "('# z')")
        close (iounit)
        return
      end if
      ls1 = [(wf_atoms(itype1)%get_angular_momentum(id), id = 1, wf_atoms(itype1)%get_nshells())]
      ls2 = [(wf_atoms(itype2)%get_angular_momentum(id), id = 1, wf_atoms(itype2)%get_nshells())]
      call indices_twocenter_set(INDICES_TWOCENTER, ls1, ls2, nindex, s1, s2, nl1, nl2, nm1, nm2, names=names)
      if (nindex /= index_max) then
        call utils_log(UTILS_LOG_ERROR, "kinetic: mk_index and indices_twocenter_set disagree")
        error stop 1
      end if
      write (iounit, "('# n = ',i8)") nddke
      write (iounit, "('# z',1000(4x,a))") (trim(names(index)), index = 1, index_max)
      do id = 1, nddke
        write (iounit, "(f14.6,1000ES16.8)") (id - 1)*dr, tkinetic(:, id)
      end do
      close (iounit)
 
        close (iounit)
        return
        end
