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
        use dimensions
        use quadrature
        use precision, only: wp
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
 
        real(kind=wp) ang_integral
        real(kind=wp) d
        real(kind=wp) dmax
        real(kind=wp) dq
        real(kind=wp) dr
        real(kind=wp) eV
        real(kind=wp) factor
        real(kind=wp), external :: jl2
        real(kind=wp) pi
        real(kind=wp), external :: psiofr
        real(kind=wp) q
        real(kind=wp) qmax
        real(kind=wp) r
        real(kind=wp) rcutoff1
        real(kind=wp) rcutoff2
        real(kind=wp) rrmax
        real(kind=wp) rmin
        real(kind=wp) xnqtop
        real(kind=wp) xtra
 
        real(kind=wp) angular (0:kmax, inter_max)
        real(kind=wp) answer (0:kmax)
 
! The 5 in the esplit is to split q into 5 ranges for comparison purposes.
! This is similarly done for xnormq and qsplit below.
        real(kind=wp) esplit (5)
        real(kind=wp) qsplit (5)
 
        real(kind=wp) rmax (2)
        real(kind=wp) rjj (nsh_max, 2, 2*nqke + 1)
        real(kind=wp) sumj (nsh_max, 2)
        real(kind=wp) tkinetic (inter_max, nddke)
        real(kind=wp) xnormq (2, 0:nsh_max, 5)
 
        character (len = 1) ang
        character (len = 40) filename
        character (len = 20) root
 
        logical skip
 
! Procedure
! ===========================================================================
! Initialize pi
        pi = 3.141592653589793238462643D0
 
! ang is the Angstrom symbol (A)
        ang = char(197)
 
! Open necessary file
        iounit = ((itype2 - 1)*nspec_max + itype1 - 1) + 36
        root = 'coutput/kinetic'
        call iofile2c (root, 'dat', nzx1, nzx2, iounit, filename, skip)
        if (skip) return
 
        factor = 7.62d0/2.0d0
        qmax = sqrt(ecutke/factor)
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' *------------------------------------------* '
         write (*,*) ' | Welcome to the kinetic-energy subroutine | '
         write (*,*) ' *------------------------------------------* '
         write (*,*) '  '
         write (*,*) ' Who and when: ', signature
         write (*,*) '  '
 
         write (*,100) itype1, itype2
 
! .............................this gets qmax...............................
         write (*,*) '  '
         write (*,*) '  '
         write (*,101) ecutke
         write (*,*) ' (10 keV could be a good choice) '
         write (*,*) '  '
         write (*,102) qmax, ang
 
! ..............................this gets nq................................
         write (*,*) '  '
         write (*,103) nqke
         write (*,*) ' (a good choice would be 400) '
         write (*,*) '  '
 
! ..............................this gets nr................................
         write (*,104) nrke
         write (*,*) ' (a good choice would be 240) '
         write (*,*) '  '
         write (*,*) '  '
        end if  !end master
 
! **************************************************************************
! This section finds the numbers and limits necessary for integration.
! Note: we use different Rc's for different atoms now;
! It is better if the atoms are really different (e.g. O & Na)
        rmax(1) = rcutoff1
        rmax(2) = rcutoff2
        nqq = 2*nqke + 1
        nrr = 2*nrke + 1
        dq = qmax/real(nqq - 1)
        if (iammaster) then
         write (*,105) rmax(1), ang, rmax(2), ang
         write (*,*) '  '
         write (*,*) ' One moment, we calculate the Fourier transforms of the '
         write (*,*) ' radial WF. '
        end if  !end master
 
! This section caculates and stores R0(q), R1(q), R2(q), R3(q).
! RL(q)=int {jL(qr) * RL(r) * r**2} dr
! Remember that psi_nlm (r) = Y_lm(Omega_r) R_nl(r), and
!               psi_nlm (q) = 4 pi (-i)^l Y_lm(Omega_q) R_nl(q)
 
        q = - dq
        rmin = 0.0d0
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
 
        if(iammaster)then
         write (*,601) index, nleft(index), lleft(index), mleft(index),      &
     &                        nright(index), lright(index), mright(index)
        end if  !end master
 
! The variable angular(k) is the integral of Ylm*Pk[cos(theta)]*Yl'm'.
! We only consider terms up to k = 6. This is for the following reason
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! For these values of l1 and l2 k is bounded from above by 6.
! See my notes or check this out for yourself if you don't believe me!
! kmax is dimensioned in parameters above.
         call Pintegral (l1, m1, l2, m2, kmax, answer)
 
         angular(0,index) = answer(0)
         angular(1:kmax,index) = 0.0d0
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
        rjj = 0.0d0
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
           factor = 4.0d0/3.0d0
           if (mod(ir,2) .eq. 1) factor = 2.0d0/3.0d0
           if (ir .eq. 1 .or. ir .eq. nrr) factor = 1.0d0/3.0d0
           do issh = 1, nssh(inj(iatom))
            rjj(issh,iatom,iq) = rjj(issh,iatom,iq)                          &
     &       + factor*jl2(lssh(inj(iatom),issh),q*r)                         &
     &               *psiofr(inj(iatom),issh,r)*r*r*dr
           end do
          end do
         end do
        end do
 
! Test the normalization
!    - we first normalize in q space.
!    - in q space we get a factor (4*pi)/((2pi)**3)=2/pi
!    - we then normalize in r space
        xtra = 2.0d0/pi
        nsplit = 5
        do isplit = 1, nsplit
         xnqtop = real(nqq)*sqrt(real(isplit)/real(nsplit))
         nqtop = int(xnqtop)
         sumj = 0.0d0
         q = - dq
         do iq = 1, nqtop
          q = q + dq
          factor = (4.0d0/3.0d0)
          if (mod(iq,2) .eq. 1) factor = (2.0d0/3.0d0)
          if (iq .eq. 1 .or. iq .eq. nqq) factor = (1.0d0/3.0d0)
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
 
        if (iammaster) then
         write (*,*)
         write (*,*) ' **** Test normalization in q space **** '
         do iatom = 1, 2
          write (*,*) '  '
          write (*,200) iatom
          do issh = 1, nssh(inj(iatom))
           write (*,*) '  '
           write (*,201) issh
           write (*,*) '  '
           write (*,202) issh
           write (*,*) '  '
           do isplit = 1, 5
            write (*,203) qsplit(isplit), esplit(isplit),                   &
     &                    xnormq(iatom,issh,isplit)
           end do
          end do
         end do
         write (*,*) '  '
        end if  !end master
 
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
 
        if (iammaster) then
         write (*,*) ' Now calculate the kinetic matrix element for '
         write (*,*) ' each d from 0 to rmax.'
         write (*,*) '  '
         write (*,300) nddke
         write (*,*) ' (107 would have been a good choice, e.g.) '
         write (*,*) '  '
         write (*,301) rrmax, ang
         write (*,*) '  '
         write (*,302) dr, ang
         write (*,*) '  '
         write (*,*) '  '
 
! The following section calculates the kinetic energy matrix elements
! for all the possible values of d
! Some preliminaries. Set up simpson rule factors and eV.
! Do a convergence test for d = 0.0
 
         write (*,*) ' Convergence tests of Kinetic energy integral '
         write (*,*) ' at d = 0.0 <Phi(atom1)| T | Phi(atom2>. '
        end if  !end master
 
! eV= ((h/2pi)**2)/m/pi (in eV times angstroms squared)
        eV = 2.4255d0
 
! Now calculate the kinetic energy for all distances.
! ***************************************************************************
        do id = 1, nddke
         d = d + dr
 
! Initialize
         tkinetic(1:index_max,id) = 0.0d0
 
! Now the integral over q loop.
         q = - dq
         do iq = 1, nqq
          q = q + dq
 
! Simpson rule for integration.
          factor = (4.0d0/3.0d0)*eV
          if (mod(iq,2) .eq. 1) factor = (2.0d0/3.0d0)*eV
          if (iq .eq. 1 .or. iq .eq. nqq) factor = (1.0d0/3.0d0)*eV
 
          do index = 1, index_max
           n1 = nleft(index)
           n2 = nright(index)
 
! We only consider terms up to k = 6. This is for the following reason.
! The selection rules are such that |l1-k|<=l2<=l1+k.  Therefore, if we
! are only doing up to f-orbitals here, l1 = 0, 1, 2, 3 and l2 = 0, 1, 2, 3.
! The maximum value of k is 6!
           ang_integral = 0.0d0
           do ktype = 0, 6
            ang_integral = ang_integral + angular(ktype,index)*jl2(ktype,q*d)
           end do
           tkinetic(index,id) = tkinetic(index,id)                           &
     &      + factor*rjj(n1,1,iq)*rjj(n2,2,iq)*ang_integral*q**4*dq
          end do
         end do
         if (id .eq. 1) then
          if (iammaster) then
           write (*,*) '  '
           do isplit = 1, 5
            write (*,400) qsplit(isplit), esplit(isplit),                    &
     &                    (tkinetic(index,id), index = 1, index_max)
           end do
          end if  !end master
         end if
        end do
 
! ***************************************************************************
!
! W R I T E   T H E   O U T P U T
! ***************************************************************************
! Set up the header strings:
        write (iounit,500)
        write (iounit,*) ' Kinetic matrix elements '
        write (iounit,*) ' created by: '
        write (iounit,501) signature
        write (iounit,502) what1
        write (iounit,502) what2
        write (iounit,500)
 
        dmax = rcutoff1 + rcutoff2
        write (iounit,503) dmax, nddke, nqke, nrke
        write (iounit,500)
 
        write (iounit,504) nzx1, rcutoff1, atom1
        write (iounit,504) nzx2, rcutoff2, atom2
 
        write (iounit,*) dmax, nddke
        write (iounit,*) '  '
 
        do id = 1, nddke
         write (iounit,505) (tkinetic(index,id),index = 1, index_max)
        end do
 
        write (*,506) filename
 
! Format Statements
! ===========================================================================
100     format (' Doing kinetic for type ', i2, ' and type ', i2)
101     format (2x, ' The energy cutoff (for q sum) is ', f9.2, ' eV.')
102     format (2x, ' The k-space cutoff qmax = ', f8.4, ' 1/', a1, '.')
103     format (2x, ' The number of q points is 2*nq+1; nq = ', i5)
104     format (2x, ' The number of r points is 2*nr+1; nr = ', i5)
105     format (2x, ' In the Fourier transform we use rmax1 = ',f7.4,1x,     &
     &          a1,' and rmax2 = ', f7.4, 1x, a1, '.')
200     format (2x, ' Atomic species #', i2, ':')
201     format (2x, ' Convergence: Shell = ', i2)
202     format (8x, 'q', 12x, 'E', 5x, 'Int[psi (l = ', i1,',q)**2]')
203     format (3x, f9.3, 1x, f11.2, 1x, f12.7)
300     format (2x, 'You chose', i4, ' grid points between 0 and rmax.')
301     format (2x, ' Rc(type1) + Rc(type2) = ', f7.4, 1x, a1, '.')
302     format (2x, ' The point separation is ', f7.4, 1x, a1, '.')
400     format (2x, ' q(max) = ', f8.2, 3x, ' E(max) = ', f10.2, /,          &
     &          5f9.4, /, 5f9.4, /, 5f9.4, /, 5f9.4, /, 5f9.4)
500     format (70('='))
501     format (2x, a45)
502     format (a70)
503     format (2x, 'R(d) = ', f8.4, 1x, 2x, ' nddke = ', i4,                &
     &          2x, ' nz_points = ', i4, 2x, ' nrhopoints = ', i4)
504     format (i5, f9.4, ' <=== ', a4)
505     format (4d20.10)
506     format (2x, ' Finished writing to file: ', a40)
601     format (' index = ', i3, 3x, '<',3i2,' |  |',3i2,'>')
 
        close (iounit)
        return
        end