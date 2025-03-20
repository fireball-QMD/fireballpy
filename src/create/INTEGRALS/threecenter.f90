! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
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


! threecenter.f90
! Program Description
! ===========================================================================
!       Generalized three center integration routines for neutral atom 
! potential and exchange correlation.
!
!       The parameter ntheta_max=5 allows five angles in the Legendre expansion.
! Should you want to test how good that is increase ntheta_max in 
! quadrature.inc.  The zeros of the (ntheta_max+1) P_l will be used then.
!
!       The index 0..nssh(in3) numbers the shells of the neutral atom and their
! sum (index = 0). The isorp-loop covers this by looping over 
! isorp = 0, ispmax, = nssh(in3).  For the exchange correlation, we have to do 
! integrals for isorp = 1, 7.  We need to define the loop limits for 
! isorp (ispmin, ispmax).  We need to redefine the array size from 0:nsh_max 
! to 0:10 (or so).  We need to pass information about the interaction type.
! bcna: interaction = 1, xc3c: interaction = 2, denS: interaction = 3.
!
!       We now employ a different strategy from the previous version of the
! code. This is called by create.f in the section. 'Real three center 
! integrals'. The call is inside of a triple loop over in1, in2, and in3 which 
! specify which atoms for the bond charge (atoms 1 and 2) and for the neutral 
! atom (atom 3).
!
!       The job of the bcna.f is to create tables for the specified atomic 
! configuration, e.g. <Si|Al|N> for all possible geometries.  The theory behind
! three center integrals is simple. We write them as a Coulomb interaction of 
! the spherical atomic density on atom 3 with an ugly bond-charge potential 
! located between atoms 1 and 2. It is possible to show (see notes) that in the
! case of the spherical atomic density, these matrix elements can be writen in 
! one of the two forms:
!
!     type A:  V(dbc,dna,theta)=SUM(l)P_l(cos(theta))*Q_l
!     type B:  V(dbc,dna,theta)=sin(theta)*SUM(l)P_l(cos(theta))*Q_l
!
! where P_l(cos(theta) is a Legendre polynimial. Sankey's theory suggests to 
! calculate these matrix elements on a 2D grid over dbc and dna, and to fit the 
! theta dependence with the above mentioned formulas. The l sum is going up to 
! l=4. Therefore, it is possible to compute the matrix element on five angles,
! and solve a matrix equation: Q=P^(-1)V, where the matrix P is the matrix of 
! five Legendre polynomials computed on the five angles.
!
!       Here is the outline of what the code is going to do:
!
!        1. First we call Kevin Schmidt's Gauss-Legendre program to generate 
!           the angles which are the zeros of the P_(l+1).
!
!        2. Set up the dna, dbc grids with a double loop.
!
!        3. Call threecenter_integral(d1,d2,theta) and get index_max
!           matrix elements.
!
!        4. Set the loop over the five (or whatever (ntheta_max+1) is)) angles.
!           Look at nabs to decide whether V(i,j) is type A or B.
!           Now store the matrix elements in a 2-D array for this isorp.
!
!        5. Use Kevin's method to generate the Q coefficients and write them 
!           out.
!
!        6. Continue on to the next (dna,dbc) pair.
!
! ===========================================================================
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
        subroutine threecenter (in1, in2, in3, index_max, iexc,         &
     &                          interaction, nzx, nssh, n1, l1, m1, n2, &
     &                          l2, m2, rcutoff1, rcutoff2, rcutoff3,   &
     &                          atom1, atom2, atom3, what1, what2,      &
     &                          what3, dbc, dna, signature,             &
     &                          ctheta, ctheta_weights, ispmin,         &
     &                          ispmax, ispnum, iammaster, ispherical)
        use dimensions
        use quadrature
        use precision
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: in3
        integer, intent (in) :: iexc         ! type of exchange-correlation
        integer, intent (in) :: index_max    ! maximum number of matrix elements
        integer, intent (in) :: interaction  ! number for the interaction
        integer, intent (in) :: ispmax
        integer, intent (in) :: ispmin
        integer, intent (in) :: ispnum

        integer, intent (in), dimension (nspec_max) :: nssh  ! number of shells
        integer, intent (in), dimension (nspec_max) :: nzx
        integer, intent (in), dimension (inter_max) :: n1  ! left atom shell 
        integer, intent (in), dimension (inter_max) :: n2  ! right atom shell
        integer, intent (in), dimension (inter_max) :: l1  ! angular momentum 
        integer, intent (in), dimension (inter_max) :: l2 
        integer, intent (in), dimension (inter_max) :: m1  ! m-value in shell
        integer, intent (in), dimension (inter_max) :: m2

        real(kind=long), intent (in) :: dbc  ! maximum bond distance: rc1 + rc2
        real(kind=long), intent (in) :: dna  ! maximum neutral atom distance: 2*rc3
        real(kind=long), intent (in) :: rcutoff1  ! largest radius of i-th atom (in Ang.)
        real(kind=long), intent (in) :: rcutoff2  ! 1, 2, 3 = left, right, neutral atom
        real(kind=long), intent (in) :: rcutoff3

        real(kind=long), intent (in), dimension (ntheta_max) :: ctheta
        real(kind=long), intent (in), dimension (ntheta_max) :: ctheta_weights

        character (len = 2), intent (in) :: atom1
        character (len = 2), intent (in) :: atom2
        character (len = 2), intent (in) :: atom3
        character (len = 70), intent (in) :: signature
        character (len = 70), intent (in) :: what1
        character (len = 70), intent (in) :: what2
        character (len = 70), intent (in) :: what3

        logical, intent (in) :: iammaster
        logical, intent (in) :: ispherical

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ibcba
        integer iam
        integer inaba
        integer index
        integer iounit
        integer isorp
        integer itheta
        integer jtheta
        integer nphi2ba
        integer nphi
        integer nrba
        integer nr
        integer ntheta
        integer nthetaba

        integer, dimension (inter_max) :: nabs

        real(kind=long) dbcx
        real(kind=long) dnax
        real(kind=long) distance_bc
        real(kind=long) cost
        real(kind=long) pl
        real(kind=long) plm
        real(kind=long) plmm
        real(kind=long) sint
        real(kind=long) temp

        real(kind=long), dimension (ntheta_max) :: answer
! Final results after sin(theta) or cos(theta) factor.
        real(kind=long), dimension (0:10, inter_max, ntheta_max) :: ggstore
! Results from the integrator
        real(kind=long), dimension (0:10, inter_max) :: gmat
        real(kind=long), dimension (ntheta_max, inter_max, 0:10) :: qpl
        real(kind=long), dimension (3) :: rna

        character (len = 40) fname
        character (len = 12) ftype

        logical skip
 
! Procedure
! ===========================================================================
        if (interaction .eq. 1 .or. interaction .eq. 3) then
         nphi2ba = nphi2ba_na
         nrba = nrba_na
         nthetaba = nthetaba_na
        else if (interaction .eq. 2)then
         nphi2ba = nphi2ba_xc
         nrba = nrba_xc
         nthetaba = nthetaba_xc
        else
         write (*,*) ' Interaction error in threecenter'
         write (*,*) ' interaction = ', interaction
        end if
 
        iam = (in3 - 1)*nspec_max*nspec_max + (in2 - 1)*nspec_max + in1 - 1
 
        if (interaction .lt. 1 .or. interaction .gt. 3) then
         stop' Interaction error in threecenter '
        end if
 
        nr = 2*nrba + 1
        ntheta = 2*nthetaba + 1
        nphi = 2*nphi2ba + 1
        if (nr .gt. 1000) stop 'error---must re-dimension vrho'
 
! Now open the data files:
! what you get looks like this: bcna_01_01.14.06.14.dat
!                           or  xc3c_05_07.06.06.14.dat
        if (interaction .eq. 1) then
         ftype = 'coutput/bcna'
        else if (interaction .eq. 2)  then
         ftype = 'coutput/xc3c'
        else if (interaction .eq. 3) then
         if (ispherical) then
          ftype = 'coutput/deS3'
         else
          ftype = 'coutput/den3'
         end if
        end if
 
        do itheta = 1, ntheta_max    ! loop over all expansion coefficients
         do isorp = ispmin, ispmax
!         10 is a guess at the largest that ispnum could ever be
!         first part gets it out of the way of the two center integrals 
          iounit = (nspec_max*nspec_max+36) + iam*10*ntheta_max
          iounit = iounit + (itheta - 1)*ispnum + isorp
          call iofile3c (ftype, 'dat', itheta, isorp, nzx(in1),         &
     &                   nzx(in2), nzx(in3), iounit, fname, skip)
          if (skip) return 
          if (iammaster) write (*,*) ' writing to: ', fname
 
! Write out basic information to files.
! Additional writes to specify what combination of atoms for the non-elemental 
! case. 
          write (iounit, 100)
          if (interaction .eq. 1) then
           write (iounit,*) 'Matrix elements for neutral atom potential'
          else if (interaction .eq. 2) then
           write (iounit,*)                                             &
     &      ' Matrix elements for exchange-correlation (Horsfield)'
          else if (interaction .eq. 3) then
           write (iounit,*)                                             &
     &      ' Matrix elements for exchange correlation (SNXC or OLSXC)'
          end if
          write (iounit, 110) itheta, isorp
          write (iounit, 100)
          write (iounit, *) ' created by: '
          write (iounit, 120) signature
          write (iounit, 120) what1
          write (iounit, 120) what2
          write (iounit, 120) what3
          write (iounit, 100)

          write (iounit, 130) nphi, nr, ntheta
          write (iounit, 140) dbc, nbcba
          write (iounit, 140) dna, nnaba
          write (iounit, 100)

          write (iounit, 150) nzx(in1), rcutoff1, atom1
          write (iounit, 150) nzx(in2), rcutoff2, atom2
          write (iounit, 150) nzx(in3), rcutoff3, atom3
          write (iounit, 100)
         end do
        end do

! ----------------------------------------------------------------------------
! 1. Call gauss-legendre routine to set up the angles.  The Legendre 
!    polynomials go from l=0 and up, P_1 has one node, P_2 has two...
!    Sankey's theory has ntheta_max = 4, so here we have, ntheta_max + 1 = 5
! ----------------------------------------------------------------------------
! Set up nabs while we are here
        do index = 1, index_max
         nabs(index) = abs(m1(index) - m2(index))
         if (nabs(index) .gt. 2) stop 'WRONG NABS IN BCNA!!!!!'
        end do
 
! ----------------------------------------------------------------------------
! 2. Begin the big loops over dbc and dna.
! ----------------------------------------------------------------------------
! Loop over all bondcharge distances.
        do ibcba = 1, nbcba
         dbcx = real(ibcba - 1, kind=long)*dbc/real(nbcba - 1, kind=long)
 
! for all bondcharges-- we set b=dbcx/2.
         distance_bc = dbcx/2.d0

! Loop over all neutral atom distances.
! The distance is measured from the bondcharge center (b=dbcx/2)
         do inaba = 1, nnaba
          dnax = real(inaba - 1, kind=long)*dna/real(nnaba - 1, kind=long)
 
! ----------------------------------------------------------------------------
! 3. Since threecenter_integral internaly loops over ispmin to ispmax.
!    The thetas are roots of P_(ntheta_max+1)
! ----------------------------------------------------------------------------
          do itheta = 1, ntheta_max
           cost = ctheta(itheta)
           sint = sqrt(1 - cost*cost)
           rna(1) = sint*dnax
           rna(2) = 0.0d0
! -------------------------------------------------------
!           rna(3)=cost*dnax + distance_bc  ! old staff
! -------------------------------------------------------
! JOM : now the integral is done with the origin at the
!       center of the bondcharge
           rna(3)=cost*dnax               ! new staff
! ------------------------------------------------------- 
! threecenter_integral computes 3-c-integrals for ispmin...ispmax potentials 
! and for all non-zero orbital combinations 1...index_max for one fixed
! dcb, rna configuration.  (index_max = index_max3c(itype1,itype2) )
! THIS IS WHERE ALL THE TIME IS SPENT
           call threecenter_integral(dbcx, rna, rcutoff1, rcutoff2,     &
     &                               nr, ntheta, nphi, in1, in2, in3,   &
     &                               gmat, index_max, n1, l1, m1, n2,   &
     &                               l2, m2, iexc, interaction, ispmin, &
     &                               ispmax, ispherical)

! ----------------------------------------------------------------------------
! 4. Correct integrals as either type A or B
! ----------------------------------------------------------------------------
           do index = 1, index_max
            if (nabs(index) .eq. 1)then
             do isorp = ispmin, ispmax
!             type B: V=sin(theta)*Sum(l)* P*Q (nabs=1)
              if (sint .lt. 0.001d0) stop 'sin theta is zero!!!'
              ggstore(isorp,index,itheta) = gmat(isorp,index)/sint
             end do
            else
!            type A: V=Sum(l) P*Q. (nabs=0,2):do nothing
             do isorp = ispmin, ispmax
              ggstore(isorp,index,itheta) = gmat(isorp,index)
             end do
            end if
           end do
          end do

! Now all qpl coefficients are computed for a fixed dbc, dna pair, for all 
! potentials  isorp = ispmin, ispmax and for all itheta(...as). The results 
! for different matrix elements of one combination of isorp and itheta
! are written sequentially into one file.
! Begin the loop over all possible potentials on the third site (in3).
          do isorp = ispmin, ispmax

! ----------------------------------------------------------------------------
! 5. Now the time has come for the Gauss-Legendre integration.  We are still 
!    in the isorp-loop.  Since the Gauss-Legendre integration has to be done 
!    separately for each kind of potential and for each kind of orbital 
!    combination, we need to loop here once more over the number of 
!    non-vanishing three center integrals
! ----------------------------------------------------------------------------
           do index = 1, index_max
! Looping over the thetas, which are the roots of a Pl
            answer = 0.0d0
            do itheta = 1, ntheta_max 
             plmm = 1.0d0
             plm = ctheta(itheta)
             temp = ctheta_weights(itheta)*ggstore(isorp,index,itheta)
             if (abs(temp) .lt. 1.0d-10) temp = 0.0d0
             answer(1) = answer(1) + plmm*temp
             answer(2) = answer(2) + plm*temp
             do jtheta = 3, ntheta_max
              pl = (plm*ctheta(itheta)                                  &
     &                 *(2.0d0*jtheta - 3.0d0) - (jtheta - 2.0d0)*plmm) &
     &             /(jtheta - 1.0d0)
              answer(jtheta) = answer(jtheta) + pl*temp
              plmm = plm
              plm = pl
             end do 
            end do 

! Normalize the coefficient, and write them out to qpl's.
            do itheta = 1, ntheta_max
             qpl(itheta,index,isorp) = answer(itheta)*(2.d0*itheta - 1.d0)*0.5d0
            end do
           end do !  end of the GL loop
          end do !  end of the isorp loop
 
! ----------------------------------------------------------------------------
! qpl's are the answer
! ----------------------------------------------------------------------------
! Write the qpl coefficients into the data files each combination in1, in2, in3,
! itheta(=1,ntheta_max), isorp gives an individual file.  The values for the 
! different non-zero matrix elements of a given combination are written out 
! after the index loop.
! ----------------------------------------------------------------------------
          do isorp = ispmin, ispmax
           do itheta = 1, ntheta_max
            iounit = (nspec_max*nspec_max + 36) + iam*10*ntheta_max
            iounit = iounit + (itheta - 1)*ispnum + isorp
            write (iounit,200) (qpl(itheta,index,isorp), index = 1, index_max)
           end do
          end do
         end do   ! end of the dna loop
        end do! the end of the dbc loop
 
! Close files
        do itheta = 1, ntheta_max
         do isorp = ispmin, ispmax
          iounit = (nspec_max*nspec_max + 36) + iam*10*ntheta_max
          iounit = iounit + (itheta - 1)*ispnum + isorp
          close (unit = iounit)
         end do
        end do

! Format Statements
! ===========================================================================
100     format (70('='))
110     format ('   itheta:  ',i2,'  isorp:  ',i2)
120     format (a70)
130     format (3(2x, i3))
140     format (2x, f9.4, 2x, i3)
150     format (i5, f9.4, ' <=== ', a4)
200     format (4d18.8)

        return
        end