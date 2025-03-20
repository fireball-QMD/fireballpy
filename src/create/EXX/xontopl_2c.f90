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

 
! xontopl_2c.f90
! Program Description
! ===========================================================================
!
!       This routine calculates the two-center integrals for the exact
! exchange ontop(left) interactions.
!
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
        subroutine xontopl_2c (nspec_max, isorp, fraction, nsshxc, &
     &                         itype1, itype2, atom1, atom2, what1, &
     &                         what2, nzx1, nzx2, rcutoff1, rcutoff2, nz,   &
     &                         nrho, ndd, index_max, inter_max, nalpha,     &
     &                         lalpha, malpha, nleft, lleft, mleft, nright, &
     &                         lright, mright, signature, iammaster)
        use precision
        use x_exact
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: index_max
        integer, intent (in) :: inter_max
        integer, intent (in) :: isorp
        integer, intent (in) :: itype1
        integer, intent (in) :: itype2
        integer, intent (in) :: lalpha
        integer, intent (in) :: malpha
        integer, intent (in) :: nalpha
 
! arrays which determine non-zero matrix elements
        integer lleft  (inter_max)
        integer lright (inter_max)
        integer mleft  (inter_max)
        integer mright (inter_max)
        integer nleft  (inter_max)
        integer nright (inter_max)
 
        integer, intent (in) :: ndd
        integer, intent (in) :: nrho
        integer, intent (in) :: nspec_max
        integer, intent (in) :: nz
        integer, intent (in) :: nzx1
        integer, intent (in) :: nzx2
 
        integer, intent (in), dimension (nspec_max) :: nsshxc
 
        real(kind=long), intent (in) :: fraction
        real(kind=long), intent (in) :: rcutoff1
        real(kind=long), intent (in) :: rcutoff2
 
        character (len=2), intent (in) :: atom1
        character (len=2), intent (in) :: atom2
        character (len=70), intent (in) :: signature
        character (len=70), intent (in) :: what1
        character (len=70), intent (in) :: what2
 
        logical iammaster
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: lmax = 3
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iemergency
        integer igrid
        integer index
        integer iounit
        integer l1, l2
        integer m1, m2
        integer n1, n2
 
        real(kind=long) d
        real(kind=long) dmax
        real(kind=long) dr
        real(kind=long) sum
 
        real(kind=long), dimension (inter_max) :: hold
 
        character (len=40) fname
 
! Procedure
! ===========================================================================
! Write out the initial message.
        if (iammaster) then
         write (*,*) '  '
         write (*,100)
         write (*,*) ' Welcome to the two-center matrix element routine which '
         write (*,*) ' calculates the exact exchange interactions for the '
         write (*,*) ' ontop (left) case.  This routine can calculate matrix '
         write (*,*) ' elements of s, p, d, and f interactions. '
         write (*,*) '  '
        end if
 
! Initialize
        iounit = ((itype2 - 1)*nspec_max + itype1 - 1) + 36
 
! Open the file to store the onecenter data.
        call iofile2c_x ('coutput/xontopl_2c', 'dat', isorp, nzx1, nzx2, fname)
        open (unit = iounit, file = fname, status = 'unknown')
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Integration is done in cylinder coordinates. '
         write (*,300) nz, nrho
        end if
 
! *************************************************************************
! Check the dimensions.
        dmax = rcutoff1 + rcutoff2
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Maximum distance d (sum of Fireball radii) between the '
         write (*,*) ' atom and the exchange potential (angstroms) is ', dmax
         write (*,*) '  '
         write (*,*) ' ndd: '
         write (*,*) ' The number of points along dna should be >= 107. '
         write (*,*) ' ndd = ', ndd
         if (ndd .lt. 107) then
          write (*,*) ' You must have at LEAST 107 points. Sorry. '
          write (*,*) ' Fixup quadrature.inc: Make ndd >= 107. '
          write (*,*) ' I will give you an emergency override. Insert -99 now! '
! No override for parallel compiled code.
!$        stop ' Fixup quadrature.inc: Make ndd >= 107. '
          read (*,*) iemergency
          if (iemergency .eq. -99) then
           write (*,*) ' Emergency (emergency) override given. '
          else
           stop ' Fixup quadrature.inc: Make ndd >= 107. '
          end if
         end if
        end if
        dr = dmax/real(ndd - 1)
        if (iammaster) then
         write (*,301) dr
         write (*,*) '  '
         write (*,*) ' nz: '
         write (*,*) ' The number of points is 2*nz + 1 (>= 96) '
         write (*,*) ' nz = ', nz
         write (*,*) '  '
         write (*,*) ' nrho: '
         write (*,*) ' The number of points is 2*nrho + 1 (>= 96) '
         write (*,*) ' nrho = ', nrho
        end if
 
! Set up the header for the output file.
        write (iounit,100)
        write (iounit,*) ' Two-center exchange ontop (left) interactions: '
        write (iounit,*) ' created by: '
        write (iounit,200) signature
        write (iounit,110) what1
        write (iounit,110) what2
        write (iounit,100)
 
        write (iounit,500) dmax, ndd, nz, nrho
        write (iounit,100)
 
        write (iounit,600) nzx1, rcutoff1, atom1
        write (iounit,600) nzx2, rcutoff2, atom2
        write (iounit,*) dmax, ndd
 
! *************************************************************************
! Do the interactions: ss, sp_sig, ps_sig, pp_pi, pp_sig, sd_sig, ...
! *************************************************************************
! Calculate the matrix elements for each point (distance between atoms)
! between 0 and dmax angstroms.
!
! *************************************************************************
! I call the angular part of an orbital sigma, pi, delta, or phi if it was
! produced by taking linear combinations of spherical harmonics with
! |m|=0, 1, 2, or 3 (m is the z angular momentum quantum number and the
! z-axis points from atom 1 to atom 2):
! *************************************************************************
!
!  orbital name and f(x,y,z)           | normalization factor
!--------------------------------------|----------------------------------*
! s (sigma)    1                       |   sqrt(1/4pi)
!--------------------------------------|----------------------------------*
! p_sigma -->  z                       |
!                                      |
!            [ x                       |   sqrt(3/4pi)/r
!    p_pi -->[                         |
!            [ y                       |     (same for all)
!--------------------------------------|----------------------------------*
! d_sigma -->  (3*z**2-r**2)/sqrt(12)  |
!                                      |
!            [ z*x                     |
!    d_pi -->[                         |
!            [ z*y                     |   sqrt(15/4pi)/r**2
!                                      |
!            [ x*y                     |     (same for all)
! d_delta -->[                         |
!            [ (x**2-y**2)/2           |
!--------------------------------------|----------------------------------*
! f_sigma -->  z*(5*z**2-3*r**2)       |   sqrt(7/16/pi)/r**3
!                                      |
!            [ x*(5*z**2-r**2)         |
!    f_pi -->[                         |   sqrt(21/32/pi)/r**3
!            [ y*(5*z**2-r**2)         |     (for both f_pi)
!                                      |
!            [ x*y*z                   |
! f_delta -->[                         |   sqrt(105/4pi)/r**3
!            [ z*(x**2-y**2)/2         |     (for both f_delta)
!                                      |
!            [ x**3-3*x*y**2           |
!   f_phi -->[                         |   sqrt(35/32/pi)/r**3
!            [ y**3-3*y*x**2           |     (for both f_phi)
!--------------------------------------|----------------------------------*
!
! These matrix elements are of the form
!
!       <psi1|v(1)|psi2> = <orbital1*R1|v(1)|orbital2*R2>,
!
! where orbital1 and orbital2 are the functions related to the angular
! components given in the table above. The functions R1 and R2 are the radial
! parts (stored elsewhere). The integration is done in cylindrical coordinates:
! we take the origin to be at the first atom with the z-axis pointing to the
! second atom and
 
!       x ---> rho*cos(phi)
!       y ---> rho*sin(phi)
!       r ---> sqrt(rho**2+z**2)
 
! Also, note that for a point with coordinates (rho,z,phi) (measured relative
! to atom 1), the value of psi1 is psi(rho,z,phi) but the value of psi2 is
! psi2(rho,z-d,phi), where d is the distance between the atoms.
 
! The integral over phi has been done by hand for each matrix element - e.g.
! for the pp_pi matrix element:
 
!       pp_pi = integral[(3/4pi)*x**2*R1*R2*rho*drho*dphi*dz]
!             = integral[(3/4pi)*(rho*cos(phi))**2*R1*R2*rho*drho*dphi*dz]
!             = integral[(3/4pi)*pi*rho**2*R1*R2*rho*drho*dz]
 
! So we're left with 3/4 times a 2-d integral (over rho and z).
! In this way we get:
!
!                >>> THE MAGIC FORMULA <<<
 
! The i'th matrix element XY_Z - where X,Y = s, p, d, or f and Z = sigma, pi,
! delta, or phi - is:
!
!      faktor(i)*integral[X_Z(1)*Y_Z(2)*R1*R2*rho*drho*dz],
!
! where the 1, 2 mean the orbital is on atom 1,2. The faktor(i) are given
! below (e.g., faktor(5) is for the matrix element pp_pi and is 3/4)
! and the X_Z and Y_Z are:
!
!          -------------------------------------------------------
!          |                Magic formula table                  |
!          -------------------------------------------------------
!          |             s(sigma)  = 1                           |
!          |             p_sigma   = z/r                         |
!          |             p_pi    = rho/r                         |
!          |             d_sigma   = (2*z**2-rho**2)/r**2        |
!          |             d_pi    = rho*z/r**2                    |
!          |             d_delta   = rho**2/r**2                 |
!          |             f_sigma   = z*(2*z**2-3*rho**2)/r**3    |
!          |             f_pi      = rho*(4*z**2-rho**2)/r**3    |
!          |             f_delta   = z*rho**2/r**3               |
!          |             f_phi     = rho**3/r**3                 |
!          -------------------------------------------------------
! *************************************************************************
! Initialize d (which is along z-axis)
        d = - dr
 
! Calculate only the non-zero matrix elements. This occurs when m = m'.
! The variables lleft, mleft contains the angular momentum of the left
! wavefunction ("bra") and lright, lright the angular momentum of the right
! wavefunction ("ket").  If and only if mleft = mright, then the matrix
! element is non-zero - call the integration routine.
 
! First calculate the r' integral. This integral can be evaluates only once
! for each distance d, thus saving a lot of computational time. This is
! possible because this integral is independent of the distance d between
! the centers.
        call xontopl_2c_rprime (nspec_max, nsshxc, nalpha, itype1, rcutoff1, &
    &                           nrho, lmax)
 
! Loop over the grid points which define the distances between the two
! orbitals, i.e. this is dna.
        do igrid = 1, ndd
         d = d + dr
         do index = 1, index_max
          n1 = nleft(index)
          l1 = lleft(index)
          m1 = mleft(index)
          n2 = nright(index)
          l2 = lright(index)
          m2 = mright(index)
          call xontopl_2c_integral (fraction, nalpha, lalpha, malpha, n1,    &
     &                              l1, m1, n2, l2, m2, nz, nrho, d, itype1, &
     &                              itype2, rcutoff1, rcutoff2, lmax, sum)
          hold(index) = sum
         end do
!        write (iounit,800) (hold(index), index = 1, index_max)
         write (iounit,800) d, (hold(index), index = 1, index_max)
        end do
 
        if (iammaster) then
         write (*,*) '  '
         write (*,700) fname
         write (*,*) '  '
         write (*,100)
         write (*,*) '  '
        end if
 
! Deallocate Arrays
! ===========================================================================
        deallocate (rpoint)
        deallocate (rprime)
 
! Format Statements
! ===========================================================================
100     format (70('='))
110     format (a70)
200     format (2x, a30)
300     format (2x, 'Number of points for nz, nrho grids =', 2i5)
301     format (2x, 'point separation is ', f9.6, 1x, '(A)')
500     format (2x, 'R(d) = ', f8.4, '(A)', 2x, ' ndd = ', i4,             &
     &          2x, ' nz_points =', i4, 2x, ' nrhopoints = ', i4)
600     format (i5, f9.4, ' <=== ', a4)
700     format (2x, 'Finished. Wrote output to ', a40)
800     format (10d20.10)
 
        close (unit = iounit)
        return
        end