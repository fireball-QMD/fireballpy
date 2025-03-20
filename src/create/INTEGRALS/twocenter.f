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
 
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
 
! twocenter.f
! Program Description
! ===========================================================================
!      This code is a general two-center integration routine for matrix
! elements of the form <psi(1)|V(1)|psi(2)>.  Thus, V(1) is located at the
! site of one of the orbitals.  The potential V(1) is something like Vxc for
! the exchange correlation potential, Vna for the neutral atom potential, or
! 1 for the overlap term.
!
! ===========================================================================
! Code written by:
! John Tomfohr
! Campus Box 871504
! Department of Physics
! Arizona State Univerity
! Tempe, AZ 85287-1504
! (602) 965-0667 (office)      email: tomfohr@asu.edu
!
! with modifications by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine twocenter (interaction, isorp, ideriv, iexc,
     1                        fraction, itype1, itype2, atom1, atom2,
     2                        what1, what2, nzx1, nzx2, rcutoff1,
     3                        rcutoff2, nssh2, nz, nrho, ndd,
     4                        index_max, nleft, lleft, mleft, nright, 
     5                        lright, mright, signature, iammaster, 
     6                        ispher)
        implicit none
 
        include '../parameters.inc'
        include '../pseudopotentials.inc'
 
! Argument Declaration and Description
! ===========================================================================
        integer ideriv          ! for xc interactions need derivative type
        integer iexc            ! type of exchange-correlation approximation
        integer index_max       ! maximum number of non-zero elements
        integer isorp
        integer itype1
        integer iemergency
        integer itype2
 
        integer interaction     ! interaction type
                                !   0 = density_OSLXC
                                !   1 = overlap
                                !   2 = vna ontop
                                !   3 = vna atom
                                !   4 = non-local
                                !   5 = xc ontop
                                !   6 = xc atom-atom
                                !   7 = xc correction
                                !   8 = z-dipole
                                !   9 = y-dipole
                                !  10 = x-dipole
                                !  11 = coulomb
                                !  12 = extended hubbard (n1*n2 dmuxc(n1+n2)/dn)
                                !  13 = extended hubbard spin dependent term
                                !  14 = olscx (charge transfer correction)
 
        integer ndd             ! number of grid points along dna
        integer nssh2
        integer nrho            ! number of rho-points on grid
        integer nz              ! number of z-points on grid
 
        integer nzx1            ! nuclear charges
        integer nzx2
 
! arrays which determine non-zero matrix elements
        integer lleft  (inter_max)
        integer lright (inter_max)
        integer mleft  (inter_max)
        integer mright (inter_max)
        integer nleft  (inter_max)
        integer nright (inter_max)
 
        real*8 fraction
        real*8 rcutoff1
        real*8 rcutoff2
 
        character*20 root
        character*70 signature
 
        character*2 atom1
        character*2 atom2
        character*70 what1
        character*70 what2
 
        logical ispher
        logical iammaster
        logical skip
 
! Local Variable Declaration and Description
! ===========================================================================
        integer igrid
        integer index
        integer index_coulomb
        integer issh
        integer l1, l2
        integer m1,m2
        integer n1, n2
        integer iounit
 
        real*8 d
        real*8 dmax
        real*8 dr
        real*8 sum
 
        real*8 coulomb_hold (inter_max)
        real*8 hold (inter_max)
        real*8 ahold (nssh2)
        integer lsh (nssh2)
 
        character*40 fname
 
! Procedure
! ===========================================================================
        iounit = ((itype2 - 1)*nspec_max + itype1 - 1) + 36
        if (iammaster) then
         write (*,*) '  '
         write (*,100)
         write (*,*) ' Welcome to the general two-center matrix element'
         write (*,*) ' calculator.  This routine can calculate matrix '
         write (*,*) ' elements of s, p, d, and f interactions. '
         write (*,*) ' Currently the flag interaction = ', interaction
         write (*,*) '  '
         write (*,*) ' This value of the interaction flag implies '
        end if ! end master
 
! *************************************************************************
! ontop => orbitals at two different sites
! atom-atom => orbitals at same site
! interaction = 0 => do density_oslxc
! interaction = 1 => do overlap
! interaction = 2 => do neutral atom/ontop function => vneutral
! interaction = 3 => do neutral atom/atom
! interaction = 4 => do non-local
! interaction = 5 => do xc ontop function => vxc
! interaction = 6 => do xc atom-atom function => dvxc
! interaction = 7 => do xc double counting correction function => dexc
! interaction = 8 => do z-dipole
! interaction = 9 => do y-dipole
! interaction = 10 => do x-dipole
! interaction = 11 => do coulomb
! interaction = 12 => do extended hubbard (n1*n2 dmuxc(n1+n2)/dn)
! interaction = 13 => do extended hubbard spin dependent term
! interaction = 14 => do McWeda charge transfer correction terms
! *************************************************************************
        if (interaction .eq. 0) then
         if (iammaster) then
          write (*,*) ' that the average density matrix elements are '
          write (*,*) ' now being calculated. '
          write (*,100)
         end if ! end master
        if (ideriv .eq. 0) then
          ! spherical approx. 
          if(ispher) then 
             root = 'coutput/denS_atom'
          else
          ! NOT spherical approx. 
             root = 'coutput/den_atom'
          endif
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit, 
     1                       fname, skip)
         endif
         if (ideriv .eq. 1) then
          ! spherical approx. 
          if(ispher) then 
             root = 'coutput/denS_ontopl'
          else
          ! NOT spherical approx. 
             root = 'coutput/den_ontopl'          
          endif
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit, 
     1                        fname, skip)
         end if
         if (ideriv .eq. 2) then
          ! spherical approx. 
          if(ispher) then 
             root = 'coutput/denS_ontopr'
          else
          ! NOT spherical approx. 
             root = 'coutput/den_ontopr'
          endif
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit,
     1                        fname, skip)
         endif
        else if (interaction .eq. 1) then
         if (iammaster) then
          write (*,*) ' that the overlap matrix elements are now '
          write (*,*) ' being calculated. '
          write (*,100)
         end if ! end master
          ! spherical approx. 
          if(ispher) then 
             root = 'coutput/overlapS'
          else
          ! NOT spherical approx. 
             root = 'coutput/overlap'
          endif
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 2) then
         if (iammaster) then
          write (*,*) ' that the (non)-neutral atom potential matrix '
          write (*,*) ' elements (ontop) are now being calculated. '
          if (isorp .eq. 0) then
           write (*,*) ' Doing neutral potential, isorp = ', isorp
          else
           write (*,*) ' Doing non-neutral potential, shell = ', isorp
          end if
          write (*,100)
         end if ! end master
         if (ideriv .eq. 1) then
          root = 'coutput/vna_ontopl'
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit, 
     1                        fname, skip)
         end if
         if (ideriv .eq. 2) then
          root = 'coutput/vna_ontopr'
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit,
     1                        fname, skip)
         end if
        else if (interaction .eq. 3) then
         if (iammaster) then
          write (*,*) ' that the (non)-neutral atom potential matrix '
          write (*,*) ' elements (atom-atom) are now being calculated. '
          write (*,100)
          if (isorp .eq. 0) then
           write (*,*) ' Doing neutral potential, isorp = ', isorp
          else
           write (*,*) ' Doing non-neutral potential, shell = ', isorp
          end if
         end if ! end master
         root = 'coutput/vna_atom' 
         call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit,
     1                       fname, skip)
        else if (interaction .eq. 4) then
         if (iammaster) then
          write (*,*) ' that the non-local potential matrix elements '
          write (*,*) ' are now being calculated. '
          write (*,100)
         end if ! end master
         root = 'coutput/vnl' 
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 5) then
         if (iammaster) then
          write (*,*) ' that the exchange-correlation matrix '
          write (*,*) ' elements (ontop) are now being calculated. '
          write (*,100)
         end if ! end master
         root = 'coutput/xc_ontop'
         call iofile2c_xcna (root, 'dat', ideriv, nzx1, nzx2, iounit, 
     1                       fname, skip)
        else if (interaction .eq. 6) then
         if (iammaster) then
          write (*,*) ' that the exchange-correlation matrix '
          write (*,*) ' elements (atom-atom) are now being calculated. '
          write (*,100)
         end if ! end master
         root = 'coutput/xc_atom'
         call iofile2c_xcna (root, 'dat', ideriv, nzx1, nzx2, iounit,
     1                       fname, skip)
        else if (interaction .eq. 7) then
         if (iammaster) then
          write (*,*) ' that the exchange-correlation correction from '
          write (*,*) ' double counting is now being calculated. '
          write (*,100)
         end if ! end master
         root = 'coutput/xc_corr'
         call iofile2c_xcna (root, 'dat', ideriv, nzx1, nzx2, iounit,
     1                       fname, skip)
        else if (interaction .eq. 8) then
         if (iammaster) then
          write (*,*) ' that the z-dipole correction is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/dipole_z'
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 9) then
         if (iammaster) then
          write (*,*) ' that the y-dipole correction is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/dipole_y'
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 10) then
         if (iammaster) then
          write (*,*) ' that the x-dipole correction is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/dipole_x'
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 11) then
         if (iammaster) then
          write (*,*) ' that the coulomb correction is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/coulomb'
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 12) then
         if (iammaster) then
          write (*,*) ' that the nu extended hubbard xc is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/nuxc' 
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
        else if (interaction .eq. 13) then
         if (iammaster) then 
          write (*,*) ' that the hubbard xc spin dependent term '     
          write (*,*) ' now being calculated.'
          write (*,100)
         end if ! end master
         root = 'coutput/nuxcs'
         call iofile2c (root, 'dat', nzx1, nzx2, iounit, fname, skip)
! jel-X
        else if (interaction .eq. 14) then
         if (iammaster) then
          write (*,*) ' that the mc-weda nuxc2crho is now being '
          write (*,*) ' calculated.'
          write (*,100)
         end if ! end master
         if (ideriv .eq. 1) then 
          root = 'coutput/dnuxc_ol' 
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit,
     1                        fname, skip)
         endif
         if (ideriv .eq. 2) then 
          root = 'coutput/dnuxc_or' 
          call iofile2c_xcna (root, 'dat', isorp, nzx1, nzx2, iounit, 
     1                        fname, skip)
         endif
! end jel-X
        end if
        if (skip) return
 
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' Integration is done in cylinder coordinates. '
         write (*,300) nz, nrho
        end if ! end master
 
! *************************************************************************
! Check the dimensions.
        dmax = rcutoff1 + rcutoff2
        if (iammaster) then
         write (*,*) '  '
         write (*,*) ' The maximum distance d (sum of Fireball radii) '
         write (*,*) ' between the atom and the neutral atom potential '
         write (*,*) ' (in angstroms) is ', dmax
         write (*,*) '  '
         write (*,*) ' ndd: '
         write (*,*) ' The number of points along dna should be >= 107.'
         write (*,*) ' ndd = ', ndd
         if (ndd .lt. 107) then
          write (*,*) ' You must have at LEAST 107 points. Sorry. '
          write (*,*) ' Fixup quadrature.inc: Make ndd >= 107. '
          stop ' Fixup quadrature.inc: Make ndd >= 107. '
         end if
        end if ! end master
        dr = dmax/dfloat(ndd - 1)
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
        end if ! end master
 
! *************************************************************************
! Depending on the interaction type, then open the appropriate output files.
! Write out initial information.
! *************************************************************************
! Set up the header strings:
        write (iounit,100)
        if (interaction .eq. 0) then
         write (iounit,*) ' Oslxc density matrix elements '
        else if (interaction .eq. 1) then
         write (iounit,*) ' Overlap matrix elements '
        else if (interaction .eq. 2) then
         write (iounit,*)
     1             ' (Non)-neutral atom potential (ontop) two-center'
        else if (interaction .eq. 3) then
         write (iounit,*)
     1             ' (Non)-neutral atom potential (atom) two-center'
        else if (interaction .eq. 4) then
         write (iounit,*) ' Non-local pseudopotential matrix elements '
        else if (interaction .eq. 5) then
         write (iounit,*)' Matrix elements for the xc potential (ontop)'
        else if (interaction .eq. 6) then
         write (iounit,*) ' Matrix elements for the xc potential (atom)'
        else if (interaction .eq. 7) then
         write (iounit,*) ' The xc correction for over counting '
        else if (interaction .eq. 8) then
         write (iounit,*) ' Matrix elements for the z-dipole '
        else if (interaction .eq. 9) then
         write (iounit,*) ' Matrix elements for the y-dipole '
        else if (interaction .eq. 10) then
         write (iounit,*) ' Matrix elements for the x-dipole '
        else if (interaction .eq. 11) then
         write (iounit,*) ' Matrix elements for the short-range coulomb'
        else if (interaction .eq. 12) then
         write (iounit,*)
     1          ' Extended Hubbard: nu12 = int[n1*n2*nu(n1+n2)]'
        else if (interaction .eq. 13) then
         write (iounit,*) 
     1          ' Spin extended Hubbard: nu12s=int[n1*n2*nus(n1+n2)]'
        else if (interaction .eq. 14) then
         write (iounit,*) 
     1          ' McWeda charge transfer: nu12 = int[n1*n2*nu(n1+n2)] ' 
        end if
 
        write (iounit,*) ' created by: '
        write (iounit,420) signature
        write (iounit,110) what1
        write (iounit,110) what2
        write (iounit,100)
 
        write (iounit,500) dmax, ndd, nz, nrho
        write (iounit,100)
 
        write (iounit,600) nzx1, rcutoff1, atom1
        write (iounit,600) nzx2, rcutoff2, atom2
 
        if (interaction .eq. 4) then
         write (iounit,*) nssh2
         write (iounit,*) (cl_pp(issh,itype2), issh = 1, nssh2)
        end if
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
!          ----------------------------------------------------------
!          |             Magic formula table                        |
!          ----------------------------------------------------------
!          |             s(sigma)  = 1                              |
!          |             p_sigma   = z/r                            |
!          |             p_pi      = rho/r                          |
!          |             d_sigma   = (2*z**2-rho**2)/r**2           |
!          |             d_pi      = rho*z/r**2                     |
!          |             d_delta   = rho**2/r**2                    |
!          |             f_sigma   = z*(2*z**2-3*rho**2)/r**3       |
!          |             f_pi      = rho*(4*z**2-rho**2)/r**3       |
!          |             f_delta   = z*rho**2/r**3                  |
!          |             f_phi     = rho**3/r**3                    |
!          ----------------------------------------------------------
! *************************************************************************
! Initialize d (which is along z-axis)
        d = - dr
 
! Calculate only the non-zero matrix elements. This occurs when m = m'.
! The variables lleft, mleft contains the angular momentum of the left
! wavefunction ("bra") and lright, lright the angular momentum of the right
! wavefunction ("ket").  If and only if mleft = mright, then the matrix
! element is non-zero - call the integration routine.
 
! Loop over the grid points which define the distances between the two
! orbitals, i.e. this is dna.
        do igrid = 1, ndd
         d = d + dr
         if (interaction .eq. 5 .or. interaction .eq. 6 .or. 
     1       interaction .eq. 7 .or. interaction .eq. 12
     2       .or. interaction .eq. 13 ) then
          call rho2c_store (iexc, itype1, itype2, rcutoff1, rcutoff2, d,
     1                      ideriv + 1)
         end if
         if (interaction .eq. 14 ) then
          call rho2c_store (iexc, itype1, itype2, rcutoff1, rcutoff2, d, 
     1                              1)
         end if
         do index = 1, index_max
          n1 = nleft(index)
          l1 = lleft(index)
          m1 = mleft(index)
          n2 = nright(index)
          l2 = lright(index)
          m2 = mright(index)
 
! For coulomb and extended hubbard, we only do those where mleft=mright=0
! (See below).
          if (interaction .eq. 11 .or. interaction .eq. 12 .or.
     1        interaction .eq. 13 ) then
           if (mleft(index) .eq. 0 .and. mright(index) .eq. 0) then
            call twocenter_integral (interaction, isorp, ideriv, iexc,
     1                               fraction, n1, l1, m1, n2, l2, m2,
     2                               nz, nrho, d, itype1, itype2,
     3                               rcutoff1, rcutoff2, sum, ispher)
           end if
          else if (interaction .ne. 11 .and. interaction .ne. 12 .and.
     1             interaction .ne. 13) then
            call twocenter_integral (interaction, isorp, ideriv, iexc,
     1                              fraction, n1, l1, m1, n2, l2, m2,
     2                              nz, nrho, d, itype1, itype2,
     3                              rcutoff1, rcutoff2, sum, ispher)
          end if
          if (abs(sum) .lt. 1.0d-10) sum = 0.0d0
          hold(index) = sum
         end do

! symmetrize diagonal terms at distance = 0

         if ((nzx1 .eq. nzx2) .and. (igrid .eq. 1)) then
!JOM
         if (interaction .ne. 11 .and. interaction .ne. 12 .and.
     1             interaction .ne. 13) then
          lsh = 0
          ahold = 0.0d0
          do index = 1, index_max

           n1 = nleft(index)
           l1 = lleft(index)
           m1 = mleft(index)
           n2 = nright(index)
           l2 = lright(index)
           m2 = mright(index)
           if (n1 .eq. n2 .and. l1 .eq. l2 .and. m1 .eq. m2) then
            ahold(n1) = ahold(n1) + hold(index)
            lsh(n1) = lsh(n1) + 1
           endif

          end do ! index
          do issh = 1,nssh2
           ahold(issh) = ahold(issh)/lsh(issh)
          enddo

          do index = 1, index_max

           n1 = nleft(index)
           l1 = lleft(index)
           m1 = mleft(index)
           n2 = nright(index)
           l2 = lright(index)
           m2 = mright(index)
           if (n1 .eq. n2 .and. l1 .eq. l2 .and. m1 .eq. m2) then
            hold(index) = ahold(n1)
           endif
          end do ! index

         endif !JOM
         endif !if (nzx1 .and. igrid)
 
         if (interaction .eq. 11 .or. interaction .eq. 12 . or.
     1       interaction .eq. 13 ) then
          index_coulomb = 0
          do index = 1, index_max
           if (mleft(index) .eq. 0 .and. mright(index) .eq. 0) then
            index_coulomb = index_coulomb + 1
            coulomb_hold(index_coulomb) = hold(index)
           end if
          end do
          write(iounit,800)(coulomb_hold(index),index = 1,index_coulomb)
         else
          write (iounit,800) (hold(index), index = 1, index_max)
         end if
        end do
 
        if (iammaster) then
         write (*,*) '  '
         write (*,700) fname
         write (*,*) ' Exchange correlation model = ', iexc
         write (*,*) '  '
         write (*,100)
         write (*,*) '  '
        end if ! end master
 
! Format Statements
! ===========================================================================
100     format (70('='))
110     format (a70)
300     format (2x, 'Number of points for nz, nrho grids =', 2i5)
301     format (2x, 'point separation is ', f9.6, 1x, '(A)')
420     format (2x, a30)
500     format (2x, 'R(d) = ', f8.4, '(A)', 2x, ' ndd = ', i4,
     1          2x, ' nz_points =', i4, 2x, ' nrhopoints = ', i4)
600     format (i5, f9.4, ' <=== ', a4)
700     format (2x, 'Finished. Wrote output to ', a40)
800     format (4d18.8)
        close (unit = iounit)
        return
        end
