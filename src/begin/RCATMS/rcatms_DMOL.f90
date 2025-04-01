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
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! rcatms_DMOL.f90
! Program Description
! ===========================================================================
!       Once the ground state wavefunctions are determined, the DMOL excited
! states are calculated (provided iexcite = 2).  The excited states are
! calculated from the pseudo-potential from the 2+ ion of the atom which
! is provided by the user.
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
        subroutine rcatms_DMOL (mesh, mesh_psirc, rc_max, rcutoff_psirc, r,  &
     &                          psi_gs)
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
        real(kind=long), intent (in), dimension (nssh, mesh) :: psi_gs

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: max_iterations = 120

! Local Variable Declaration and Description
! ===========================================================================
        integer ienergy
        integer inum
        integer ioptionpp
        integer ipoint, jpoint
        integer iremainder
        integer iteration
        integer issh
        integer jssh
        integer k
        integer l, l1, l2
        integer mesh_new

        real(kind=long) difference
        real(kind=long) dr
        real(kind=long) ebs
        real(kind=long) eee
        real(kind=long) ekin
        real(kind=long) enew
        real(kind=long) eold
        real(kind=long) eout
        real(kind=long) etot
        real(kind=long) exc
        real(kind=long) exmix
        real(kind=long) integral
        real(kind=long) total
        real(kind=long) slope
        real(kind=long) xnz
        real(kind=long) xnorm
        real(kind=long) zero

        real(kind=long), dimension (0:3,0:3,0:6) :: cgcoeff
        real(kind=long), dimension (:), allocatable :: eigenvalue
        real(kind=long), dimension (:, :), allocatable :: psi
        real(kind=long), dimension (:, :), allocatable :: psi_new
        real(kind=long), dimension (mesh) :: psi0
        real(kind=long), dimension (nssh, mesh) :: Rall
        real(kind=long), dimension (mesh) :: rho
        real(kind=long), dimension (mesh) :: rhop
        real(kind=long), dimension (mesh) :: rhopp
        real(kind=long), dimension (mesh) :: rmult
        real(kind=long), dimension (mesh) :: sigma
        real(kind=long), dimension (mesh) :: sigma_old
        real(kind=long), dimension (4) :: terms
        real(kind=long), dimension (mesh) :: uxc
        real(kind=long), dimension (:), allocatable :: uxc_local
        real(kind=long), dimension (mesh) :: v
        real(kind=long), dimension (mesh) :: vc
        real(kind=long), dimension (mesh) :: vee
        real(kind=long), dimension (nssh, mesh) :: vexx
        real(kind=long), dimension (:), allocatable :: vkin_local
        real(kind=long), dimension (:, :), allocatable :: vl
        real(kind=long), dimension (:), allocatable :: v_local
        real(kind=long), dimension (:), allocatable :: vna_local
        real(kind=long), dimension (:), allocatable :: xx
        real(kind=long), dimension (:), allocatable :: yy

        real(kind=long), external :: clebsch_gordon
! optimized w.f.
        real(kind=long), dimension (:,:), allocatable :: vXe

! Allocate Arrays
! ===========================================================================
        allocate (eigenvalue(nssh))
        allocate (psi(nssh,mesh))
        allocate (psi_new(nssh,mesh))
        allocate (uxc_local(nssh))
        allocate (vkin_local(nssh))
        allocate (vl(nssh,mesh))
        allocate (v_local(nssh))
        allocate (vna_local(nssh))
        allocate (vXe(nssh,mesh))

! Procedure
! ===========================================================================
        !write (*,*) ' Determining the wavefunctions for the 2+ ion. Once '
        !write (*,*) ' the wavefunctions are found then orthogonalize to the '
        !write (*,*) ' ground state wavefunctions to obtain a double '
        !write (*,*) ' numerical basis set. '
        !write (*,*) '  '

        !write (*,*) ' mesh = ', mesh

! Initialize the wavefunctions and the exact exchange potential.
        psi = 0.0d0
        psi_new = 0.0d0
        vexx = 0.0d0

! Get the Clebsch Gordan Coefficients.
        do l1 = 0, 3
         do l2 = 0, 3
          do k = 0, 6
           cgcoeff(l1,l2,k) = clebsch_gordon(l1,0,k,0,l2,0)**2
          end do
         end do
        end do

! Construct r at each mesh point.
        dr = rc_max/real(mesh - 1)
        !write (*,*) '  '
        !write (*,*) ' Number of points in mesh = ', mesh
        !write (*,*) ' r(mesh) = ', r(mesh), ' must = rc_max = ', rc_max
        if (abs(r(mesh) - rc_max) .gt. 1.0d-6) stop ' bad agreement between r(mesh) and rc '

        !write (*,*) '  '
        !write (*,*) ' We average the output and input charge densities at each '
        !write (*,*) ' iteration as: '
        !write (*,*) ' RHO(next iteration) = (1-BETA)*RHO(out) + BETA*RHO(in) '
        !write (*,*) ' BETA = 0.5 is almost always the best choice.'
        !write (*,*) ' You have BETA = ', beta
        !write (*,*) '  '

        !write (*,*) ' Tolerance checks closeness to self-consistentcy. It is '
        !write (*,*) ' the average difference between the old and new Etot and '
        !write (*,*) ' sum of eigenvalues, divided by the new average. '
        !write (*,*) ' The standard tolerance is 1.0d-5. '
        !write (*,*) ' You have a tolerance of ', tolerance
        !write (*,*) '  '

! Determine the non-local pseudopotential.
        xnz = real(nzval_ion, kind=long)
!        call pp(ppionfile, nssh, mesh, xnz, r, vc, vl, ioptionpp, exmix)
        call pp(1, mesh, r, vc, vl, ioptionpp, exmix)

        if (ioptionpp .ne. ioption) then
         !write (*,*) ' The exchange-correlation potential option in your '
         !write (*,*) ' pseudopotential file does not match the option chosen '
         !write (*,*) ' by running initial.x! Are you positive that you want '
         !write (*,*) ' to proceed? Enter 1-yes, 0-no. '
        ! read (*,*) ioptionpp
        ! if (ioptionpp .eq. 0) then
          !write (*,*) ' You have chosen to stop! '
          stop
        ! end if
        end if

        !write (*,*) ' This program will write out data files containing the '
        !write (*,*) ' s, p, d, or f excited wavefunctions created from the '
        !write (*,*) ' DMOL formalism. File-name convention: '
        !write (*,*) ' The name begins with the chemical symbol Z (e.g. 014 for '
        !write (*,*) ' Si or 006 for C), this is then followed by the value for '
        !write (*,*) ' Rc (for example, 500 for 5.00), and finally it is '
        !write (*,*) ' appended by .ewf#, where # is the shell number '
        !write (*,*) ' (not l quantum number). An example for Si is 014_500.ewf1'
        !write (*,*) '  '
        !write (*,*) ' You have chosen: '
        !write (*,101) filename_ewf
        !write (*,*) '  '

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
                    vXe(issh,ipoint) =                                      &
     &                 v0(jssh)*exp( -1.0d0*(rcutoff(issh) - r0(jssh))      &
     &                 / (r(ipoint)-r0(jssh)) ) /(rcutoff(issh)+0.01-r(ipoint))
                 else
                    vXe(issh,ipoint) = 0.0d0
                 endif
              end do ! enddo ipoint
           enddo ! enddo  issh
! write vX        
!           do ipoint = 1, mesh
!              write (100,*) r(ipoint),(vXe(issh,ipoint),issh=1,nssh)
!           enddo
        endif ! if(ioptim)


! Construct approximate starting potential. We do this by forming hydrogenic
! wavefunctions using an effective bohr radius a0(m).
        sigma = 0.0d0
        do issh = 1, nssh
         l = lam(issh)
         do ipoint = 1, mesh - 1
          if (r(ipoint) .lt. rcutoff_psirc(issh)) then
           psi0(ipoint) = r(ipoint)**(l+1)*exp(-r(ipoint)/a0(issh))
          else
           psi0(ipoint) = 0.0d0
          end if
         end do
         psi0(mesh) = 0.0d0

! Grossly enforce the boundary condition at rc.
! Normalize and add to total charge density.
         xnorm = sqrt(dot_product(psi0,psi0)*dr)
         psi(issh,:) = psi0/xnorm
         sigma = sigma + xocc_ion(issh)*psi(issh,:)**2
        end do

! Check that int[0,rc] sig(r) dr = z
        difference = abs(dr*sum(sigma) - real(nzval_ion, kind=long))
!dani.        if (nznuc .ne. 1 .and. difference .gt. tolerance)                    &
!dani.     &   stop ' Check: Sum of charge density - not equal to nzval_ion '

! Now find rho dependent potentials, vee and uxc.
! The subroutine get_uxc depends on the actual density rho, calculate rho
! from sigma - rho = sigma/r**2. Also, since we are doing GGA's for some
! options, calculate the first and second derivatives of rho.
        rho(2:mesh-1) = sigma(2:mesh-1)/r(2:mesh-1)**2

! endpoints
        rho(1) = 2.0d0*rho(2) - rho(3)
        rho(mesh) = 2.0d0*rho(mesh-1) - rho(mesh-2)

! derivatives
        rhop(2:mesh-1) = (rho(3:mesh) - rho(1:mesh-2))/(2.0d0*dr)
        rhopp(3:mesh-1) = (rho(4:mesh) - 2.0d0*rho(3:mesh-1)                 &
     &                                 + rho(2:mesh-2))/(dr**2)

! endpoints
        rhop(1) = 2.0d0*rhop(2) - rhop(3)
        rhopp(1) = 2.0d0*rhopp(2) - rhopp(3)
        rhop(mesh) = 2.0d0*rhop(mesh-1) - rhop(mesh-2)
        rhopp(mesh) = 2.0d0*rhopp(mesh-1) - rhopp(mesh-2)

        ienergy = 0     ! do not calculate energies
        call get_uxc (ioption, mesh, r, dr, rho, rhop, rhopp, uxc, exc,      &
     &                ienergy, exmix)
        call get_vee (mesh, dr, sigma, vee)

! We do not have wave functions at this point. So do the first iteration with
! only local exchange correlation. Initialize exchange correlation potential
! to zero. - Juergen

! *****************************************************************************
! Start iteration
! Solve Schroedinger equation for each state in turn. Also calculate core and
! valence charge density.
!
! We do each orbital issh. The wave function is PNL, which is the solution of
! the radial equation. The real wave function is PSI = RNL YLM CHI(s), where
! RNL = PNL/R, YLM is spherical harmonic, and CHI(s) is spin. The integral over
! all space of PSISTAR*PSI is 1, which translates into
! integral (0,infinity) PNL**2 dR =1.
!
! The result of PNL for orbital issh is stored in psi0, and later put into psi.
        !write (*,*) ' *--------------------* '
        !write (*,*) ' | START OF ITERATION | '
        !write (*,*) ' *--------------------* '
        !write (*,*) '  '
        !write (*,200)
        !write (*,201)

        eold = 0.0d0
        do iteration = 1, max_iterations  ! ::::::: SCF LOOP ::::::::
         sigma_old = sigma
         sigma = 0.0d0

         do issh = 1, nssh  ! -------- LOOP OVER ALL SHELLS ---------
          l = lam(issh)

! put potential in v
! Note that the factor of 2 is conversion from Hartree units to Rydberg units.
          if (ioptionpp .ne. 12) then
!            v = 2.0d0*(vc + vl(issh,:)) + vee + uxc
            v = 2.0d0*(vc + vl(issh,:)) + vee + uxc + vXe(issh,:)
          else
!           v = 2.0d0*(vc + vl(issh,:)) + vee + uxc + exmix*vexx(issh,:)
           v = 2.0d0*(vc + vl(issh,:)) + vee + uxc + exmix*vexx(issh,:)     &
      &         + vXe(issh,:)
          end if

! Solve Schroedinger equation for potential v
          psi0 = psi(issh,:)
          call psirc (mesh_psirc(issh), 0, l, rcutoff_psirc(issh), v,        &
      &               eout, psi0)
          eigenvalue(issh) = eout

! Add to charge density
          psi(issh,:) = psi0
          sigma = sigma + xocc_ion(issh)*psi(issh,:)**2

          do ipoint = 1, mesh_psirc(issh)
           Rall(issh,ipoint) = 0.0d0
           if (r(ipoint) .gt. 1.0d-4) Rall(issh,ipoint) = psi0(ipoint)/r(ipoint)
          end do
         end do  ! ---------------- END LOOP OVER SHELLS -----------------

! Mix in old and new sigma's
         sigma = (1.0d0 - beta)*sigma + beta*sigma_old

! Finished with all states of angular momentum l.
! Calculate exchange correlation and Hartree potential for new sigma
! Now find rho dependent potentials, vee and uxc.
! The subroutine get_uxc depends on the actual density rho, calculate rho
! from sigma - rho = sigma/r**2. Also, since we are doing GGA's for some
! options, calculate the first and second derivatives of rho.
         rho(2:mesh-1) = sigma(2:mesh-1)/(4.0d0*pi*r(2:mesh-1)**2)

! endpoints
         rho(1) = 2.0d0*rho(2) - rho(3)
         rho(mesh) = 2.0d0*rho(mesh-1) - rho(mesh-2)

! derivatives
         rhop(2:mesh-1) = (rho(3:mesh) - rho(1:mesh-2))/(2.0d0*dr)
         rhopp(3:mesh-1) = (rho(4:mesh) - 2.0d0*rho(3:mesh-1)                &
     &                                  + rho(2:mesh-2))/(dr**2)

! endpoints
         rhop(1) = 2.0d0*rhop(2) - rhop(3)
         rhopp(1) = 2.0d0*rhopp(2) - rhopp(3)
         rhop(mesh) = 2.0d0*rhop(mesh-1) - rhop(mesh-2)
         rhopp(mesh) = 2.0d0*rhopp(mesh-1) - rhopp(mesh-2)

         ienergy = 1     ! calculate energies
         call get_uxc (ioption, mesh, r, dr, rho, rhop, rhopp, uxc, exc,     &
     &                 ienergy, exmix)
         call get_vee (mesh, dr, sigma, vee)

! The xc energy is in units of 4*pi, convert back.
         exc = 4.0d0*pi*exc

         if (ioptionpp .eq. 12) then
          do issh = 1, nssh
           call ExxPot (mesh, nssh, mesh_psirc(issh), dr, r, Rall, lam,      &
     &                  issh, 1, nssh, xocc_ion, cgcoeff, vexx)
          end do
         end if

! Calculate the total energy:  Etot=sum(eigenval)-1/2*Vee-1/4*Vxc
! (Minus sign to convert binding energies into real energies)
! Sum over the eigenvalues:
         ebs = dot_product(xocc_ion,eigenvalue)

! Evaluate the double counted electron-electron repulsion and the exchange-
! correlation energy.
! Do the integral using a trapezoidal rule.
         eee = 0.0d0

! eee = vee/2.
! exc = integral dr*sigma*(vxc-uxc)
! Contribution from origin is zero. Find terms in the sum of eigenvalues
         v_local = 0.0d0
         uxc_local = 0.0d0
         vna_local = 0.0d0
         vkin_local = 0.0d0
         do ipoint = 2, mesh - 1
          eee = eee + 0.5d0*dr*sigma(ipoint)*vee(ipoint)

! Remember that vc and vl are in hartrees so to get rydbergs multiply by 2.0
          v_local = v_local + 2.0d0*dr*psi(:,ipoint)**2*vl(:,ipoint)
          uxc_local = uxc_local + dr*psi(:,ipoint)**2*uxc(ipoint)
          vna_local = vna_local + dr*psi(:,ipoint)**2                        &
     &                              *(2.0d0*vc(ipoint) + vee(ipoint))

! Kinetic energy
! The kinetic energy is given by integral(PSI Delta PSI d^3r) with
! PSI_nlm = R_nl*Y_lm (radial times spherical harmonic, with R = u/r := PNL/r).
! Delta is the Laplace operator. In spherical coordinates, Delta is given by
! Delta = [d^2/dr^2 + 2/r d/dr] + 1/r^2 [1/sin(theta) d/d(theta) sin(theta)
!         d/d(theta) + 1/sin^2(theta) d^2/d(phi)^2].
!
! The first [...] contains only radial derivatives and always gives
! [...]RY = 1/r d^2u/dr^2 Y =: u''/r Y.
!
! The second [...] gives different results for different Y's, depending on the
! angular momentum, in general -n u^2/r^3. As we kill the angular integral by
! normalization, we only need to know n relative to the first term for our radial
! integration.
! We get for the different kinetic energies:
! T(ss) = -integral {u(r) u''(r) dr}.
! T(pp) = -integral {u(r) [u''(r) - 2u(r)/r^2] dr}.
! T(dd) = -integral {u(r) [u''(r) - 6u(r)/r^2] dr}.
! T(ff) = -integral {u(r) [u''(r) -12u(r)/r^2] dr}.
!
! Therefore, we must first find the second derivative of u(r) with respect to r.
! This is KINTERM(l) for the different orbitals.
!
! We use a four-point interpolation formula to get the second derivative (this
! can be done by function d2u by setting norder = 3)
!
! Note: In our units (Rydbergs) hbar**2/(2*mass) = 1.d0
! It is faster to not use function d2u to get the second derivative of u(r)
! if you only want to use the 4-point formula.
          do issh = 1, nssh
           ekin = (psi(issh,ipoint-1) - 2.0d0*psi(issh,ipoint)                 &
     &                                + psi(issh,ipoint+1))/dr**2
           ekin = ekin - issh*(issh-1)*psi(issh,ipoint)/r(ipoint)**2
           vkin_local(issh) = vkin_local(issh) - dr*psi(issh,ipoint)*ekin
          end do
         end do
! end integration
         etot = ebs - eee + exc

! term 1 = sum_is {vnl(is)}
! term 2 = uxc
! term 3 = vna = vcore + vhartree
! term 4 = kinetic
         terms = 0.0d0
         do issh = 1, nssh
          terms(1) = terms(1) + xocc_ion(issh)*v_local(issh)
          terms(2) = terms(2) + xocc_ion(issh)*uxc_local(issh)
          terms(3) = terms(3) + xocc_ion(issh)*vna_local(issh)
          terms(4) = terms(4) + xocc_ion(issh)*vkin_local(issh)
         end do

         enew = abs(etot) + abs(ebs)
         difference = abs((enew - eold)/enew)

         !write (*,202) iteration, etot, ebs, eee, exc, difference
         if (difference .le. tolerance) exit
         eold = enew
        end do  ! :::::::::::::::: END OF SCF LOOP :::::::::::::::::::::
! *****************************************************************************

        do issh = 1, nssh
         !write (*,300) lam(issh), xocc_ion(issh), eigenvalue(issh)
        end do

        if (difference .gt. tolerance) then
         !write (*,*) '  '
         !write (*,*) ' No self consistency after maximum iterations = ',    &
!     &   max_iterations
         stop
        end if

        !write (*,*) '  '
        !write (*,*) ' Terms in the sum of eigenvalues: '
        !write (*,*) '  '
        !write (*,*) '  in Rydbergs:'
        total = terms(1) + terms(2) + terms(3) + terms(4)
        !write (*,*) '  '
        !write (*,*) '  vnl      = ', terms(1)
        !write (*,*) '  uxc      = ', terms(2)
        !write (*,*) '  vc + vh  = ', terms(3)
        !write (*,*) '  kinetic  = ', terms(4)
        !write (*,*) '  total    = ', total
        !write (*,*) '  '

        terms = terms*ryd
        !write (*,*) '  in eV:'
        total = terms(1) + terms(2) + terms(3) + terms(4)
        !write (*,*) '  '
        !write (*,*) '  vnl      = ', terms(1)
        !write (*,*) '  uxc      = ', terms(2)
        !write (*,*) '  vc + vh  = ', terms(3)
        !write (*,*) '  kinetic  = ', terms(4)
        !write (*,*) '  total    = ', total
        !write (*,*) '  '

! Now write out wavefunction
        !write (*,*) ' This program writes out the following data files: '
        !write (*,*) '  '
        !write (*,*) '  '
        !write (*,*) ' s.wf [p.wf...]  u_s(r) [u_p...] r=100 points '
        !write (*,*) '  '

! ***************************************************************************
!       O R T H O G O N A L I Z E
! ***************************************************************************
! First find the integral of u1(x)*psi0(x) for each shell.
        do issh = 1, nssh
         mesh_new = mesh_psirc(issh)

! Set up integration factors.
         rmult(1) = dr/3.0d0
         rmult(mesh_new) = dr/3.0d0
         do ipoint = 2, mesh_new - 1, 2
          rmult(ipoint) = 4.0d0*dr/3.0d0
         end do
         do ipoint = 3, mesh_new - 2, 2
          rmult(ipoint) = 2.0d0*dr/3.0d0
         end do

         integral = 0.0d0
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_gs(issh,ipoint)*psi(issh,ipoint)
         end do

! Now orthogonalize by taking psi1_new = psi1_old - sum*psi0
         do ipoint = 1, mesh_new
          psi_new(issh,ipoint) = psi(issh,ipoint) - integral*psi_gs(issh,ipoint)
         end do

! Check normalization
         integral = 0.0d0
         !write (*,*) '  '
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)**2
         end do
         psi_new(issh,1:mesh_new) = psi_new(issh,1:mesh_new)/sqrt(integral)

         integral = 0.0d0
         !write (*,*) '  '
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)**2
         end do
         !write (*,*) ' Checking normalization for issh = ', issh
         !write (*,*) ' integral = ', integral

! Check orthoganalization
         integral = 0.0d0
         !write (*,*) '  '
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)*psi_gs(issh,ipoint)
         end do
         !write (*,*) ' Checking orthogonalization for issh = ', issh
         !write (*,*) ' integral = ', integral
        end do

! ***************************************************************************

! ***************************************************************************
!       W R I T E    O U T    W A V E F U N C T I O N S
! ***************************************************************************
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
         zero = 0.0d0
         if(sav(issh)) write(14,104) rcutoff_psirc(issh), rc_max, zero
         if(sav(issh)) write(14,105) lam(issh)

! Find (approximately) the value at r=0
         xx(2) = r(2)*abohr
         xx(3) = r(3)*abohr

         yy(2) = psi_new(issh,2)/(abohr15*r(2))
         yy(3) = psi_new(issh,3)/(abohr15*r(3))

         slope = (yy(3) - yy(2))/(xx(3) - xx(2))
         xx(1) = 0.0d0
         yy(1) = - slope*xx(2) + yy(2)

         xx(4:mesh_new-1) = r(4:mesh_new-1)*abohr
         yy(4:mesh_new-1) = psi_new(issh,4:mesh_new-1)/(abohr15*r(4:mesh_new-1))

         xx(mesh_new) = r(mesh_new)*abohr
         yy(mesh_new) = 0.0d0

         inum = idint(real(mesh_new)/4.0d0)
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
      !dani.JOM
      !    call vnn_excite

        !write (*,*) ' Bye from RCATMS_DMOL! '
        !write (*,*) '  '


! Deallocate Arrays
! ===========================================================================
        deallocate (eigenvalue)
        deallocate (psi)
        deallocate (psi_new)
        deallocate (uxc_local)
        deallocate (vkin_local)
        deallocate (vl)
        deallocate (v_local)
        deallocate (vna_local)

! Format Statements
! ===========================================================================
100     format (2x, 4f6.2)
101     format (2x, 4(2x,a11))
102     format (2x, i5, 2x, a11)
103     format (2x, i8)
104     format (2x, 2(2x,f7.4), 2x, f7.2)
105     format (2x, i3)
200     format (2x, ' iteration ', 2x, ' e(tot) ', 2x, ' eigenvalue ', 3x, &
     &          ' u(e-e) ', 6x, ' u(xc) ', 4x, ' tolerance ')
201     format (2x, 75('='))
202     format (5x, i3, 5x, 1e11.4, 4(2x,e11.5))
300     format (2x, ' l = ', i1, ' contains ', f6.3, ' electrons --- ' ,   &
     &              ' eigenvalue = ', f12.7)
400     format (4d18.10)
401     format (2x, f6.3, 2x, f12.6)

        return
        end
