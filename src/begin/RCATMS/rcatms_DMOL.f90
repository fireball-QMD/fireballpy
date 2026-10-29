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
        use, intrinsic :: iso_fortran_env, only: dp => real64
	use :: utils, only: utils_itos, utils_open, utils_log, UTILS_LOG_ERROR
        use :: xc, only: xc_calc, xc_isgga
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: mesh

        integer, intent (in), dimension (nssh) :: mesh_psirc

        real(kind=dp), intent (in) :: rc_max

        real(kind=dp), intent (in), dimension (mesh) :: r
        real(kind=dp), intent (in), dimension (nssh) :: rcutoff_psirc
        real(kind=dp), intent (in), dimension (nssh, mesh) :: psi_gs

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
        integer k, io
        integer l, l1, l2
        integer mesh_new

        real(kind=dp) difference
        real(kind=dp) dr
        real(kind=dp) ebs
        real(kind=dp) eee
        real(kind=dp) ekin
        real(kind=dp) enew
        real(kind=dp) eold
        real(kind=dp) eout
        real(kind=dp) etot
        real(kind=dp) exc
        real(kind=dp) exmix
        real(kind=dp) integral
        real(kind=dp) total
        real(kind=dp) slope
        real(kind=dp) xnz
        real(kind=dp) xnorm
        real(kind=dp) zero
        real(kind=dp) dexcrho, dexcsigma, duxcrho, duxcsigma, duxclapl, duxccross

        real(kind=dp), dimension (0:3,0:3,0:6) :: cgcoeff
        real(kind=dp), dimension (:), allocatable :: eigenvalue
        real(kind=dp), dimension (:, :), allocatable :: psi
        real(kind=dp), dimension (:, :), allocatable :: psi_new
        real(kind=dp), dimension (mesh) :: psi0
        real(kind=dp), dimension (nssh, mesh) :: Rall
        real(kind=dp), dimension (mesh) :: rho
        real(kind=dp), dimension (mesh,1) :: rhop
        real(kind=dp), dimension (mesh,1,1) :: rhopp
        real(kind=dp), dimension (mesh) :: rmult
        real(kind=dp), dimension (mesh) :: sigma
        real(kind=dp), dimension (mesh) :: sigma_old
        real(kind=dp), dimension (4) :: terms
        real(kind=dp), dimension (mesh) :: uxc, eexc
        real(kind=dp), dimension (:), allocatable :: uxc_local
        real(kind=dp), dimension (mesh) :: v
        real(kind=dp), dimension (mesh) :: vc
        real(kind=dp), dimension (mesh) :: vee
        real(kind=dp), dimension (nssh, mesh) :: vexx
        real(kind=dp), dimension (:), allocatable :: vkin_local
        real(kind=dp), dimension (:, :), allocatable :: vl
        real(kind=dp), dimension (:), allocatable :: v_local
        real(kind=dp), dimension (:), allocatable :: vna_local
        real(kind=dp), dimension (:), allocatable :: xx
        real(kind=dp), dimension (:), allocatable :: yy

        real(kind=dp), external :: clebsch_gordon
! optimized w.f.
        real(kind=dp), dimension (:,:), allocatable :: vXe

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
! Initialize the wavefunctions and the exact exchange potential.
        psi = 0.0_dp
        psi_new = 0.0_dp
        vexx = 0.0_dp

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
	  
! Determine the non-local pseudopotential.
        xnz = real(nzval_ion, kind=dp)
        call pp(1, mesh, r, vc, vl, ioptionpp, exmix)

! create additional X-potential for excited orbitals to optimize 
! pseudoatomic w.f.
        vXe = 0.0_dp
        if (ioptim .eq. 1) then 
           do issh = 1,nssh
              jssh = issh + nssh
              do ipoint = 1, mesh
! (r/r0)**6.0
!                  vXe(issh,ipoint) = (r(ipoint)/r0(issh))**2.0_dp
                 if (r(ipoint) .gt. r0(jssh)) then 
! siesta formulae, PRB 64, 235111 (2001)
                    vXe(issh,ipoint) =                                      &
     &                 v0(jssh)*exp( -1.0_dp*(rcutoff(issh) - r0(jssh))      &
     &                 / (r(ipoint)-r0(jssh)) ) /(rcutoff(issh)+0.01-r(ipoint))
                 else
                    vXe(issh,ipoint) = 0.0_dp
                 endif
              end do ! enddo ipoint
           end do ! enddo  issh
        endif ! if(ioptim)


! Construct approximate starting potential. We do this by forming hydrogenic
! wavefunctions using an effective bohr radius a0_eff(m).
        sigma = 0.0_dp
        do issh = 1, nssh
         l = lam(issh)
         do ipoint = 1, mesh - 1
          if (r(ipoint) .lt. rcutoff_psirc(issh)) then
           psi0(ipoint) = r(ipoint)**(l+1)*exp(-r(ipoint)/a0_eff(issh))
          else
           psi0(ipoint) = 0.0_dp
          end if
         end do
         psi0(mesh) = 0.0_dp

! Grossly enforce the boundary condition at rc.
! Normalize and add to total charge density.
         xnorm = sqrt(dot_product(psi0,psi0)*dr)
         psi(issh,:) = psi0/xnorm
         sigma = sigma + xocc_ion(issh)*psi(issh,:)**2
        end do

! Now find rho dependent potentials, vee and uxc.
! The subroutine get_uxc depends on the actual density rho, calculate rho
! from sigma - rho = sigma/r**2. Also, since we are doing GGA's for some
! options, calculate the first and second derivatives of rho.
        rho(2:mesh-1) = sigma(2:mesh-1)/r(2:mesh-1)**2

! endpoints
        rho(1) = 2.0_dp*rho(2) - rho(3)
        rho(mesh) = 2.0_dp*rho(mesh-1) - rho(mesh-2)

! derivatives
        rhop(2:mesh-1,1) = (rho(3:mesh) - rho(1:mesh-2))/(2.0_dp*dr)
        rhopp(3:mesh-1,1,1) = (rho(4:mesh) - 2.0_dp*rho(3:mesh-1)                 &
     &                                 + rho(2:mesh-2))/(dr**2)

! endpoints
        rhop(1,1) = 2.0_dp*rhop(2,1) - rhop(3,1)
        rhopp(1,1,1) = 2.0_dp*rhopp(2,1,1) - rhopp(3,1,1)
        rhop(mesh,1) = 2.0_dp*rhop(mesh-1,1) - rhop(mesh-2,1)
        rhopp(mesh,1,1) = 2.0_dp*rhopp(mesh-1,1,1) - rhopp(mesh-2,1,1)

        ienergy = 0     ! do not calculate energies
          do ipoint = 1, mesh
            call xc_calc(rho(ipoint)/abohr3, rhopp(ipoint,1,1)/abohr5, rhop(ipoint,:)/abohr4, rhopp(ipoint,:,:)/abohr5, &
              &          eexc(ipoint), uxc(ipoint), dexcrho, dexcsigma, duxcrho, duxcsigma, duxclapl, duxccross)
            eexc(ipoint) = eexc(ipoint)/ryd
            uxc(ipoint) = uxc(ipoint)/ryd
          end do
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

        eold = 0.0_dp
        do iteration = 1, max_iterations  ! ::::::: SCF LOOP ::::::::
         sigma_old = sigma
         sigma = 0.0_dp

         do issh = 1, nssh  ! -------- LOOP OVER ALL SHELLS ---------
          l = lam(issh)

! put potential in v
! Note that the factor of 2 is conversion from Hartree units to Rydberg units.
!            v = 2.0_dp*(vc + vl(issh,:)) + vee + uxc
            v = 2.0_dp*(vc + vl(issh,:)) + vee + uxc + vXe(issh,:)

! Solve Schroedinger equation for potential v
          psi0 = psi(issh,:)
          call psirc (mesh_psirc(issh), 0, l, rcutoff_psirc(issh), v,        &
      &               eout, psi0)
          eigenvalue(issh) = eout

! Add to charge density
          psi(issh,:) = psi0
          sigma = sigma + xocc_ion(issh)*psi(issh,:)**2

          do ipoint = 1, mesh_psirc(issh)
           Rall(issh,ipoint) = 0.0_dp
           if (r(ipoint) .gt. 1.0e-4_dp) Rall(issh,ipoint) = psi0(ipoint)/r(ipoint)
          end do
         end do  ! ---------------- END LOOP OVER SHELLS -----------------

! Mix in old and new sigma's
         sigma = (1.0_dp - beta)*sigma + beta*sigma_old

! Finished with all states of angular momentum l.
! Calculate exchange correlation and Hartree potential for new sigma
! Now find rho dependent potentials, vee and uxc.
! The subroutine get_uxc depends on the actual density rho, calculate rho
! from sigma - rho = sigma/r**2. Also, since we are doing GGA's for some
! options, calculate the first and second derivatives of rho.
         rho(2:mesh-1) = sigma(2:mesh-1)/(4.0_dp*pi*r(2:mesh-1)**2)

! endpoints
         rho(1) = 2.0_dp*rho(2) - rho(3)
         rho(mesh) = 2.0_dp*rho(mesh-1) - rho(mesh-2)

! derivatives
         rhop(2:mesh-1,1) = (rho(3:mesh) - rho(1:mesh-2))/(2.0_dp*dr)
         rhopp(3:mesh-1,1,1) = (rho(4:mesh) - 2.0_dp*rho(3:mesh-1)                &
     &                                  + rho(2:mesh-2))/(dr**2)

! endpoints
         rhop(1,1) = 2.0_dp*rhop(2,1) - rhop(3,1)
         rhopp(1,1,1) = 2.0_dp*rhopp(2,1,1) - rhopp(3,1,1)
         rhop(mesh,1) = 2.0_dp*rhop(mesh-1,1) - rhop(mesh-2,1)
         rhopp(mesh,1,1) = 2.0_dp*rhopp(mesh-1,1,1) - rhopp(mesh-2,1,1)

         ienergy = 1     ! calculate energies
           do ipoint = 1, mesh
             call xc_calc(rho(ipoint)/abohr3, rhopp(ipoint,1,1)/abohr5, rhop(ipoint,:)/abohr4, rhopp(ipoint,:,:)/abohr5, &
               &          eexc(ipoint), uxc(ipoint), dexcrho, dexcsigma, duxcrho, duxcsigma, duxclapl, duxccross)
             eexc(ipoint) = eexc(ipoint)/ryd
             uxc(ipoint) = uxc(ipoint)/ryd
           end do
           exc = 0.0_dp
           do ipoint = 1, mesh
             exc = exc + dr*r(ipoint)**2*rho(ipoint)*eexc(ipoint)
           end do
         call get_vee (mesh, dr, sigma, vee)

! The xc energy is in units of 4*pi, convert back.
         exc = 4.0_dp*pi*exc

! Calculate the total energy:  Etot=sum(eigenval)-1/2*Vee-1/4*Vxc
! (Minus sign to convert binding energies into real energies)
! Sum over the eigenvalues:
         ebs = dot_product(xocc_ion,eigenvalue)

! Evaluate the double counted electron-electron repulsion and the exchange-
! correlation energy.
! Do the integral using a trapezoidal rule.
         eee = 0.0_dp

! eee = vee/2.
! exc = integral dr*sigma*(vxc-uxc)
! Contribution from origin is zero. Find terms in the sum of eigenvalues
         v_local = 0.0_dp
         uxc_local = 0.0_dp
         vna_local = 0.0_dp
         vkin_local = 0.0_dp
         do ipoint = 2, mesh - 1
          eee = eee + 0.5_dp*dr*sigma(ipoint)*vee(ipoint)

! Remember that vc and vl are in hartrees so to get rydbergs multiply by 2.0
          v_local = v_local + 2.0_dp*dr*psi(:,ipoint)**2*vl(:,ipoint)
          uxc_local = uxc_local + dr*psi(:,ipoint)**2*uxc(ipoint)
          vna_local = vna_local + dr*psi(:,ipoint)**2                        &
     &                              *(2.0_dp*vc(ipoint) + vee(ipoint))

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
! Note: In our units (Rydbergs) hbar**2/(2*mass) = 1._dp
! It is faster to not use function d2u to get the second derivative of u(r)
! if you only want to use the 4-point formula.
          do issh = 1, nssh
           ekin = (psi(issh,ipoint-1) - 2.0_dp*psi(issh,ipoint)                 &
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
         terms = 0.0_dp
         do issh = 1, nssh
          terms(1) = terms(1) + xocc_ion(issh)*v_local(issh)
          terms(2) = terms(2) + xocc_ion(issh)*uxc_local(issh)
          terms(3) = terms(3) + xocc_ion(issh)*vna_local(issh)
          terms(4) = terms(4) + xocc_ion(issh)*vkin_local(issh)
         end do

         enew = abs(etot) + abs(ebs)
         difference = abs((enew - eold)/enew)

         if (difference .le. tolerance_scf) exit
         eold = enew
        end do  ! :::::::::::::::: END OF SCF LOOP :::::::::::::::::::::
! *****************************************************************************

        if (difference .gt. tolerance_scf) then
	 call utils_log(UTILS_LOG_ERROR, "self-consistency not reached (DMOL)")   
         error stop 1
        end if

! ***************************************************************************
!       O R T H O G O N A L I Z E
! ***************************************************************************
! First find the integral of u1(x)*psi0(x) for each shell.
        do issh = 1, nssh
         mesh_new = mesh_psirc(issh)

! Set up integration factors.
         rmult(1) = dr/3.0_dp
         rmult(mesh_new) = dr/3.0_dp
         do ipoint = 2, mesh_new - 1, 2
          rmult(ipoint) = 4.0_dp*dr/3.0_dp
         end do
         do ipoint = 3, mesh_new - 2, 2
          rmult(ipoint) = 2.0_dp*dr/3.0_dp
         end do

         integral = 0.0_dp
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_gs(issh,ipoint)*psi(issh,ipoint)
         end do

! Now orthogonalize by taking psi1_new = psi1_old - sum*psi0
         do ipoint = 1, mesh_new
          psi_new(issh,ipoint) = psi(issh,ipoint) - integral*psi_gs(issh,ipoint)
         end do

! Check normalization
         integral = 0.0_dp
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)**2
         end do
         psi_new(issh,1:mesh_new) = psi_new(issh,1:mesh_new)/sqrt(integral)

         integral = 0.0_dp
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)**2
         end do

! Check orthoganalization
         integral = 0.0_dp
         do ipoint = 1, mesh_new
          integral =                                                         &
           integral + rmult(ipoint)*psi_new(issh,ipoint)*psi_gs(issh,ipoint)
         end do
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

! Find (approximately) the value at r=0
         xx(2) = r(2)*abohr
         xx(3) = r(3)*abohr

         yy(2) = psi_new(issh,2)/(abohr1_5*r(2))
         yy(3) = psi_new(issh,3)/(abohr1_5*r(3))

         slope = (yy(3) - yy(2))/(xx(3) - xx(2))
         xx(1) = 0.0_dp
         yy(1) = - slope*xx(2) + yy(2)

         xx(4:mesh_new-1) = r(4:mesh_new-1)*abohr
         yy(4:mesh_new-1) = psi_new(issh,4:mesh_new-1)/(abohr1_5*r(4:mesh_new-1))

         xx(mesh_new) = r(mesh_new)*abohr
         yy(mesh_new) = 0.0_dp

	 if (sav(issh)) then
	   io = utils_open(trim(outpath)//trim(utils_itos(nznuc, lpad=3))//"_"//lnames(l)//".ewf", "w")
	   write (io, "('# Z, l = ',i3,2x,i2)") nznuc, l
	   write (io, "('# q = ',f14.6)") 0.0_dp
	   write (io, "('# n = ',i8)") mesh_new
	   do jpoint = 1, mesh_new
	     write (io, "(2ES16.8)") xx(jpoint), yy(jpoint)
	   end do
	   close(io)
	 end if
        end do
        deallocate (xx)
        deallocate (yy)

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

        return
        end
