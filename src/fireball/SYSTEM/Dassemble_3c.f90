subroutine Dassemble_3c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_constants, only: eq2
  use M_system, only: natoms, ratom, imass, neigh_comb, neigh_comj, neigh_comm, neigh_comn, numorb_max, rho, xl, f3naa, f3nab, f3nac
  use M_fdata, only: num_orb
  implicit none
  integer iatom
  integer ibeta
  integer imu
  integer in1
  integer in2
  integer indna
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer ix
  integer jatom
  integer jbeta
  integer mneigh
  integer ialp
  real(double) cost
  real(double) x
  real(double) y
  real(double), dimension (numorb_max, numorb_max) :: bcnax
  real(double), dimension (3, 3, 3) :: depsA
  real(double), dimension (3, 3, 3) :: depsB
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3, numorb_max, numorb_max) :: f3naXa
  real(double), dimension (3, numorb_max, numorb_max) :: f3naXb
  real(double), dimension (3, numorb_max, numorb_max) :: f3naXc
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r21
  real(double), dimension (3) :: rhat
  real(double), dimension (3) :: rna
  real(double), dimension (3) :: rnabc
  real(double), dimension (3) :: sighat
  f3naa = 0.0d0
  f3nab = 0.0d0
  f3nac = 0.0d0
  do ialp = 1, natoms
   rna(:) = ratom(:,ialp)
   indna = imass(ialp)
   do ineigh = 1, neigh_comn(ialp)
    mneigh = neigh_comm(ineigh,ialp)
    if (mneigh .ne. 0) then
     iatom = neigh_comj(1,ineigh,ialp)
     ibeta = neigh_comb(1,ineigh,ialp)
     r1(:) = ratom(:,iatom) + xl(:,ibeta)
     in1 = imass(iatom)
     jatom = neigh_comj(2,ineigh,ialp)
     jbeta = neigh_comb(2,ineigh,ialp)
     r2(:) = ratom(:,jatom) + xl(:,jbeta)
     in2 = imass(jatom)
     r21(:) = r2(:) - r1(:)
     y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
     if (y .lt. 1.0d-05) then
      sighat(1) = 0.0d0
      sighat(2) = 0.0d0
      sighat(3) = 1.0d0
     else
      sighat(:) = r21(:)/y
     end if
     rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
     x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
     if (x .lt. 1.0d-03) then
      rhat(1) = 0.0d0
      rhat(2) = 0.0d0
      rhat(3) = 0.0d0
     else
      rhat(:) = rnabc(:)/x
     end if
     cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
     call epsilon (rhat, sighat, eps)
     call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
     isorp = 0
     interaction = 1
     call Dtrescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, bcnax, f3naXa, f3naXb, f3naXc) 
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        f3naa(ix,ialp) = f3naa(ix,ialp) + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXa(ix,imu,inu)*eq2
        f3nab(ix,iatom) = f3nab(ix,iatom) + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXb(ix,imu,inu)*eq2
        f3nac(ix,jatom) = f3nac(ix,jatom) + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXc(ix,imu,inu)*eq2
       end do
      end do
     end do
    end if
   end do
  end do
end subroutine Dassemble_3c
