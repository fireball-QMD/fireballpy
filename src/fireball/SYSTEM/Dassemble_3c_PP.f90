subroutine Dassemble_3c_PP () 
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, ratom, imass, rhoPP, numorb_max, sVNL, spVNL, neighPP_comn, neighPP_comm, neighPP_comj, neighPP_comb, &
    & xl, f3nla, f3nlb, f3nlc
  use M_fdata, only: num_orb,num_orbPP
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer imu
  integer in1
  integer in2
  integer indna
  integer ineigh
  integer inu
  integer ix
  integer jatom
  integer jbeta
  integer m31
  integer m32
  integer mneigh
  integer ncc
  integer, external :: mpairnay
  real(double), dimension (numorb_max, numorb_max) :: bcnlx
  real(double), dimension (numorb_max) :: cl
  real(double), dimension (3, numorb_max, numorb_max) :: f3nlXa
  real(double), dimension (3, numorb_max, numorb_max) :: f3nlXb
  real(double), dimension (3, numorb_max, numorb_max) :: f3nlXc
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r31
  real(double), dimension (3) :: r32
  real(double), dimension (3) :: rna
  f3nla = 0.0d0
  f3nlb = 0.0d0
  f3nlc = 0.0d0
  do ialp = 1, natoms
   rna(:) = ratom(:,ialp)
   indna = imass(ialp)
   call cl_value (indna, cl)
   do ineigh = 1, neighPP_comn(ialp)
    mneigh = neighPP_comm(ineigh,ialp)
    if (mneigh .ne. 0) then
     iatom = neighPP_comj(1,ineigh,ialp)
     ibeta = neighPP_comb(1,ineigh,ialp)
     r1(:) = ratom(:,iatom) + xl(:,ibeta)
     in1 = imass(iatom)
     jatom = neighPP_comj(2,ineigh,ialp)
     jbeta = neighPP_comb(2,ineigh,ialp)
     r2(:) = ratom(:,jatom) + xl(:,jbeta)
     in2 = imass(jatom)
     r31(:) = rna(:) - r1(:)
     r32(:) = rna(:) - r2(:)
     m31 = mpairnay (iatom, ialp, r31)
     m32 = mpairnay (jatom, ialp, r32)
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       f3nlXb(:,imu,inu) = 0.0d0
       f3nlXc(:,imu,inu) = 0.0d0
       bcnlx(imu,inu) = 0.0d0
       do ncc = 1, num_orbPP(indna)
        do ix = 1, 3
         f3nlXb(ix,imu,inu) = f3nlXb(ix,imu,inu)  - cl(ncc)*spVNL(ix,imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
         f3nlXc(ix,imu,inu) = f3nlXc(ix,imu,inu)  - cl(ncc)*sVNL(imu,ncc,m31,iatom)*spVNL(ix,inu,ncc,m32,jatom)
        end do
        bcnlx(imu,inu) = bcnlx(imu,inu)  + cl(ncc)*sVNL(imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
       end do
       f3nlXa(:,imu,inu) = - f3nlXb(:,imu,inu) - f3nlXc(:,imu,inu)
      end do
     end do
     do inu = 1, num_orb(in2)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
        f3nla(ix,ialp) = f3nla(ix,ialp)  + rhoPP(imu,inu,mneigh,iatom)*f3nlXa(ix,imu,inu)
        f3nlb(ix,iatom) = f3nlb(ix,iatom)+ rhoPP(imu,inu,mneigh,iatom)*f3nlXb(ix,imu,inu)
        f3nlc(ix,jatom) = f3nlc(ix,jatom)+ rhoPP(imu,inu,mneigh,iatom)*f3nlXc(ix,imu,inu)
       end do
      end do
     end do
    end if
   end do
  end do
end subroutine Dassemble_3c_PP
