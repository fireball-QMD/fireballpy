subroutine assemble_3c_PP ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: natoms, ratom, imass, numorb_max, sVNL, vnl, neighPP_comn, neighPP_comm, neighPP_comj, neighPP_comb, xl
  use M_fdata, only: num_orb, num_orbPP 
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
  integer jatom
  integer jbeta
  integer m31
  integer m32
  integer mneigh
  integer ncc
 
  integer, external :: mpairnay

  real(double), dimension (numorb_max, numorb_max) :: bcnlx
  real(double), dimension (numorb_max) :: cl
  real(double), dimension (3) :: r1
  real(double), dimension (3) :: r2
  real(double), dimension (3) :: r31
  real(double), dimension (3) :: r32
  real(double), dimension (3) :: rna
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
            bcnlx(imu,inu) = 0.0d0
            do ncc = 1, num_orbPP(indna)
              bcnlx(imu,inu) = bcnlx(imu,inu) + cl(ncc)*sVNL(imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
            end do ! do ncc
          end do ! do imu
        end do ! do inu
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vnl(imu,inu,mneigh,iatom) =  vnl(imu,inu,mneigh,iatom) + bcnlx(imu,inu)
          end do ! do imu
        end do ! do inu
      end if
    end do ! do ineigh
  end do ! do ialp
  return
end subroutine assemble_3c_PP
