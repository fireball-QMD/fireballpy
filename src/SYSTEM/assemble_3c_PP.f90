subroutine assemble_3c_PP ()
  use iso_c_binding
  use M_system
  use M_fdata, only: num_orb, num_orbPP 
  implicit none
  integer(c_long) ialp
  integer(c_long) iatom
  integer(c_long) ibeta
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) indna
  integer(c_long) ineigh
  integer(c_long) inu
  integer(c_long) jatom
  integer(c_long) jbeta
  integer(c_long) m31
  integer(c_long) m32
  integer(c_long) mneigh
  integer(c_long) ncc
 
  integer(c_long), external :: mpairnay

  real(c_double), dimension (numorb_max, numorb_max) :: bcnlx
  real(c_double), dimension (numorb_max) :: cl
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r31
  real(c_double), dimension (3) :: r32
  real(c_double), dimension (3) :: rna
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
