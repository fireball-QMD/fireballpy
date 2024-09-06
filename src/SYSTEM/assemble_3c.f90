subroutine assemble_3c ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system
  use M_fdata, only: num_orb
  implicit none

  integer(c_long) ialp
  integer(c_long) iatom
  integer(c_long) ibeta
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) indna
  integer(c_long) ineigh
  integer(c_long) interaction
  integer(c_long) inu
  integer(c_long) isorp
  integer(c_long) jatom
  integer(c_long) jbeta
  integer(c_long) mneigh
  integer(c_long) jneigh       

  real(c_double) cost
  real(c_double) distance_13
  real(c_double) distance_23
  real(c_double) x
  real(c_double) y

  real(c_double), dimension (numorb_max, numorb_max) :: bcnax
  real(c_double), dimension (3, 3) :: eps
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r21
  real(c_double), dimension (3) :: r13
  real(c_double), dimension (3) :: r23
  real(c_double), dimension (3) :: rhat
  real(c_double), dimension (3) :: rna
  real(c_double), dimension (3) :: rnabc
  real(c_double), dimension (3) :: sighat

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
        jneigh = neigh_back(iatom,mneigh)
        r21(:) = r2(:) - r1(:)
        y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3)) 
        r13(:) = r1(:) - rna(:)
        r23(:) = r2(:) - rna(:)
        distance_13 = sqrt(r13(1)**2 + r13(2)**2 + r13(3)**2)
        distance_23 = sqrt(r23(1)**2 + r23(2)**2 + r23(3)**2)
        if (y .lt. 1.0d-05) then
          sighat(1) = 0.0d0
          sighat(2) = 0.0d0
          sighat(3) = 1.0d0
        else
          sighat(:) = r21(:)/y
        end if
        rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
        x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
        if (x .lt. 1.0d-05) then
          rhat(1) = 0.0d0
          rhat(2) = 0.0d0
          rhat(3) = 0.0d0
        else
          rhat(:) = rnabc(:)/x
        end if
        cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
        call epsilon (rhat, sighat, eps)
        ! CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
        isorp = 0
        interaction = 1
        call trescentros (interaction, isorp, in1, in2, indna, x, y, cost, eps, bcnax)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vna(imu,inu,mneigh,iatom) = vna(imu,inu,mneigh,iatom) + bcnax(imu,inu)*eq2
            vna(inu,imu,jneigh,jatom) = vna(imu,inu,mneigh,iatom)
          end do
        end do
      end if
    end do ! do ineigh
  end do ! do ialp
  return
end subroutine assemble_3c
