subroutine dipole_proyection()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: natoms, ratom, imass, neigh_max, neigh_b, neigh_j, neighn, neigh_self, dq_DP, rho, s_mat, dipc, xl, &
    & Q0_TOT, dipTot_x, dipTot_y, dipTot_z
  use M_fdata, only: num_orb
  implicit none
  real(dp), parameter ::  klambda = 1.4_dp
  integer iatom
  integer imu
  integer in1, in2
  integer ineigh
  integer inu
  integer jatom
  integer mbeta
  real(dp) dip_x, dip_y, dip_z, dip_proy
  real(dp), dimension(3) :: r1,r2,u21,Rbc
  real(dp), dimension(3) :: rmedio, raux
  real(dp)      w_suma
  real(dp), dimension(3,3) :: bwrr, bwrr_inv, u_bwrr, ut_bwrr, v_bwrr, vt_bwrr, zero_bwrr
  real(dp), dimension (natoms) :: c_k
  real(dp), dimension (neigh_max) :: w_k
  real(dp), dimension(3,natoms) :: intra_dip, res_dip
  real(dp), dimension(3,1) :: intra_dip_aux, delta_ck
  integer :: n_bwrr = 3
  integer :: lda_bwrr = 3
  integer :: lwork = 15 !MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) 
  integer i
  integer :: info
  real(dp), dimension(3) :: s_bwrr
  real(dp), dimension(15) :: dummy
  dip_x = 0.0_dp
  dip_y = 0.0_dp
  dip_z = 0.0_dp
  do iatom = 1, natoms
    r1(:) = ratom(:,iatom)
    dip_x = dip_x-Q0_TOT(iatom)*r1(1)
    dip_y = dip_y-Q0_TOT(iatom)*r1(2)
    dip_z = dip_z-Q0_TOT(iatom)*r1(3)
  end do
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      Rbc(:)=(r1(:)+r2(:))/2.0_dp
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
        end do
      end do
    end do
  end do
  dipTot_x = dip_x
  dipTot_y = dip_y
  dipTot_z = dip_z
  do iatom = 1, natoms
    dip_x = 0.0_dp
    dip_y = 0.0_dp
    dip_z = 0.0_dp
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)    
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          if ((iatom .ne. jatom) .or. (mbeta .ne. 0)) then
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1) + dipc(2,imu,inu,ineigh,iatom)*u21(2) + dipc(3,imu,inu,ineigh,iatom)*u21(3)
            dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) - dip_proy*u21(1))
            dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) - dip_proy*u21(2))
            dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) - dip_proy*u21(3))
          end if
        end do
      end do
    end do
    res_dip(1,iatom) = dip_x
    res_dip(2,iatom) = dip_y
    res_dip(3,iatom) = dip_z
  end do
  do iatom = 1, natoms
    dip_x = 0.0_dp
    dip_y = 0.0_dp
    dip_z = 0.0_dp
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    ineigh = neigh_self(iatom)
    do imu = 1,num_orb(in1)
      do inu = 1,num_orb(in1)
        if (imu .ne. inu) then
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + r1(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + r1(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + r1(3)*s_mat(imu,inu,ineigh,iatom))
        end if
      end do
    end do
    intra_dip(1,iatom) = dip_x + res_dip(1,iatom)
    intra_dip(2,iatom) = dip_y + res_dip(2,iatom)
    intra_dip(3,iatom) = dip_z + res_dip(3,iatom)
  end do
  dq_DP = 0.0_dp
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    rmedio = 0.0_dp
    w_suma = 0.0_dp
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      w_k(ineigh)=exp(-klambda*sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      rmedio = rmedio + w_k(ineigh)*r2
      w_suma = w_suma + w_k(ineigh)
    end do !end ineigh = 1,natoms
    rmedio = rmedio / w_suma
    bwrr = 0.0_dp
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      do imu = 1,3 !xyz
        do inu = 1,3 !xyz
          bwrr(inu,imu) = bwrr(inu,imu) + w_k(ineigh)*r2(imu)*(r2(inu)-rmedio(inu))
        enddo !inu
      enddo !imu
    end do !end ineigh = 1,natoms
    !inversa de bwrr
    bwrr_inv = bwrr
    call dgesvd('A', 'S', n_bwrr, n_bwrr, bwrr_inv, lda_bwrr, s_bwrr, u_bwrr, lda_bwrr, vt_bwrr, lda_bwrr, dummy, lwork, info)
    zero_bwrr = 0.00
    do i = 1,3
      if (abs(s_bwrr(i)) .gt. 0.000001) then
        zero_bwrr(i,i) = 1/s_bwrr(i)
      endif
    enddo
    v_bwrr = transpose(vt_bwrr)
    ut_bwrr = transpose(u_bwrr)
    bwrr_inv=matmul(v_bwrr,zero_bwrr)
    bwrr_inv=matmul(bwrr_inv,ut_bwrr)
    intra_dip_aux(1,1)=intra_dip(1,iatom)
    intra_dip_aux(2,1)=intra_dip(2,iatom)
    intra_dip_aux(3,1)=intra_dip(3,iatom)
    delta_ck = matmul(bwrr_inv,intra_dip_aux)
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      raux = r2(:) - rmedio(:)
      c_k(jatom) = w_k(ineigh) * (raux(1)*delta_ck(1,1)+raux(2)*delta_ck(2,1)+raux(3)*delta_ck(3,1))
      dq_DP(jatom) = dq_DP(jatom) + c_k(jatom)
    end do
  end do ! end do iatom = 1,natoms
end subroutine dipole_proyection

function get_min_parabola(x1,x2,x3,y1,y2,y3)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
  real(dp) a,b,c,get_min_parabola
  b=(y2-y3)-(y1-y2)*(x2**2-x3**2)/(x1**2-x2**2)
  b=b/( (x2-x3)-(x1-x2)*(x2**2-x3**2)/(x1**2-x2**2) )
  a=(y1-y2)/(x1**2-x2**2)-b*(x1-x2)/(x1**2-x2**2)
  c=y1-a*x1**2-b*x1
  get_min_parabola=-b/(2*a)
  ! print*,a,'*x**2+',b,'*x+',c
end function get_min_parabola
