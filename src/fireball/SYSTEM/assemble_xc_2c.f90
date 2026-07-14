subroutine assemble_xc_2c ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, natoms, neigh_self, ratom, imass, neighn, neigh_b, neigh_j, xl, vxc, rho_off, rhoij_off, s_mat, Kscf, g_xc, get_shell_ofatom_issh, Qin, iqout, &
    & vxc_ca0, iscf_fast, nssh_tot, get_issh_ofshell, get_iatom_ofshell
  use M_fdata, only: num_orb, nssh, nsh_max, Qneutral, TWOCENTER_VXC_0, TWOCENTER_VXC_L, TWOCENTER_VXC_R
  implicit none
  integer iatom
  integer iatomstart
  integer imu
  integer in1, in2, in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer jatom
  integer kforce
  integer matom
  integer mbeta
  integer natomsp
  integer issh
  real(double)  y 
  real(double) dxn
  real(double), dimension (numorb_max, numorb_max) :: bcxcx
  real(double), dimension (numorb_max, numorb_max) :: denmx
  real(double), dimension (numorb_max, numorb_max) :: den1x
  real(double), dimension (numorb_max, numorb_max) :: rhomx
  real(double), dimension (3, numorb_max, numorb_max) :: rhompx
  real(double), dimension (3, 3) :: eps
  real(double), dimension (3, 3, 3) :: deps
  real(double), dimension (3) :: r1, r2, r21
  real(double), dimension (3) :: sighat
  real(double), dimension (numorb_max, numorb_max) :: sx
  real(double), dimension (nsh_max) :: dqi
  real(double), dimension (nsh_max) :: dqj
  ! Speed up SCF loop: en Kscf>1 el termino de cargas estaticas se reconstruye como
  ! vxc_ca0 + g_xc.dQ (exactamente lineal en dQ), sin llamadas a doscentros.
  logical :: fast
  integer :: alpha
  real(double), dimension (nssh_tot) :: dQg

  fast = (Kscf .gt. 1) .and. (iscf_fast .eq. 1)
  if (fast) then
    do alpha = 1, nssh_tot
      dQg(alpha) = Qin(get_issh_ofshell(alpha), get_iatom_ofshell(alpha)) &
        & - Qneutral(get_issh_ofshell(alpha), imass(get_iatom_ofshell(alpha)))
    end do
  end if
  if (Kscf .eq. 1) vxc_ca0 = 0.0d0

  kforce = 0
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    dqi = 0.0d0
    do issh = 1, nssh(in1)
      dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
    end do

    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
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
      dqj = 0.0d0
      do issh = 1, nssh(in2)
        dqj(issh) = (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      if (.not. fast) then
        call epsilon (r2, sighat, eps)
        call deps2cent (r1, r2, eps, deps)
      end if
      if (iatom .ne. jatom .or. mbeta .ne. 0) then
        if (fast) then
          ! termino de cargas estaticas congelado: VXC_0 (vxc_ca0) + suma_alpha g_xc.dQ
          in3 = in2
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + vxc_ca0(imu,inu,ineigh,iatom) &
                & + dot_product(g_xc(:,imu,inu,ineigh,iatom), dQg)
            end do
          end do
        else
        isorp = 0
        interaction = TWOCENTER_VXC_0
        in3 = in2
        call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)
            if (Kscf .eq. 1) vxc_ca0(imu,inu,ineigh,iatom) = rhomx(imu,inu)
          end do
        end do
        interaction = TWOCENTER_VXC_L
        in3 = in2
        do isorp = 1, nssh(in1)
          call doscentros (interaction, isorp, kforce, in1, in1, in3, y, eps, deps, rhomx, rhompx)
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dqi(isorp)
              if (Kscf .eq. 1) then
                g_xc(get_shell_ofatom_issh(iatom,isorp),imu,inu,ineigh,iatom) = g_xc(get_shell_ofatom_issh(iatom,isorp),imu,inu,ineigh,iatom) + rhomx(imu,inu)
              end if
            end do
          end do
        end do
        interaction = TWOCENTER_VXC_R
        in3 = in2
        do isorp = 1, nssh(in2)
          call doscentros (interaction, isorp, kforce, in1, in2, in3, y, eps, deps, rhomx, rhompx)
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dqj(isorp)
              if (Kscf .eq. 1) then
                g_xc(get_shell_ofatom_issh(jatom,isorp),imu,inu,ineigh,iatom) = g_xc(get_shell_ofatom_issh(jatom,isorp),imu,inu,ineigh,iatom) + rhomx(imu,inu)
              end if
            end do
          end do
        end do
        end if   ! end if (fast)
        in3 = in2

        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            denmx(imu,inu) = rho_off(imu,inu,ineigh,iatom)
            den1x(imu,inu) = rhoij_off(imu,inu,ineigh,iatom)
            sx(imu,inu) = s_mat(imu,inu,ineigh,iatom)
          end do
        end do

        ! Calculate <i| V_xc(n) |j> and <i|V_xc(n_i+n_j)|j>        
        call build_olsxc_off (in1, in2, den1x, denmx, sx, ineigh, iatom, bcxcx)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
          end do
        end do

      end if
    end do
  end do
  ! vxc = 0.0d0
  return
end subroutine assemble_xc_2c
