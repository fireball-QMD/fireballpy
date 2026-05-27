subroutine Dassemble_ca_olsxc_3c()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_system, only: xc_overtol, natoms, ratom, imass, neigh_comb, neigh_comj, neigh_comm, &
  & neigh_comn, neigh_back, numorb_max, Qin, rho, rho_off, rhomp_2c, neigh_com_ng, &
  & s_mat, arho_off, xl, f3xca_ca, f3xcb_ca, f3xcc_ca, f3xca, f3xcb, f3xcc
  use M_fdata, only: nssh, num_orb, lssh, nsh_max
  implicit none
  integer :: ialp, iatom, ibeta, imu, inu, in1, in2, indna, index1, index2, ineigh, interaction, &
  &        isorp, issh, jssh, jatom, jbeta, jneigh, l1, l2, mneigh, n1, n2, iback, jback
  real(dp) :: cost, x, y, muxc, dmuxc, d2muxc, exc, dexc, d2exc, sx, sm, rho_av
  real(dp), dimension(3, 3, 3) :: depsA
  real(dp), dimension(3, 3, 3) :: depsB
  real(dp), dimension(3, 3) :: eps
  real(dp), dimension(3, numorb_max, numorb_max) :: rhoxpa, rhoxpb, rhoxpc
  real(dp), dimension(3, nsh_max, nsh_max) :: rhompa, rhompb, rhompc
  real(dp), dimension(3, numorb_max, numorb_max) :: rhoinpb, rhoinpc
  real(dp), dimension(3, numorb_max, numorb_max) :: avrhop_a, avrhop_b, avrhop_c
  real(dp), dimension(3, numorb_max, numorb_max) :: mxca, mxcb, mxcc
  real(dp), dimension(numorb_max, numorb_max) :: rhoin
  real(dp), dimension(nsh_max, nsh_max) :: rhomm
  real(dp), dimension(3) :: spm, rhop_avb, rhop_avc, r1, r2, r21, rhat, rna, rnabc, sighat
  f3xca = 0.0_dp
  f3xcb = 0.0_dp
  f3xcc = 0.0_dp
  f3xca_ca = 0.0_dp
  f3xcb_ca = 0.0_dp
  f3xcc_ca = 0.0_dp
  do ialp = 1, natoms
    rna(:) = ratom(:, ialp)
    indna = imass(ialp)
    do ineigh = 1, neigh_comn(ialp)
      mneigh = neigh_comm(ineigh, ialp)
      if (mneigh /= 0) then
        iatom = neigh_comj(1, ineigh, ialp)
        ibeta = neigh_comb(1, ineigh, ialp)
        iback = neigh_com_ng(1, ineigh, ialp)
        r1 = ratom(:, iatom) + xl(:, ibeta)
        in1 = imass(iatom)
        jatom = neigh_comj(2, ineigh, ialp)
        jbeta = neigh_comb(2, ineigh, ialp)
        jback = neigh_com_ng(2, ineigh, ialp)
        jneigh = neigh_back(iatom, mneigh)
        r2 = ratom(:, jatom) + xl(:, jbeta)
        in2 = imass(jatom)
        r21 = r2 - r1
        y = norm2(r21)
        if (y < 1.0e-05_dp) then
          sighat(1) = 0.0_dp
          sighat(2) = 0.0_dp
          sighat(3) = 1.0_dp
        else
          sighat = r21/y
        end if
        rnabc = rna - (r1 + 0.5_dp*r21)
        x = norm2(rnabc)
        if (x < 1.0e-03_dp) then
          rhat(1) = 0.0_dp
          rhat(2) = 0.0_dp
          rhat(3) = 0.0_dp
        else
          rhat = rnabc/x
        end if
        cost = dot_product(sighat, rhat)
        call epsilon(rhat, sighat, eps)
        call deps3center(r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
        interaction = 3
        rhoinpb = 0.0_dp
        rhoinpc = 0.0_dp
        avrhop_a = 0.0_dp
        avrhop_b = 0.0_dp
        avrhop_c = 0.0_dp
        do isorp = 1, nssh(indna)
          call Dtrescentros(interaction, isorp, in1, in2, indna, x, y, cost, eps, &
          &               depsA, depsB, rhat, sighat, rhoin, rhoxpa, rhoxpb, rhoxpc)
          call DtrescentrosS(isorp, in1, in2, indna, x, y, cost, rhat, sighat, rhomm, rhompa, rhompb, rhompc)
          rhoin = rho_off(:, :, mneigh, iatom)
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              rhoinpb(:, imu, inu) = rhoinpb(:, imu, inu) - rhoxpb(:, imu, inu)*Qin(isorp, ialp)
              rhoinpc(:, imu, inu) = rhoinpc(:, imu, inu) - rhoxpc(:, imu, inu)*Qin(isorp, ialp)
            end do
          end do
        end do
        do jssh = 1, nssh(in2)
          do issh = 1, nssh(in1)
            ! Minus as iback is the neighbour of ialp corresponding to iatom
            avrhop_b(:, issh, jssh) = avrhop_b(:, issh, jssh) - 0.5_dp*rhomp_2c(:, issh, iback, iatom)
            avrhop_c(:, issh, jssh) = avrhop_c(:, issh, jssh) - 0.5_dp*rhomp_2c(:, jssh, jback, jatom)
          end do
        end do
        n1 = 0
        do issh = 1, nssh(in1)
          l1 = lssh(issh, in1)
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, nssh(in2)
            l2 = lssh(jssh, in2)
            n2 = n2 + l2 + 1
            rho_av = arho_off(issh, jssh, mneigh, iatom)
            rhop_avb = avrhop_b(:, issh, jssh)
            rhop_avc = avrhop_c(:, issh, jssh)
            call cepal(rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
            do index1 = -l1, l1
              imu = n1 + index1
              do index2 = -l2, l2
                inu = n2 + index2
                sx = s_mat(imu, inu, mneigh, iatom)
                mxcb(:, imu, inu) = rhop_avb*d2muxc*(rhoin(imu, inu) - rho_av*sx) + rhoinpb(:, imu, inu)*dmuxc
                mxcc(:, imu, inu) = rhop_avc*d2muxc*(rhoin(imu, inu) - rho_av*sx) + rhoinpc(:, imu, inu)*dmuxc
                mxca(:, imu, inu) = -mxcb(:, imu, inu) - mxcc(:, imu, inu)
              end do
            end do
            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            f3xca_ca(:, ialp) = f3xca_ca(:, ialp) - 2*rho(imu, inu, mneigh, iatom)*mxca(:, imu, inu)
            f3xcb_ca(:, iatom) = f3xcb_ca(:, iatom) - 2*rho(imu, inu, mneigh, iatom)*mxcb(:, imu, inu)
            f3xcc_ca(:, jatom) = f3xcc_ca(:, jatom) - 2*rho(imu, inu, mneigh, iatom)*mxcc(:, imu, inu)
          end do
        end do
      end if
    end do
  end do
end subroutine Dassemble_ca_olsxc_3c
