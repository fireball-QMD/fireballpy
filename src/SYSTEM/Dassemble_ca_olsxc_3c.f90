subroutine Dassemble_ca_olsxc_3c ()
  use M_constants, only: wp
  use M_system
  use M_fdata, only: nssh,isorpmax,nspecies,num_orb,lssh,nsh_max
  implicit none
  integer ialp
  integer iatom
  integer ibeta
  integer ierror
  integer imu
  integer inu
  integer in1
  integer in2
  integer indna
  integer index1
  integer index2
  integer ineigh
  integer interaction
  integer isorp
  integer issh
  integer ix
  integer jssh
  integer jatom
  integer jbeta
  integer l1
  integer l2
  integer mneigh
  integer n1
  integer n2
  real(wp) cost
  real(wp) x
  real(wp) y
  real(wp) muxc
  real(wp) dmuxc
  real(wp) d2muxc
  real(wp) exc
  real(wp) dexc
  real(wp) d2exc
  real(wp) sx
  real(wp) sm
  real(wp) rho_av
  real(wp), dimension (3, 3, 3) :: depsA
  real(wp), dimension (3, 3, 3) :: depsB
  real(wp), dimension (3, 3) :: eps
  real(wp), dimension (nsh_max, nsh_max) :: rho_3c
  real(wp), dimension (3, nsh_max, nsh_max) :: rhop_3ca
  real(wp), dimension (3, nsh_max, nsh_max) :: rhop_3cb
  real(wp), dimension (3, nsh_max, nsh_max) :: rhop_3cc
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoxpa
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoxpb
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoxpc
  real(wp), dimension (3, nsh_max, nsh_max) :: rhompa
  real(wp), dimension (3, nsh_max, nsh_max) :: rhompb
  real(wp), dimension (3, nsh_max, nsh_max) :: rhompc
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoinpa
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoinpb
  real(wp), dimension (3, numorb_max, numorb_max) :: rhoinpc
  real(wp), dimension (3, numorb_max, numorb_max) :: avrhop_a
  real(wp), dimension (3, numorb_max, numorb_max) :: avrhop_b
  real(wp), dimension (3, numorb_max, numorb_max) :: avrhop_c
  real(wp), dimension (3, numorb_max, numorb_max) :: mxca
  real(wp), dimension (3, numorb_max, numorb_max) :: mxcb
  real(wp), dimension (3, numorb_max, numorb_max) :: mxcc
  real(wp), dimension (numorb_max, numorb_max) :: rhoin
  real(wp), dimension (nsh_max, nsh_max) :: rhomm
  real(wp), dimension (3) :: spm
  real(wp), dimension (3) :: rhop_avb
  real(wp), dimension (3) :: rhop_avc
  real(wp), dimension (3) :: r1
  real(wp), dimension (3) :: r2
  real(wp), dimension (3) :: r21
  real(wp), dimension (3) :: rhat
  real(wp), dimension (3) :: rna
  real(wp), dimension (3) :: rnabc
  real(wp), dimension (3) :: sighat
  f3xca = 0.0d0
  f3xcb = 0.0d0
  f3xcc = 0.0d0
  f3xca_ca = 0.0d0
  f3xcb_ca = 0.0d0
  f3xcc_ca = 0.0d0
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
          write (*,*) ' There is an error here in assemble_olsxc_off.f '
          write (*,*) ' r1 = r2!!!! BAD!!!! '
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
        interaction = 3
        rho_3c = 0.0d0
        rhop_3ca = 0.0d0
        rhop_3cb = 0.0d0
        rhop_3cc = 0.0d0
        rhoinpa = 0.0d0
        rhoinpb = 0.0d0
        rhoinpc = 0.0d0
        avrhop_a = 0.0d0
        avrhop_b = 0.0d0
        avrhop_c = 0.0d0
        do isorp = 1, nssh(indna)
          call Dtrescentros (interaction, isorp, isorpmax, in1, in2, indna, x, y, cost, eps, depsA, depsB, rhat, sighat, rhoin, rhoxpa, rhoxpb, rhoxpc)
          call DtrescentrosS (isorp, isorpmax, in1, in2, indna, x, y, cost, rhat, sighat, rhomm, rhompa, rhompb, rhompc)
          rhoin(:,:) = rho_off(:,:,mneigh,iatom)
          do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
              rhoinpa(:,imu,inu) =  rhoinpa(:,imu,inu) - rhoxpa(:,imu,inu)*Qin(isorp,ialp)
              rhoinpb(:,imu,inu) =  rhoinpb(:,imu,inu) - rhoxpb(:,imu,inu)*Qin(isorp,ialp)
              rhoinpc(:,imu,inu) =  rhoinpc(:,imu,inu) - rhoxpc(:,imu,inu)*Qin(isorp,ialp)
            end do
          end do
          do inu = 1, nssh(in2)
            do imu = 1, nssh(in1)
              rho_3c(imu,inu) =          rho_3c(imu,inu)  + rhomm(imu,inu)*Qin(isorp,ialp)
              rhop_3ca(:,imu,inu) =  rhop_3ca(:,imu,inu) - rhompa(:,imu,inu)*Qin(isorp,ialp)
              rhop_3cb(:,imu,inu) =  rhop_3cb(:,imu,inu) - rhompb(:,imu,inu)*Qin(isorp,ialp)
              rhop_3cc(:,imu,inu) =  rhop_3cc(:,imu,inu) - rhompc(:,imu,inu)*Qin(isorp,ialp)
            end do
          end do
        end do
        do issh = 1, nssh(in1)
          do jssh = 1, nssh(in2)
            sm = sm_mat(issh,jssh,mneigh,iatom)
            if (abs(sm) .lt. xc_overtol) then
              if(sm .gt. 0.0d0) then
                sm = xc_overtol
              else
                sm = -1.0d0*xc_overtol
              endif
            endif
            spm(:) = spm_mat(:,issh,jssh,mneigh,iatom)
            avrhop_b (:,issh,jssh) = avrhop_b (:,issh,jssh) + (sm*rhop_3cb(:,issh,jssh) - rho_3c(issh,jssh)*spm(:))/(sm*sm)
            avrhop_c (:,issh,jssh) = avrhop_c (:,issh,jssh) + (sm*rhop_3cc(:,issh,jssh) + rho_3c(issh,jssh)*spm(:))/(sm*sm)
          end do
        end do
        n1 = 0
        do issh = 1, nssh(in1)
          l1 = lssh(issh,in1)
          n1 = n1 + l1 + 1
          n2 = 0
          do jssh = 1, nssh(in2)
            l2 = lssh(jssh,in2)
            n2 = n2 + l2 + 1
            rho_av =  arho_off(issh,jssh,mneigh,iatom)
            rhop_avb(:) =  avrhop_b(:,issh,jssh)
            rhop_avc(:) =  avrhop_c(:,issh,jssh)
            call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
            do index1 = -l1, l1
              imu = n1 + index1
              do index2 = -l2, l2
                inu = n2 + index2
                sx = s_mat(imu,inu,mneigh,iatom)
                mxcb(:,imu,inu) = rhop_avb(:)*d2muxc*(rhoin(imu,inu) - rho_av*sx) + rhoinpb(:,imu,inu)*dmuxc
                mxcc(:,imu,inu) = rhop_avc(:)*d2muxc*(rhoin(imu,inu) - rho_av*sx) + rhoinpc(:,imu,inu)*dmuxc
                mxca(:,imu,inu) = - mxcb(:,imu,inu) - mxcc(:,imu,inu)
              end do
            end do
            n2 = n2 + l2
          end do
          n1 = n1 + l1
        end do
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            do ix = 1, 3
              f3xca_ca(ix,ialp)  = f3xca_ca(ix,ialp)    - 2*rho(imu,inu,mneigh,iatom)*mxca(ix,imu,inu)
              f3xcb_ca(ix,iatom) = f3xcb_ca(ix,iatom)  - 2*rho(imu,inu,mneigh,iatom)*mxcb(ix,imu,inu)
              f3xcc_ca(ix,jatom) = f3xcc_ca(ix,jatom)  - 2*rho(imu,inu,mneigh,iatom)*mxcc(ix,imu,inu)
            end do
          end do
        end do
      end if
    end do
  end do
end subroutine Dassemble_ca_olsxc_3c

