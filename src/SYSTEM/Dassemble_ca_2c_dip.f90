subroutine Dassemble_ca_2c_dip ()
  use iso_c_binding
  use M_constants, only: eq2
  use M_system
  use M_fdata, only: nssh,rcutoff,Qneutral,num_orb
  implicit none
  integer(c_long) iatom
  integer(c_long) icount
  integer(c_long) icount_sav
  integer(c_long) ierror
  integer(c_long) imu
  integer(c_long) in1
  integer(c_long) in2
  integer(c_long) in3
  integer(c_long) ineigh
  integer(c_long) interaction
  integer(c_long) inu
  integer(c_long) isorp
  integer(c_long) issh
  integer(c_long) ix
  integer(c_long) jatom
  integer(c_long) jcount
  integer(c_long) jcount_sav
  integer(c_long) jssh
  integer(c_long) kforce
  integer(c_long) matom
  integer(c_long) mbeta
  real(c_double) dq1
  real(c_double) dq2
  real(c_double) dstn_temp
  real(c_double) dterm_1
  real(c_double) dterm_2
  real(c_double) dterm
  real(c_double) dxn
  real(c_double) rcutoff_j
  real(c_double) rend
  real(c_double) rend1
  real(c_double) rend2
  real(c_double) sterm_1
  real(c_double) sterm_2
  real(c_double) y
  real(c_double) rcutoff_i
  real(c_double), dimension (numorb_max, numorb_max) :: bcca
  real(c_double), dimension (3, numorb_max, numorb_max) :: bccap
  real(c_double), dimension (3, numorb_max, numorb_max) :: bccapx
  real(c_double), dimension (numorb_max, numorb_max) :: bccax
  real(c_double), dimension (3,numorb_max, numorb_max) :: demnpl
  real(c_double), dimension (3, 3, 3) :: deps
  real(c_double), dimension (3, numorb_max, numorb_max) :: dewaldsr
  real(c_double), dimension (3) :: dpterm_1
  real(c_double), dimension (3) :: dpterm_2
  real(c_double) dstn1
  real(c_double) dstn2
  real(c_double), dimension (numorb_max, numorb_max) :: emnpl
  real(c_double), dimension (3, 3) :: eps
  real(c_double), dimension (3) :: r1
  real(c_double), dimension (3) :: r2
  real(c_double), dimension (3) :: r21
  real(c_double), dimension (3) :: sighat
  real(c_double), dimension (3) :: spterm_1
  real(c_double), dimension (3) :: spterm_2
  real(c_double) stn1
  real(c_double) stn2 
  faca = 0.0d0
  fotca = 0.0d0
  do iatom = 1, natoms
   matom = neigh_self(iatom)
   r1(:) = ratom(:,iatom)
   in1 = imass(iatom)
    rcutoff_i = 0
    do imu = 1, nssh(in1)
     if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
    end do
   dq1 = 0.0
   do issh = 1, nssh(in1)
    dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
   end do
   do ineigh = 1, neighn(iatom)    ! <==== loop 2 over iatom's neighbors
    mbeta = neigh_b(ineigh,iatom)
    jatom = neigh_j(ineigh,iatom)
    r2(:) = ratom(:,jatom) + xl(:,mbeta)
    in2 = imass(jatom)
    rcutoff_j = 0
    do imu = 1, nssh(in2)
     if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
    end do
    dq2 = 0.0d0
    do issh = 1, nssh(in2)
     dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
    end do
    r21(:) = r2(:) - r1(:)
    y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
    if (y .lt. 1.0d-05) then
     sighat(1) = 0.0d0
     sighat(2) = 0.0d0
     sighat(3) = 1.0d0
    else
     sighat(:) = r21(:)/y
    end if
    call epsilon (r2, sighat, eps)
    call deps2cent (r1, r2, eps, deps)
    stn1 = 1.0d0
    dstn1 = 0.0d0
    stn2 = 0.0d0
    dstn2 = 0.0d0
    emnpl = 0.0d0
    demnpl = 0.0d0
    dewaldsr = 0.0d0
    if (y .gt. 1.0d-4) then
       rend = rcutoff_i + rcutoff_j
       call smoother (y, rend, smt_elect, stn1, dstn_temp)
       stn2 = 1.0d0 - stn1
       dstn1 = dstn_temp
       dstn2 = - dstn_temp
     do inu = 1, num_orb(in1)
      do imu = 1, num_orb(in1)
       ! ==============  NEW STUFF : COMPUTE THE TRUE DIPOLE
       ! ===========
       dterm = (dipc(1,imu,inu,matom,iatom)*r21(1)     &
          &   + dipc(2,imu,inu,matom,iatom)*r21(2)     &
          &   + dipc(3,imu,inu,matom,iatom)*r21(3))
       emnpl(imu,inu) =  dq2*(s_mat(imu,inu,matom,iatom)/y)   &
       &               + dq2*(dterm/(y*y*y))
       demnpl(:,imu,inu)= (dq2*s_mat(imu,inu,matom,iatom)/(y*y*y))*r21(:)  &
       &                 - dq2*dipc(:,imu,inu,matom,iatom)/(y*y*y)        &
       &                 + dq2*3*dterm*r21(:)/(y*y*y*y*y)
       dewaldsr(:,imu,inu) = eq2*demnpl(:,imu,inu)
      end do !end do imu = 1,num_orb(in1)
     end do  !end do inu = 1,num_orb(in1)
    end if
    bcca = 0.0d0
    bccap = 0.0d0
    kforce = 1
    interaction = 4
    in3 = in1
    do isorp = 1, nssh(in2)
     call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  &
     &                eps, deps, bccax, bccapx)
     dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
     do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
       bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
       bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
      end do
     end do
    end do
    do inu = 1, num_orb(in3)
     do imu = 1, num_orb(in1)
      bccap(:,imu,inu) =                                               &
     &    stn1*bccap(:,imu,inu)                               &
     &  - dstn1*bcca(imu,inu)*sighat(:)                       &
     &  + stn2*demnpl(:,imu,inu)                              &
     &  - dstn2*emnpl(imu,inu)*sighat(:)
     end do
    end do
    do inu = 1, num_orb(in3)
     do imu = 1, num_orb(in1)
      do ix = 1, 3
       faca(ix,ineigh,iatom) = faca(ix,ineigh,iatom)                   &
     &  - rho(imu,inu,matom,iatom)*bccap(ix,imu,inu)*eq2               &
     &  + rho(imu,inu,matom,iatom)*dewaldsr(ix,imu,inu)
      end do
     end do
    end do
    if (iatom .eq. jatom .and. mbeta .eq. 0) then
    else
     bccap = 0.0d0
     interaction = 2
     in3 = in2
     do isorp = 1, nssh(in1)
      call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
     &                 eps, deps, bccax, bccapx)
      dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
      do imu = 1, num_orb(in1)
       do inu = 1, num_orb(in3)
        bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
       end do
      end do
     end do
     do inu = 1, num_orb(in3)
      do imu = 1, num_orb(in1)
       do ix = 1, 3
       fotca(ix,ineigh,iatom) = fotca(ix,ineigh,iatom)                 &
     &  - rho(imu,inu,ineigh,iatom)*bccap(ix,imu,inu)*eq2              
       end do
      end do
     end do
    end if
   end do
  end do
end

