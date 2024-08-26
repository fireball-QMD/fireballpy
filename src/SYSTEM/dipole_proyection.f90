subroutine dipole_proyection()
  use M_constants, only: wp
  use M_system
  use M_fdata, only: num_orb,nssh,lssh,Qneutral
  implicit none
  real(wp), parameter ::  Debye = 0.208194
  real(wp), parameter ::  klambda = 1.0
  integer iatom
  integer imu
  integer in1, in2
  integer ineigh
  integer inu
  integer issh
  integer jatom
  integer mbeta
  real(wp)      Qtot, Qtot1, Qtot2
  real(wp), dimension(3) :: r1,r2,Rbc,u21
  real(wp), dimension(3) :: rmedio, raux
  real(wp)      w_suma  
  real(wp), dimension(3,3) :: bwrr, bwrr_inv, u_bwrr, ut_bwrr, v_bwrr, vt_bwrr, zero_bwrr
  real(wp), dimension (natoms) :: c_k
  real(wp), dimension (neigh_max) :: w_k
  real(wp), dimension(3,natoms) :: intra_dip, res_dip
  real(wp), dimension(3,1) :: intra_dip_aux, delta_ck
  integer :: n_bwrr = 3
  integer :: lda_bwrr = 3
  integer :: lwork = 15 !MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) 
  integer info,i
  integer, dimension(3) :: ipiv
  real(wp), dimension(3) :: s_bwrr
  real(wp), dimension(15) :: dummy
  do iatom = 1, natoms
    Q0_TOT(iatom) = 0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
    end do
  end do
  dip_x=0.0d0
  dip_y=0.0d0
  dip_z=0.0d0
  do iatom = 1, natoms
    Qtot=-Q0_TOT(iatom)
    in1 = imass(iatom)
    do issh = 1,nssh(in1)
      Qtot = Qtot+Qout(issh,iatom)
    end do
    dip_x = dip_x+Qtot*ratom(1,iatom)
    dip_y = dip_y+Qtot*ratom(2,iatom)
    dip_z = dip_z+Qtot*ratom(3,iatom)      
  enddo !end do iatom = 1,natoms
  dip_tot      = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )    
  dipQin_x    = dip_x
  dipQin_y    = dip_y
  dipQin_z    = dip_z
  dipQin_tot = dip_tot
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do issh = 1,nssh(in1)
     Qtot = Qtot+Qin(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    dip_x = dip_x-Q0_TOT(iatom)*r1(1)
    dip_y = dip_y-Q0_TOT(iatom)*r1(2)    
    dip_z = dip_z-Q0_TOT(iatom)*r1(3)
  end do !end do iatom = 1,natoms
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)    
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)      
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
  end do ! end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  dipTot_x    = dip_x
  dipTot_y    = dip_y
  dipTot_z    = dip_z
  dipTot_tot = dip_tot
  ! write the intraatomic dipole
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0      
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          if ((iatom .eq. jatom) .and. (imu .ne. inu)) then
            dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
            dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
            dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
          end if !end if
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
  end do ! end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  dipIntra_x    = dip_x
  dipIntra_y    = dip_y
  dipIntra_z    = dip_z
  dipIntra_tot = dip_tot
  !DIP RES = DIP_TOT - PROY
  dip_res_x = 0.0d0
  dip_res_y = 0.0d0
  dip_res_z = 0.0d0
  do iatom = 1, natoms
    dip_res_x = 0.0d0
    dip_res_y = 0.0d0
    dip_res_z = 0.0d0
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)    
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      Rbc(:)=(r1(:)+r2(:))/2.0d0
      do imu = 1,num_orb(in1)
        do inu = 1,num_orb(in2)
          if (iatom .ne. jatom) then !.or. imu .eq. inu) then
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1) + dipc(2,imu,inu,ineigh,iatom)*u21(2) + dipc(3,imu,inu,ineigh,iatom)*u21(3)
            dip_res_x = dip_res_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) - dip_proy*u21(1))
            dip_res_y = dip_res_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) - dip_proy*u21(2))
            dip_res_z = dip_res_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) - dip_proy*u21(3))
          end if !end if
        end do !end do inu
      end do !end do imu
    end do !end ineigh = 1,natoms
    dip_res_tot = sqrt (dip_res_x**2 + dip_res_y**2 + dip_res_z**2 )
    dipRes_x = dipRes_x + dip_res_x
    dipRes_y = dipRes_y + dip_res_y
    dipRes_z = dipRes_z + dip_res_z
    res_dip(1,iatom) = dip_res_x
    res_dip(2,iatom) = dip_res_y
    res_dip(3,iatom) = dip_res_z
  end do ! end do iatom = 1,natoms
  dipRes_tot = sqrt (dipRes_x**2 + dipRes_y**2 + dipRes_z**2 ) 
  do iatom = 1, natoms
    dip_x = 0.0d0
    dip_y = 0.0d0
    dip_z = 0.0d0
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    Qtot=0.0d0
    Qtot1=0.0d0
    do issh = 1,nssh(in1)
      Qtot1 = Qtot1+Qin(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    Qtot=0.0d0
    Qtot2=0.0d0
    do issh = 1,nssh(in1)
      Qtot2 = Qtot2+Qout(issh,iatom)
    end do !end do issh = 1,nssh(in1)
    jatom=iatom
    ineigh=neigh_self(iatom)
    in2=in1
    r2(:)=r1(:)
    Rbc(:)=(r1(:)+r2(:))/2.0d0
    do imu = 1,num_orb(in1)
      do inu = 1,num_orb(in2)
        if ( (imu .ne. inu)) then
          dip_x = dip_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom) + Rbc(1)*s_mat(imu,inu,ineigh,iatom))
          dip_y = dip_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom) + Rbc(2)*s_mat(imu,inu,ineigh,iatom))
          dip_z = dip_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom) + Rbc(3)*s_mat(imu,inu,ineigh,iatom))
        end if !end if
      end do !end do inu
    end do !end do imu
    intra_dip(1,iatom) = dip_x + res_dip(1,iatom)!/Debye
    intra_dip(2,iatom) = dip_y + res_dip(2,iatom)!/Debye
    intra_dip(3,iatom) = dip_z + res_dip(3,iatom)!/Debye
    dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  end do !end do iatom = 1,natoms
  dq_DP = 0.0d0
  do iatom = 1, natoms
    in1 = imass(iatom)
    r1(:) = ratom(:,iatom)
    rmedio = 0.0d0
    w_suma = 0.0d0
    do ineigh = 1,neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      w_k(ineigh)=exp(-klambda*((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
      rmedio = rmedio + w_k(ineigh)*r2
      w_suma = w_suma + w_k(ineigh)
    end do !end ineigh = 1,natoms
    rmedio = rmedio / w_suma
    bwrr = 0.0d0
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
  do iatom = 1, natoms
    Q0_TOT(iatom) = 0
    in1 = imass(iatom)
    do issh = 1, nssh(in1)
      Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
    end do
  end do
  dip_x = 0.0d0
  dip_y = 0.0d0
  dip_z = 0.0d0
  do iatom = 1, natoms
    Qtot=-Q0_TOT(iatom)+dq_DP(iatom)
    in1 = imass(iatom)
    do imu = 1,nssh(in1)
      Qtot = Qtot+Qout(imu,iatom)
    end do
    dip_x = dip_x+Qtot*ratom(1,iatom)
    dip_y = dip_y+Qtot*ratom(2,iatom)
    dip_z = dip_z+Qtot*ratom(3,iatom)
  end do !end do iatom = 1,natoms
  dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )    
  dipQout_x    = dip_x
  dipQout_y    = dip_y
  dipQout_z    = dip_z
  dipQout_tot = dip_tot
end subroutine dipole_proyection
  
function get_min_parabola(x1,x2,x3,y1,y2,y3)
  use M_constants, only: wp
  real(wp),intent(in) :: x1,x2,x3,y1,y2,y3
  real(wp) a,b,c,get_min_parabola
  b=(y2-y3)-(y1-y2)*(x2**2-x3**2)/(x1**2-x2**2)
  b=b/( (x2-x3)-(x1-x2)*(x2**2-x3**2)/(x1**2-x2**2) )
  a=(y1-y2)/(x1**2-x2**2)-b*(x1-x2)/(x1**2-x2**2)
  c=y1-a*x1**2-b*x1
  get_min_parabola=-b/(2*a)
  ! print*,a,'*x**2+',b,'*x+',c
end function get_min_parabola

