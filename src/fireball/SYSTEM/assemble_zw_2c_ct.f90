subroutine assemble_zw_2c_ct ()
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_system, only: numorb_max, nsh_max, natoms, neigh_self, neighn, neigh_b, neigh_j, imass, ratom, xl, nssh, Qneutral, Qin, orb2shell, s_mat, dip, vxc_ca, g2nu,

  implicit none
  integer iatom
  integer icount
  integer icount_sav
  integer ierror
  integer imu
  integer in1
  integer in2
  integer in3
  integer ineigh
  integer interaction
  integer inu
  integer isorp
  integer issh
  integer jatom
  integer jcount
  integer jcount_sav
  integer jssh
  integer kforce
  integer matom
  integer matom2
  integer mbeta
  integer my_proc
  integer natomsp
  integer ix
  integer iy
  integer igamma
  integer :: count_l, count_l_ini, issh1, issh2
  integer l
  real dq1
  real dq2
  real dterm
  real dterm_1
  real dterm_2
  real dstn_temp
  real dxn
  real dqdc
  real rcutoff_j
  real rend
  real rend1
  real rend2
  real sterm_1
  real sterm_2
  real y
  real rcutoff_i
  real :: qmu, q0mu, dqmu
  real :: A,B
  real, dimension (numorb_max, numorb_max) :: bcca
  real, dimension (3, nsh_max, nsh_max) :: bccapx
  real, dimension (nsh_max, nsh_max) :: bccax
  real, dimension (3, 3, 3) :: deps
  real, dimension (3, 3) :: eps
  real, dimension (3) :: r1
  real, dimension (3) :: r2
  real, dimension (3) :: r21
  real, dimension (3) :: sighat
  real stn1
  real stn2
  dxcdcc_zw = 0.0d0
  bcca = 0.0d0
  bccax = 0.0d0
  bccapx = 0.0d0
  vxc_ca = 0.0d0
  if (Kscf .eq. 1) then   
    g2nu = 0.0d0
    g2nup = 0.0d0
    do iatom = 1, natoms
      r1(:) = ratom(:,iatom)
      in1 = imass(iatom)
      matom = neigh_self(iatom)
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
        call epsilon (r2, sighat, eps)
        call deps2cent (r1, r2, eps, deps)
        bcca = 0.0d0
        kforce = 1
        if (matom .ne. ineigh) then
          interaction =  14   
          in3 = in2
          isorp = 0
          call doscentrosS (interaction, isorp, kforce, in1, in2, in3,y, eps, bccax, bccapx)
          do issh1 = 1, nssh(in1)
            do issh2 = 1,nssh(in2)
              g2nu(issh1,issh2,ineigh,iatom)=bccax(issh1,issh2)  
              g2nup(:,issh1,issh2,ineigh,iatom)=bccapx(:,issh1,issh2) 
            end do 
          end do 
        else 
          do issh1 = 1, nssh(in1)
            do issh2 = 1, nssh(in2)
              g2nu(issh1,issh2,matom,iatom) = xcnu1c(issh1,issh2,in1)
            end do 
          end do 
        end if 
      end do 
    end do 
  end if 
  do iatom = 1, natoms
    matom = neigh_self(iatom)
    r1(:) = ratom(:,iatom)
    in1 = imass(iatom)
    dq1 = 0.0d0
    do issh = 1, nssh(in1)
       dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
    end do
    do ineigh = 1, neighn(iatom)   
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)
      matom2 = neigh_self(jatom)
      r2(:) = ratom(:,jatom) + xl(:,mbeta)
      in2 = imass(jatom)
      dq2 = 0.0d0
      do issh = 1, nssh(in2)
        dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
      end do
      bcca = 0.0d0
      do isorp = 1, nssh(in2)
        dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
        in3=in1
        do imu = 1,num_orb(in1)
          issh1=orb2shell(imu,in1)
          bcca(imu,imu) = bcca(imu,imu) + g2nu(issh1,isorp,ineigh,iatom)*dxn
          if (Kscf .eq. 1 .and. iqout .eq. 6) then 
            gvhxc(imu,imu,isorp,jatom,matom,iatom) = gvhxc(imu,imu,isorp,jatom,matom,iatom) + g2nu(issh1,isorp,ineigh,iatom)
          end if 
        end do 
        do issh1 = 1,nssh(in1)
          dqdc = (Qin(issh1,iatom) - Qneutral(issh1,in1))
          uxcdcc_zw = uxcdcc_zw - (Qin(issh1,iatom)-0.50*dqdc)*g2nu(issh1,isorp,ineigh,iatom)*dxn
          dxcdcc_zw(:,ineigh,iatom) = dxcdcc_zw(:,ineigh,iatom) + (Qin(issh1,iatom)-0.50*dqdc)*g2nup(:,issh1,isorp,ineigh,iatom)*dxn
        end do 
      end do 
      in3 = in1
      do inu = 1, num_orb(in3)
        do imu = 1, num_orb(in1)
          vxc_ca(imu,inu,matom,iatom) = vxc_ca(imu,inu,matom,iatom) +  bcca(imu,inu)
        end do 
      end do 
      if (iatom .ne. jatom .or. mbeta .ne. 0) then 
        r21(:) = r2(:) - r1(:)
        y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
        bcca = 0.0d0
        in3 = in2
        do isorp = 1, nssh(in1)
          dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              issh1=orb2shell(imu,in1)
              issh2=orb2shell(inu,in3)
              A=0.5*s_mat(imu,inu,ineigh,iatom)-dip(imu,inu,ineigh,iatom)/y
              B=0.5*s_mat(imu,inu,ineigh,iatom)+dip(imu,inu,ineigh,iatom)/y  
              bcca(imu,inu) = bcca(imu,inu)+ &
              & A*g2nu(issh1,isorp,matom,iatom)*dxn+B*g2nu(isorp,issh2,ineigh,iatom)*dxn
              if (Kscf .eq. 1 .and. iqout .eq. 6) then 
                gvhxc(imu,inu,isorp,iatom,ineigh,iatom) = &  
                & gvhxc(imu,inu,isorp,iatom,ineigh,iatom) + &
                & A*g2nu(issh1,isorp,matom,iatom)+B*g2nu(isorp,issh2,ineigh,iatom)
              end if 
            end do 
          end do 
        end do 
        in3 = in2
        do isorp = 1, nssh(in2)
          dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
          do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
              issh1=orb2shell(imu,in1)
              issh2=orb2shell(inu,in3)
              A=0.5*s_mat(imu,inu,ineigh,iatom)-dip(imu,inu,ineigh,iatom)/y
              B=0.5*s_mat(imu,inu,ineigh,iatom)+dip(imu,inu,ineigh,iatom)/y  
              bcca(imu,inu) = bcca(imu,inu)+ &
              & A*g2nu(issh1,isorp,ineigh,iatom)*dxn+B*g2nu(issh2,isorp,matom2,jatom)*dxn     
              if (Kscf .eq. 1 .and. iqout .eq. 6) then 
                gvhxc(imu,inu,isorp,jatom,ineigh,iatom) =  gvhxc(imu,inu,isorp,jatom,ineigh,iatom) + &
                & A*g2nu(issh1,isorp,ineigh,iatom)+B*g2nu(issh2,isorp,matom2,jatom)
              end if 
            end do 
          end do 
        end do 
        do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
            vxc_ca(imu,inu,ineigh,iatom) = vxc_ca(imu,inu,ineigh,iatom) + bcca(imu,inu)
          end do
        end do
      end if
    end do 
  end do 
  return
end subroutine assemble_zw_2c_ct
