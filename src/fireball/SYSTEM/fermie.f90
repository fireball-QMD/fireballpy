subroutine fermie ()
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use M_constants, only: spin, kconvert
  use M_system, only: tempfe, efermi, eigen_k, norbitals, nkpoints, ioccupy_k, foccupy, ztot, weight_k
  implicit none
  integer, parameter :: imax = 1000 ! maximum sc iterations
  integer, parameter :: nmax = 5000 ! cutoff for degeneracy check
  real(dp), parameter :: tol = 1.0e-10_dp
  integer ikpoint
  integer imu
  integer inu
  integer iter
  integer jkpoint
  real(dp) delta
  real(dp) emin
  real(dp) emax
  real(dp) qcharge
  real(dp) qztot
  real(dp) temp
  qztot = ztot 
  ioccupy_k = 0
  foccupy = 0.0_dp
  temp = tempfe/kconvert

  if (norbitals**2*nkpoints**2 .lt. nmax) then
   emin = eigen_k(1,1)
   emax = eigen_k(norbitals,1)
   do ikpoint = 1, nkpoints
    do imu = 1, norbitals
     if (eigen_k(imu,ikpoint) .lt. emin) emin = eigen_k(imu,ikpoint)
     if (eigen_k(imu,ikpoint) .gt. emax) emax = eigen_k(imu,ikpoint)
     do jkpoint = ikpoint, nkpoints
      do inu = imu, norbitals
       if (abs(eigen_k(imu,ikpoint) - eigen_k(inu,jkpoint)) .lt. tol) then
        eigen_k(inu,jkpoint) = (eigen_k(imu,ikpoint) + eigen_k(inu,jkpoint))/2.0_dp
        eigen_k(imu,ikpoint) = eigen_k(inu,jkpoint)
       end if
      end do
     end do
    end do
   end do
  else
   emin = eigen_k(1,1)
   emax = eigen_k(norbitals,1)
   do ikpoint = 1, nkpoints
    do imu = 1, norbitals
     if (eigen_k(imu,ikpoint) .lt. emin) emin = eigen_k(imu,ikpoint)
     if (eigen_k(imu,ikpoint) .gt. emax) emax = eigen_k(imu,ikpoint)
    end do
   end do
  end if
  emax=emax+0.20_dp 
  iter = 0
  qcharge = 0.0_dp
  do while (abs(qcharge - qztot) .gt. tol .and. iter .le. imax)
    iter = iter + 1
    efermi = (emax + emin)/2.0_dp
    qcharge = 0.0_dp
    do ikpoint = 1, nkpoints
      do imu = 1, norbitals
        delta = (eigen_k(imu,ikpoint) - efermi)/temp
        if (delta .gt. 10.0_dp) then
          foccupy(imu,ikpoint) = 0.0_dp
          ioccupy_k(imu,ikpoint) = 0
        else if (delta .lt. -10.0_dp) then
          foccupy(imu,ikpoint) = 1.0_dp
          ioccupy_k(imu,ikpoint) = 1
        else
          foccupy(imu,ikpoint) = 1.0_dp/(1.0_dp + exp(delta))
          if (foccupy(imu,ikpoint) .gt. 1.0e-5_dp) then
            ioccupy_k(imu,ikpoint) = 1
          else
            ioccupy_k(imu,ikpoint) = 0
          end if
        end if
        qcharge = qcharge + spin*foccupy(imu,ikpoint)*weight_k(ikpoint)
      end do ! do imu
    end do ! do ikpoint
    if (qcharge .gt. qztot) then
      emax = efermi
    else
      emin = efermi
    end if
  end do

  return
  end subroutine fermie
