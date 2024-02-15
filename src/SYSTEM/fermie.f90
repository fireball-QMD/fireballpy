subroutine fermie ()
  use M_system
  use M_constants
  implicit none
  integer, parameter :: imax = 1000 ! maximum sc iterations
  integer, parameter :: nmax = 5000 ! cutoff for degeneracy check
  real, parameter :: tol = 1.0d-10
  integer ikpoint
  integer imu
  integer inu
  integer iter
  integer jkpoint
  real delta
  real emin
  real emax
  real qcharge
  real qztot
  real temp
  qztot = ztot 
  ioccupy_k = 0
  foccupy = 0.0d0
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
        eigen_k(inu,jkpoint) = (eigen_k(imu,ikpoint) + eigen_k(inu,jkpoint))/2.0d0
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
  emax=emax+0.20d0 
  iter = 0
  qcharge = 0.0d0
  do while (abs(qcharge - qztot) .gt. tol .and. iter .le. imax)
    iter = iter + 1
    efermi = (emax + emin)/2.0d0
    qcharge = 0.0d0
    do ikpoint = 1, nkpoints
      do imu = 1, norbitals
        delta = (eigen_k(imu,ikpoint) - efermi)/temp
        if (delta .gt. 10.0d0) then
          foccupy(imu,ikpoint) = 0.0d0
          ioccupy_k(imu,ikpoint) = 0
        else if (delta .lt. -10.0d0) then
          foccupy(imu,ikpoint) = 1.0d0
          ioccupy_k(imu,ikpoint) = 1
        else
          foccupy(imu,ikpoint) = 1.0d0/(1.0d0 + exp(delta))
          if (foccupy(imu,ikpoint) .gt. 1.0d-5) then
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

  if (iter .gt. imax) then
   write (*,*) '  '
   write (*,*) ' ************ WARNING ******** WARNING ************* '
   write (*,*) '  not under tolerance (toll) after ',imax
   write (*,*) '    iterations in subroutine fermie'
   write (*,*) '  '
   write (*,*) '    qcharge = ', qcharge
   write (*,*) '    qztot = ', qztot
   write (*,*) '    emax = ', emax
   write (*,*) '    emin = ', emin
   write (*,*) '    tol = ', tol, abs(qcharge - qztot),qcharge,qztot
   write (*,*) ' *************************************************** '
   write (*,*) '  '
  end if
  return
  end

