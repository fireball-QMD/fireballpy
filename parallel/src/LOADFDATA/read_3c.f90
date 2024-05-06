subroutine read_3c (interaction)
use M_fdata, only: nsh_max, numXmax, numYmax, ME3c_max, nspecies, icon3c, nssh, nzx, index_max3c, index_maxS, &
  & bcna_01, bcna_02, bcna_03, bcna_04, bcna_05, numx3c_bcna, numy3c_bcna, x3cmax_bcna, y3cmax_bcna, hx_bcna, hy_bcna, &
  & den3_01, den3_02, den3_03, den3_04, den3_05, numx3c_den3, numy3c_den3, x3cmax_den3, y3cmax_den3, hx_den3, hy_den3, &
  & den3S_01, den3S_02, den3S_03, den3S_04, den3S_05, fdataLocation, ntheta, threecfname, errno3c
implicit none
integer, intent (in) :: interaction
integer iounit, in1, in2, in3, index, isorp, itheta, maxtype, maxisorp, &
  &  minisorp, numx, numy, nz1, nz2, nz3, num_nonzero
real*8 xmax, ymax
real*8, dimension (:,:,:,:,:), allocatable :: xintegral
character (len=1000) extension, filename, root, root1, root2

iounit = 71
threecfname(1) = "bcna"
threecfname(2) = "xc3c"
threecfname(3) = "den3"
threecfname(4) = "deS3"
root1 = trim(fdataLocation) // trim(threecfname(interaction))

if (interaction .eq. 1) then
  maxtype = nsh_max
  bcna_01 = 0.0d0
  bcna_02 = 0.0d0
  bcna_03 = 0.0d0
  bcna_04 = 0.0d0
  bcna_05 = 0.0d0
  numx3c_bcna = 0
  numy3c_bcna = 0
  x3cmax_bcna = 0.0d0
  y3cmax_bcna = 0.0d0
!else if (interaction .eq. 2) then
!  maxtype = ideriv_max
!  xc3c_01 = 0.0d0
!  xc3c_02 = 0.0d0
!  xc3c_03 = 0.0d0
!  xc3c_04 = 0.0d0
!  xc3c_05 = 0.0d0
!  numx3c_xc3c = 0
!  numy3c_xc3c = 0
!  x3cmax_xc3c = 0.0d0
!  y3cmax_xc3c = 0.0d0
else if (interaction .eq. 3) then
  maxtype = nsh_max
  den3_01 = 0.0d0
  den3_02 = 0.0d0
  den3_03 = 0.0d0
  den3_04 = 0.0d0
  den3_05 = 0.0d0
  numx3c_den3 = 0
  numy3c_den3 = 0
  x3cmax_den3 = 0.0d0
  y3cmax_den3 = 0.0d0
else if (interaction .eq. 4 ) then 
  maxtype = nsh_max
  den3S_01 = 0.0d0
  den3S_02 = 0.0d0
  den3S_03 = 0.0d0
  den3S_04 = 0.0d0
  den3S_05 = 0.0d0
end if

do in1 = 1, nspecies
  do in2 = 1, nspecies
    num_nonzero = index_max3c(in1,in2)
    do in3 = 1, nspecies
      index = icon3c(in1,in2,in3)
      maxisorp = 0
      minisorp = 0
      if (interaction .eq. 1) maxisorp = nssh(in3)
      !if (interaction .eq. 2) maxtype = ideriv_max
      if (interaction .eq. 3 .or. interaction .eq. 4) then
        maxisorp = nssh(in3)
        minisorp = 1
      end if
      do isorp = minisorp, maxisorp
        do itheta = 1, ntheta
          write(extension,'(''_'',i2.2)') itheta
          root2 = trim(root1) // trim(extension) !append_string(root1,extension)
          root = root2
          write(extension,'(''_'',i2.2)') isorp
          root = trim(root2) // trim(extension) !append_string(root2,extension)
          nz1 = nzx(in1)
          nz2 = nzx(in2)
          nz3 = nzx(in3)
          write(extension,'(''.'',i2.2,''.'',i2.2,''.'',i2.2)') nz1, nz2, nz3
          filename = trim(root) // trim(extension) !append_string(root,extension)
          write(extension,'(''.dat'')')
          filename = trim(filename) // trim(extension) !append_string(filename,extension)
          open (unit = iounit, file = filename, status = 'old')
          call readheader_3c (iounit, numx, numy, xmax, ymax)
          if (errno3c .ne. 0) return
          if (itheta .eq. 1) then
            if (interaction .eq. 1) then
              x3cmax_bcna(isorp,index) = xmax
              y3cmax_bcna(isorp,index) = ymax
              numx3c_bcna(isorp,index) = numx
              numy3c_bcna(isorp,index) = numy
            !else if (interaction .eq. 2) then
            !  x3cmax_xc3c(isorp,index) = xmax
            !  y3cmax_xc3c(isorp,index) = ymax
            !  numx3c_xc3c(isorp,index) = numx
            !  numy3c_xc3c(isorp,index) = numy
            else if (interaction .eq. 3 .or. interaction .eq. 4) then
              x3cmax_den3(isorp,index) = xmax
              y3cmax_den3(isorp,index) = ymax
              numx3c_den3(isorp,index) = numx
              numy3c_den3(isorp,index) = numy
            end if
          end if
          xintegral = 0.0d0
          call readdata_3c (iounit,numx,numy,num_nonzero,isorp,maxtype,index,xintegral)
          if (interaction .eq. 1) then
            if (itheta .eq. 1) bcna_01 = xintegral
            if (itheta .eq. 2) bcna_02 = xintegral
            if (itheta .eq. 3) bcna_03 = xintegral
            if (itheta .eq. 4) bcna_04 = xintegral
            if (itheta .eq. 5) bcna_05 = xintegral
          !else if (interaction .eq. 2) then
          !  if (itheta .eq. 1) xc3c_01 = xintegral
          !  if (itheta .eq. 2) xc3c_02 = xintegral
          !  if (itheta .eq. 3) xc3c_03 = xintegral
          !  if (itheta .eq. 4) xc3c_04 = xintegral
          !  if (itheta .eq. 5) xc3c_05 = xintegral
          else if (interaction .eq. 3) then
            if (itheta .eq. 1) den3_01 = xintegral
            if (itheta .eq. 2) den3_02 = xintegral
            if (itheta .eq. 3) den3_03 = xintegral
            if (itheta .eq. 4) den3_04 = xintegral
            if (itheta .eq. 5) den3_05 = xintegral
          else if (interaction .eq. 4) then 
            if (itheta .eq. 1) den3S_01 = xintegral
            if (itheta .eq. 2) den3S_02 = xintegral
            if (itheta .eq. 3) den3S_03 = xintegral
            if (itheta .eq. 4) den3S_04 = xintegral
            if (itheta .eq. 5) den3S_05 = xintegral
          end if
          close (iounit)
        end do
      end do
    end do
  end do
end do
deallocate(xintegral)
end subroutine read_3c
