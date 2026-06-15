subroutine read_2c (interaction)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: fdataLocation, twocfname, nspecies, nssh, nsshPP, nzx, nfofx, cl_PP, &
    & ind2c, index_max2c, index_maxS, index_maxPP, index_max2cDipX, index_max2cDipY, z2cmax, numz2c, nsh_max, initype, &
    & TWOCENTER_VNA_L, TWOCENTER_VNA_R, TWOCENTER_VNA_A, TWOCENTER_VXC_L, TWOCENTER_VXC_R, TWOCENTER_DEN_L, TWOCENTER_DEN_R, &
    & TWOCENTER_DEN_A, TWOCENTER_DENS_L, TWOCENTER_DENS_R, TWOCENTER_DENS_A, TWOCENTER_VNL, TWOCENTER_COULOMB, TWOCENTER_OVERLAPS, &
    & TWOCENTER_DIP_Y, TWOCENTER_DIP_X
  implicit none
  integer, intent (in) :: interaction
  integer :: in1, in2, iounit, isorp, issh, isub2c, itype, npseudo, num_nonzero, numz, nzx1, nzx2
  real(double) :: rc1, rc2, zmax, zmin
  real(double), dimension (nsh_max) :: cl_pseudo
  character (len = 1000) :: extension, filename, root, root_isorp

  iounit = 71
  root = trim(fdataLocation) // trim(twocfname(interaction))
  isub2c = 0
  !if (interaction .ge. 6) .and. (interaction .le. 8)) isub2c = 4

  do in1 = 1, nspecies
    do in2 = 1, nspecies

      if (interaction .eq. TWOCENTER_VNA_L) isub2c = nssh(in1)
      if (interaction .eq. TWOCENTER_VNA_R) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_VNA_A) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_VXC_L) isub2c = nssh(in1)
      if (interaction .eq. TWOCENTER_VXC_R) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_DEN_L) isub2c = nssh(in1)
      if (interaction .eq. TWOCENTER_DEN_R) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_DEN_A) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_DENS_L) isub2c = nssh(in1)
      if (interaction .eq. TWOCENTER_DENS_R) isub2c = nssh(in2)
      if (interaction .eq. TWOCENTER_DENS_A) isub2c = nssh(in2)
      !write(*,'("Read interaction ",I0," in read2c  isub2c=",I0, " twocf=",A)') interaction, isub2c, trim(twocfname(interaction))
      do isorp = initype(interaction), isub2c
        root_isorp = root
        if (isub2c .ge. 1) then  
          write (extension,'(''_'',i2.2)') isorp
          root_isorp = trim(root) // trim(extension) !append_string (root,extension)
        end if
        nzx1 = nzx(in1)
        nzx2 = nzx(in2)
        write (extension,'(''.'',i2.2,''.'',i2.2)') nzx1, nzx2
        filename = trim(root_isorp) // trim(extension) !append_string (root_isorp, extension)
        write (extension,'(''.dat'')')
        filename = trim(filename) // trim(extension) !append_string (filename, extension)
        open(iounit, file=filename, status='old')
        call readheader_2c (interaction, iounit, nsh_max, numz, rc1, rc2, zmin, zmax, npseudo, cl_pseudo)
        if (interaction .eq. TWOCENTER_VNL) then
          do issh = 1, nsshPP(in2)
            cl_PP(issh,in2) = cl_pseudo(issh)
          end do
        end if
        itype = ind2c(interaction,isorp)
        z2cmax(itype,in1,in2) = zmax
        numz2c(itype,in1,in2) = numz
        num_nonzero = index_max2c(in1,in2)
        if ((interaction .eq. TWOCENTER_VNA_A) .or. (interaction .eq. TWOCENTER_DEN_A))   num_nonzero = index_max2c(in1,in1)
        if (interaction .eq. TWOCENTER_VNL)   num_nonzero = index_maxPP(in1,in2)
        if (interaction .eq. TWOCENTER_DIP_Y)  num_nonzero = index_max2cDipY(in1,in2)
        if (interaction .eq. TWOCENTER_DIP_X)  num_nonzero = index_max2cDipX(in1,in2)
        if (interaction .eq. TWOCENTER_DENS_A)  num_nonzero = index_maxS(in1,in1)
        if (interaction .eq. TWOCENTER_COULOMB)  num_nonzero = nssh(in1)*nssh(in2)
        if ((interaction .eq. TWOCENTER_DENS_L) .or. (interaction .eq. TWOCENTER_DENS_R)) num_nonzero = index_maxS(in1,in2)
        if (interaction .eq. TWOCENTER_OVERLAPS) num_nonzero = index_maxS(in1,in2)

        call readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
        close(iounit)
      end do
    end do
  end do
end subroutine read_2c
