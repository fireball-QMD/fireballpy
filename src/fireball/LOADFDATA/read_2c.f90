subroutine read_2c (interaction)
  use, intrinsic :: iso_fortran_env, only: double => real64
  use M_fdata, only: fdataLocation, twocfname, nspecies, nssh, nsshPP, nzx, nfofx, cl_PP, &
    & ind2c, index_max2c, index_maxS, index_maxPP, index_max2cDipX, index_max2cDipY, z2cmax, numz2c, nsh_max, initype
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
      if (interaction .eq. 2) isub2c = nssh(in1)
      if (interaction .eq. 3) isub2c = nssh(in2)
      if (interaction .eq. 4) isub2c = nssh(in2)
      if (interaction .eq. 7) isub2c = nssh(in1)
      if (interaction .eq. 8) isub2c = nssh(in2)
      if (interaction .eq. 14) isub2c = nssh(in1)
      if (interaction .eq. 15) isub2c = nssh(in2)
      if (interaction .eq. 16) isub2c = nssh(in2)
      if (interaction .eq. 17) isub2c = nssh(in1)
      if (interaction .eq. 18) isub2c = nssh(in2)
      if (interaction .eq. 19) isub2c = nssh(in2)
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
        if (interaction .eq. 5) then
          do issh = 1, nsshPP(in2) 
            cl_PP(issh,in2) = cl_pseudo(issh)
          end do
        end if
        itype = ind2c(interaction,isorp)
        z2cmax(itype,in1,in2) = zmax
        numz2c(itype,in1,in2) = numz
        num_nonzero = index_max2c(in1,in2)
        if (interaction .eq. 4) num_nonzero = index_max2c(in1,in1)
        if (interaction .eq. 5) num_nonzero = index_maxPP(in1,in2)
        if (interaction .eq. 10)  num_nonzero = index_max2cDipY(in1,in2)
        if (interaction .eq. 11)  num_nonzero = index_max2cDipX(in1,in2)
        if (interaction .eq. 16)  num_nonzero = index_max2c(in1,in1)
        if (interaction .eq. 19)  num_nonzero = index_maxS(in1,in1)
        if (interaction .eq. 12) num_nonzero = nssh(in1)*nssh(in2)
        if (interaction .eq. 17 .or. interaction .eq. 18) num_nonzero = index_maxS(in1,in2)
        if (interaction .eq. 20) num_nonzero = index_maxS(in1,in2)

        call readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
        close(iounit)
      end do
    end do
  end do
end subroutine read_2c
