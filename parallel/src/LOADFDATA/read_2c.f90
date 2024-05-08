subroutine read_2c (interaction)
  use M_fdata, only: fdataLocation, twocfname, nspecies, nssh, nsshPP, nzx, nfofx, cl_PP, errno2c, &
    & ind2c, index_max2c, index_maxS, index_maxPP, index_max2cDipX, index_max2cDipY, z2cmax, numz2c, nsh_max
  implicit none
  integer, intent (in) :: interaction
  integer in1, in2, initype, iounit, isorp, issh, isub2c, itype, npseudo, num_nonzero, numz, nzx1, nzx2
  real*8 rc1, rc2, zmax, zmin
  real*8, dimension (nsh_max) :: cl_pseudo
  character (len = 1000) extension, filename, root, root_isorp

  iounit = 71
  twocfname(1) = "overlap    "
  twocfname(2) = "vna_ontopl "
  twocfname(3) = "vna_ontopr "
  twocfname(4) = "vna_atom   "
  twocfname(5) = "vnl        "
  twocfname(6) = "xc_ontop   "
  twocfname(7) = "xc_atom    "
  twocfname(8) = "xc_corr    "
  twocfname(9) = "dipole_z   "
  twocfname(10) = "dipole_y   "
  twocfname(11) = "dipole_x   "
  twocfname(12) = "coulomb    "
  twocfname(13) = "kinetic    "
  twocfname(14) = "nuxc       "
  twocfname(15) = "den_ontopl "
  twocfname(16) = "den_ontopr "
  twocfname(17) = "den_atom   "
  twocfname(18) = "dnuxc_ol   "
  twocfname(19) = "dnuxc_or   "
  twocfname(20) = "denS_ontopl"
  twocfname(21) = "denS_ontopr"
  twocfname(22) = "denS_atom  "
  twocfname(23) = "overlapS   "
  root = trim(fdataLocation)//trim(twocfname(interaction))
  isub2c = 0
  if ((interaction .ge. 6) .and. (interaction .le. 8)) isub2c = 4

  do in1 = 1, nspecies
    do in2 = 1, nspecies
      initype = 0
      if(interaction .ge. 15 ) initype = 1
      if(interaction .eq. 23 ) initype = 0
      if (interaction .eq. 2) isub2c = nssh(in1)
      if (interaction .eq. 3) isub2c = nssh(in2)
      if (interaction .eq. 4) isub2c = nssh(in2)
      if (interaction .eq. 15) isub2c = nssh(in1)
      if (interaction .eq. 16) isub2c = nssh(in2)
      if (interaction .eq. 17) isub2c = nssh(in2)
      if (interaction .eq. 18) isub2c = nssh(in1)
      if (interaction .eq. 19) isub2c = nssh(in2)
      if (interaction .eq. 20) isub2c = nssh(in1)
      if (interaction .eq. 21) isub2c = nssh(in2)
      if (interaction .eq. 22) isub2c = nssh(in2)
      do isorp = initype, isub2c
        root_isorp = root
        if (isub2c .ge. 1) then  
          write (extension,"(""_"",i2.2)") isorp
          root_isorp = trim(root)//trim(extension) !append_string (root,extension)
        end if
        nzx1 = nzx(in1)
        nzx2 = nzx(in2)
        write (extension,"(""."",i2.2,""."",i2.2)") nzx1, nzx2
        filename = trim(root_isorp)//trim(extension) !append_string (root_isorp, extension)
        write (extension,"("".dat"")")
        filename = trim(filename)//trim(extension) !append_string (filename, extension)
        open(iounit, file=filename, status="old")
        call readheader_2c (interaction, iounit, numz, rc1, rc2, zmin, zmax, npseudo, cl_pseudo)
        if (numz .gt. nfofx) then
          !write (*,*) " numz = ", numz, " in read_2c.f90"
          !write (*,*) " nfofx = ",nfofx
          !write (*,*) " Fix this parameter and recompile! "
          !stop
          errno2c = 1
          return
        end if
        if (interaction .eq. 5) then
          do issh = 1, nsshPP(in2) 
            cl_PP(issh,in2) = cl_pseudo(issh)
          end do
        end if
        itype = ind2c(interaction,isorp)
        z2cmax(itype,in1,in2) = zmax
        numz2c(itype,in1,in2) = numz
        num_nonzero = index_max2c(in1,in2)
        if (interaction .eq. 5) num_nonzero = index_maxPP(in1,in2)
        if (interaction .eq. 4 .or. interaction .eq. 7) num_nonzero = index_max2c(in1,in1)
        if (interaction .eq. 10)  num_nonzero = index_max2cDipY(in1,in2)
        if (interaction .eq. 11)  num_nonzero = index_max2cDipX(in1,in2)
        if (interaction .eq. 17)  num_nonzero = index_max2c(in1,in1)
        if (interaction .eq. 22)  num_nonzero = index_maxS(in1,in1)
        if (interaction .eq. 12 .or. interaction .eq. 14) num_nonzero = nssh(in1)*nssh(in2)
        if (interaction .eq. 20 .or. interaction .eq. 21) num_nonzero = index_maxS(in1,in2)
        if (interaction .eq. 23) num_nonzero = index_maxS(in1,in2)

        call readdata_2c (interaction, iounit, num_nonzero, numz, zmax,itype, in1, in2)
        close(iounit)
      end do
    end do
  end do
end subroutine read_2c
