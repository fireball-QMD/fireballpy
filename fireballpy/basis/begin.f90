subroutine generate_wavefunctions(ioption_in, nexcite_in, nssh_in, nznuc_in, nzval_in, nzval_ion_in, &
    & ioptim_in, atomname_in, ppfile_in, ppionfile_in, outpath_in, sav_in, lam_in, a0_in, rcutoff_in, &
    & xocc_in, xocc0_in, xocc_ion_in, cmix_in, r0_in, v0_in, filename_wf_in, filename_ewf_in)
  use iso_c_binding
  use begin_input, only: ioption, nexcite, nssh, nznuc, nzval, nzval_ion, &
    & ioptim, atomname, ppfile, ppionfile, outpath, sav, lam, a0, rcutoff, &
    & xocc, xocc0, xocc_ion, cmix, r0, v0, filename_wf, filename_ewf
  implicit none
  integer, intent(in) :: ioption_in, nexcite_in, nssh_in, nznuc_in, &
    & nzval_in, nzval_ion_in, ioptim_in
  character(len=10), intent(in) :: atomname_in
  character(len=6), intent(in) :: ppfile_in
  character(len=8), intent(in) :: ppionfile_in
  character(len=1000), intent(in) :: outpath_in
  logical, dimension(nssh_in), intent(in) :: sav_in
  integer, dimension(nssh_in), intent(in) :: lam_in
  real(c_double), dimension(nssh_in), intent(in) :: a0_in, rcutoff_in, &
    & xocc_in, xocc0_in, xocc_ion_in, cmix_in, r0_in, v0_in
  character(len=11), dimension(nssh_in), intent(in) :: filename_wf_in
  character(len=12), dimension(nssh_in), intent(in) :: filename_ewf_in

  ! Allocate module variables
  allocate(sav(nssh_in))
  allocate(lam(nssh_in))
  allocate(a0(nssh_in), rcutoff(nssh_in), xocc(nssh_in), xocc0(nssh_in), xocc_ion(nssh_in), cmix(nssh_in))
  allocate(filename_wf(nssh_in), filename_ewf(nssh_in))
  if (nexcite_in .eq. 3) then
    allocate(r0(2*nssh_in), v0(2*nssh_in))
  else
    allocate(r0(nssh_in), v0(nssh_in))
  end if

  ! Assign module variables
  ioption = ioption_in
  nexcite = nexcite_in
  nssh = nssh_in
  nznuc = nznuc_in
  nzval = nzval_in
  nzval_ion = nzval_ion_in
  ioptim = ioptim_in
  atomname = atomname_in
  ppfile = ppfile_in
  ppionfile = ppionfile_in
  sav = sav_in
  outpath = outpath_in
  lam = lam_in
  a0 = a0_in
  rcutoff = rcutoff_in
  xocc0 = xocc0_in
  xocc = xocc_in
  cmix = cmix_in
  filename_wf = filename_wf_in
  filename_ewf = filename_ewf_in
  r0(1:nssh_in) = r0_in(1:nssh_in)
  v0(1:nssh_in) = v0_in(1:nssh_in)
  if (nexcite_in .eq. 3) then
    r0(nssh_in+1:2*nssh_in) = r0_in(1:nssh_in)
    v0(nssh_in+1:2*nssh_in) = v0_in(1:nssh_in)
  end if

  ! Generate the files
  call rcatms ()

  ! Deallocate module variables
  deallocate(sav)
  deallocate(lam)
  deallocate(a0, rcutoff, xocc, xocc0, xocc_ion, cmix)
  deallocate(filename_wf, filename_ewf)
  deallocate(r0, v0)
end subroutine generate_wavefunctions


subroutine generate_vnn(nexcite_in, nssh_in, nznuc_in, ppfile_in, ppionfile_in, outpath_in, lam_in, rcutoff_in, &
    & filename_wf_in, filename_ewf_in, filename_na_in, filename_ena_in)
  use iso_c_binding
  use begin_input, only: nexcite, nssh, nznuc, ppfile, ppionfile, outpath, lam, rcutoff, &
    & filename_wf, filename_ewf, filename_na, filename_ena
  implicit none
  integer, intent(in) :: nexcite_in, nssh_in, nznuc_in
  character(len=6), intent(in) :: ppfile_in
  character(len=8), intent(in) :: ppionfile_in
  character(len=1000), intent(in) :: outpath_in
  integer, dimension(nssh_in), intent(in) :: lam_in
  real(c_double), dimension(nssh_in), intent(in) :: rcutoff_in
  character(len=11), dimension(nssh_in), intent(in) :: filename_wf_in
  character(len=11), dimension(0:nssh_in), intent(in) :: filename_na_in
  character(len=12), dimension(nssh_in), intent(in) :: filename_ewf_in, filename_ena_in

  ! Allocate module variables
  allocate(lam(nssh_in))
  allocate(rcutoff(nssh_in))
  allocate(filename_wf(nssh_in), filename_ewf(nssh_in), filename_ena(nssh_in))
  allocate(filename_na(0:nssh_in))

  ! Assign module variables
  nexcite = nexcite_in
  nssh = nssh_in
  nznuc = nznuc_in
  ppfile = ppfile_in
  ppionfile = ppionfile_in
  outpath = outpath_in
  lam = lam_in
  rcutoff = rcutoff_in
  filename_wf = filename_wf_in
  filename_ewf = filename_ewf_in
  filename_na(0:nssh_in) = filename_na_in(0:nssh_in)
  filename_ena = filename_ena_in

  ! Generate the files
  call vnn ()

  ! Deallocate module variables
  deallocate(lam)
  deallocate(rcutoff)
  deallocate(filename_wf, filename_ewf, filename_ena)
  deallocate(filename_na)
end subroutine generate_vnn
