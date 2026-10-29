subroutine generate_wavefunctions(nssh_in, nxc_in, nznuc_in, nexcite_in, nzval_in, nzval_ion_in, nzval_pp_in, &
                                  ioptim_in, ppfile_in, ppionfile_in, outpath_in, sav_in, lam_in, rcutoff_in, &
                                  xocc_in, xocc0_in, xocc_ion_in, cmix_in, exc_names_in, exc_weights_in, r0_in, v0_in)
  use, intrinsic :: iso_c_binding, only: c_double
  use begin_input, only: nssh, nexcite, nznuc, nzval, nzval_ion, nzval_pp, &
      ioptim, ppfile, ppionfile, outpath, sav, lam, rcutoff, &
      xocc, xocc0, xocc_ion, cmix, exc_names, exc_weights, r0, v0
  implicit none
  integer, intent(in) :: nssh_in, nxc_in, nexcite_in, nznuc_in, &
    & nzval_in, nzval_ion_in, ioptim_in
  real(c_double), intent(in) :: nzval_pp_in
  character(len=6), intent(in) :: ppfile_in
  character(len=8), intent(in) :: ppionfile_in
  character(len=1000), intent(in) :: outpath_in
  logical, intent(in) :: sav_in(nssh_in)
  integer, intent(in) :: lam_in(nssh_in)
  real(c_double), intent(in) :: rcutoff_in(nssh_in), xocc_in(nssh_in), xocc0_in(nssh_in), &
                                xocc_ion_in(nssh_in), cmix_in(nssh_in), r0_in(nssh_in), v0_in(nssh_in), &
                                exc_weights_in(nxc_in)
  character(len=32), intent(in) :: exc_names_in(nxc_in)

  ! Allocate module variables
  allocate(sav(nssh_in), lam(nssh_in), rcutoff(nssh_in), xocc(nssh_in), xocc0(nssh_in), &
      xocc_ion(nssh_in), cmix(nssh_in), exc_names(nxc_in), exc_weights(nxc_in))
  if (nexcite_in > 0) then
    allocate(r0(2*nssh_in), v0(2*nssh_in))
  else
    allocate(r0(nssh_in), v0(nssh_in))
  end if

  ! Assign module variables
  nexcite = nexcite_in
  nssh = nssh_in
  nznuc = nznuc_in
  nzval = nzval_in
  nzval_ion = nzval_ion_in
  nzval_pp = nzval_pp_in
  ioptim = ioptim_in
  ppfile = ppfile_in
  ppionfile = ppionfile_in
  sav = sav_in
  outpath = outpath_in
  lam = lam_in
  rcutoff = rcutoff_in
  xocc = xocc_in
  xocc0 = xocc0_in
  xocc_ion = xocc_ion_in
  cmix = cmix_in
  exc_names = exc_names_in
  exc_weights = exc_weights_in
  r0(:nssh_in) = r0_in
  v0(:nssh_in) = v0_in
  if (nexcite_in > 0) then
    r0(nssh_in+1:) = r0_in
    v0(nssh_in+1:) = v0_in
  end if

  ! Generate the files
  call rcatms ()

  ! Deallocate module variables
  deallocate(sav, lam, rcutoff, xocc, xocc0, xocc_ion, cmix, exc_names, exc_weights, r0, v0)
end subroutine generate_wavefunctions

subroutine generate_vnn(nssh_in, nexcite_in, nznuc_in, ppfile_in, ppionfile_in, outpath_in, lam_in, rcutoff_in)
  use, intrinsic :: iso_c_binding, only: c_double
  use begin_input, only: nssh, nexcite, nznuc, nzval_pp, ppfile, ppionfile, outpath, lam, rcutoff
  implicit none
  integer, intent(in) :: nssh_in, nexcite_in, nznuc_in
  character(len=6), intent(in) :: ppfile_in
  character(len=8), intent(in) :: ppionfile_in
  character(len=1000), intent(in) :: outpath_in
  integer, intent(in) :: lam_in(nssh_in)
  real(c_double), intent(in) :: rcutoff_in(nssh_in)

  ! Allocate module variables
  allocate(lam(nssh_in), rcutoff(nssh_in))

  ! Assign module variables
  nexcite = nexcite_in
  nssh = nssh_in
  nznuc = nznuc_in
  ppfile = ppfile_in
  ppionfile = ppionfile_in
  outpath = outpath_in
  lam = lam_in
  rcutoff = rcutoff_in
  nzval_pp = -1.0_c_double

  ! Generate the files
  call vnn ()

  ! Deallocate module variables
  deallocate(lam, rcutoff)
end subroutine generate_vnn
