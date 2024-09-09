subroutine getforces ()
  use iso_c_binding
  use M_system
  implicit none
  call Dassemble_2c ()
  call Dassemble_2c_PP ()
  call Dassemble_ca_olsxc_on ()
  call Dassemble_ca_olsxc_2c ()
  if (idipole .eq. 0) call Dassemble_ca_2c ()
  if (idipole .eq. 1) call Dassemble_ca_2c_dip ()
  call Dassemble_3c ()
  call Dassemble_3c_PP ()
  if (idipole .eq. 0) call Dassemble_ca_3c ()
  if (idipole .eq. 1) call Dassemble_ca_3c_dip ()
  if (idipole .eq. 0) call Dassemble_lr ()
  if (idipole .eq. 1) call Dassemble_lr_dip ()
  if (iqmmm .eq. 1) then
    if (idipole .eq. 0) call Dassemble_qmmm ()
    if (idipole .eq. 1) call Dassemble_qmmm_dip ()
  else
    flrew_qmmm = 0.0d0
  end if
  call Dassemble_ca_olsxc_3c ()
  call assemble_F ()
 end subroutine getforces
