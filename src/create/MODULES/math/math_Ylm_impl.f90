submodule (math) math_lstsq_impl
  implicit none

contains

  module real(kind=sp) function sYlm(l, m, prec)
    integer, intent(in) :: l, m
    real(kind=sp), intent(in) :: prec

    if (l == 0) then
      sYlm = 1.0_sp
    else if (l == 1) then
      if (m == -1 .or. m == 1 ) then
        sYlm = 1.22474487139158904909_sp
      else if (m == 0) then
        sYlm = 1.73205080756887729352_sp
      else
        sYlm = 0.0_sp
      end if
    else if (l == 2) then
      if (m == -2 .or. m == 2) then
        sYlm = 1.36930639376291528364_sp
      else if (m == -1 .or. m == 1 ) then
        sYlm = 2.73861278752583056728_sp
      else if (m == 0) then
        sYlm = 1.11803398874989484820_sp
      else
        sYlm = 0.0_sp
      end if
    else
      sYlm = 0.0_sp
    end if
  end function sYlm

  module real(kind=dp) function dYlm(l, m, prec)
    integer, intent(in) :: l, m
    real(kind=dp), intent(in) :: prec

    if (l == 0) then
      dYlm = 1.0_dp
    else if (l == 1) then
      if (m == -1 .or. m == 1 ) then
        dYlm = 1.22474487139158904909_dp
      else if (m == 0) then
        dYlm = 1.73205080756887729352_dp
      else
        dYlm = 0.0_dp
      end if
    else if (l == 2) then
      if (m == -2 .or. m == 2) then
        dYlm = 1.36930639376291528364_dp
      else if (m == -1 .or. m == 1 ) then
        dYlm = 2.73861278752583056728_dp
      else if (m == 0) then
        dYlm = 1.11803398874989484820_dp
      else
        dYlm = 0.0_dp
      end if
    else
      dYlm = 0.0_dp
    end if
  end function dYlm

end submodule math_lstsq_impl
