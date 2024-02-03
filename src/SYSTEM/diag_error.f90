subroutine diag_error (info, istyle)
  implicit none
  integer, intent (in) :: info
  integer, intent (in) :: istyle
  write (*,*) '  '
  write (*,*) ' Diagonalization not successful, info = ', info
  if (info .lt. 0) then
    write (*,*) ' The ', info, '-th argument had an illegal '
    write (*,*) ' value. '
  else if(istyle .eq. 0)then
    write (*,*) ' It failed to converge.'
    write (*,*) info, ' off-diagonal elements of an intermediate'
    write (*,*) ' tridiagonal form did not converge to zero. '
  else if (mod(info,2) .ne. 0) then
    write (*,*) ' one or more eigenvectors failed to converge.'
    write (*,*) ' This should not have occured.  Send e-mail to'
    write (*,*) ' scalapack@cs.utk.edu if you feel like it'
  else if (mod(info/2,2) .ne. 0) then
    write (*,*) ' DARNGER Will Robinson.  Mr. Smith is in the house. '
    write (*,*) ' eigenvectors corresponding to one or more clusters '
    write (*,*) ' of eigenvalues could not be reorthogonalized'
    write (*,*) ' because of insufficient workspace.'
    write (*,*) ' We will blindly go on and hope for the best. '
    return
  else if (mod(info/4,2) .ne. 0) then
    write (*,*) ' space limit prevented computing all of the eigenvectors '
  else if (mod(info/8,2) .ne. 0) then
    write (*,*) ' PCSTEBZ failed to compute eigenvalues '
    write (*,*) ' This is very strange indeed.  It should not happen'
    write (*,*) ' Send e-mail to scalapack@cs.utk.edu if you feel like it'
  end if
  stop
end

