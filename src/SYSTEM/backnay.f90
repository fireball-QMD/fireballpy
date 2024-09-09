subroutine backnay ()
  use iso_c_binding
  use M_system, only: natoms, ratom, neigh_b, neigh_j, neighn, neigh_back, xl
  implicit none
  integer(c_long) :: iatom
  integer(c_long) :: icount
  integer(c_long) :: ineigh
  integer(c_long) :: jatom
  integer(c_long) :: jbeta
  integer(c_long) :: jneigh
  integer(c_long) :: katom
  integer(c_long) :: mbeta
  integer(c_long) :: iloop
  727   continue
  do iatom = 1, natoms
   do ineigh = 1, neighn(iatom)
    mbeta = neigh_b(ineigh,iatom)
    jatom = neigh_j(ineigh,iatom)
    icount = 0
    do jneigh = 1, neighn(jatom)
     jbeta = neigh_b(jneigh,jatom)
     katom = neigh_j(jneigh,jatom)
     if (katom .ne. iatom .or.          &  
     &   abs(xl(1,jbeta) + xl(1,mbeta)) .gt. 0.001d0 .or.       &
     &   abs(xl(2,jbeta) + xl(2,mbeta)) .gt. 0.001d0 .or.       &
     &   abs(xl(3,jbeta) + xl(3,mbeta)) .gt. 0.001d0) then
     else
      icount = icount + 1
      neigh_back (iatom,ineigh) = jneigh
     end if
    end do
    if (icount .eq. 1) then
    else if (icount .eq. 0) then
     neighn(iatom) = neighn(iatom) - 1
     do iloop = ineigh, neighn(iatom)
      neigh_b(iloop,iatom) = neigh_b(iloop+1,iatom)
      neigh_j(iloop,iatom) = neigh_j(iloop+1,iatom)
     end do
     goto 727
    else if (icount .gt. 1) then ! BAD-Lots of debug information
     write (*,*) ' icount =' , icount
     write (*,*) ' The variable icount MUST be ONE. NO EXCEPTIONS! '
     write (*,*) ' Bad icount in backnay. Must abort.'
     write (*,*) ' iatom =' , iatom 
     write (*,*) ' ineigh =' , ineigh 
     write (*,*) ' mbeta =' , mbeta 
     write (*,*) ' jatom =' , jatom
     write (*,*) ' ratom(iatom) = ', ratom(:,iatom)
     write (*,*) ' ratom(jatom) = ', ratom(:,jatom)
     write (*,*) ' xl(mbeta) = ', xl(:,mbeta)
     do jneigh = 1, neighn(jatom) 
      jbeta = neigh_b(jneigh,jatom) 
      katom = neigh_j(jneigh,jatom)
      if (katom .ne. iatom .or.          &
     &    abs(xl(1,jbeta) + xl(1,mbeta)) .gt. 0.001d0 .or.       &
     &    abs(xl(2,jbeta) + xl(2,mbeta)) .gt. 0.001d0 .or.       &
     &    abs(xl(3,jbeta) + xl(3,mbeta)) .gt. 0.001d0) then
      else
       write (*,*) '  '
       write (*,*) ' A back neighbor jneigh = ', jneigh
       write (*,*) ' jbeta = ', jbeta
       write (*,*) ' xl(jbeta) = ', xl(:,jbeta) 
      end if
     end do
     stop
    end if
   end do
  end do
  return
end subroutine backnay
