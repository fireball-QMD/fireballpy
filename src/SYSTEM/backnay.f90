! This is a neighbor routine, that we need for the nonlocal 
! pseudopotential.  For now we will worry about ontops.  Here is the problem.
! Suppose we need <i | VNL(i) |j> = <i|V(i)><V(i)|j>.  Now, suppose we are
! looping over iatom and jatom is the ineigh'th neighbor.  Well, the second 
! term of above is <V(iatom)|jatom>.  So we need 
! overlap(mu,nu,jatom,neighbor of j), but how do we know it.  We have 
! overlap(mu,nu,iatom,ineigh) and we are looping over iatom, ineigh, but
! what neighbor number is jatom to iatom. You see we know the neighbor
! number of jatom to iatom, but what is the neighbor number of iatom to jatom?
! Confusing isn't it.
subroutine backnay ()
  use M_system
  implicit none
  integer iatom
  integer icount
  integer ineigh
  integer jatom
  integer jbeta
  integer jneigh
  integer katom
  integer mbeta
  integer iloop
   
  727   continue
  do iatom = 1, natoms
    do ineigh = 1, neighn(iatom)
      mbeta = neigh_b(ineigh,iatom)
      jatom = neigh_j(ineigh,iatom)

      ! The question is: what neighbor number is iatom to ineigh?  We call this
      ! neighbor number neigh_back(ineigh,iatom).  We now loop over all neighbors
      ! of ineigh, to spot iatom.
      icount = 0
      do jneigh = 1, neighn(jatom)
        jbeta = neigh_b(jneigh,jatom)
        katom = neigh_j(jneigh,jatom)

        ! But just because katom is iatom does not prove we have the correct one.
        ! A specific atom may have several atoms of specific basis if the cell is 
        ! really small.  For instance, diamond has 4 neighbors all of the same basis
        ! in a 2 atom cell calculation.

        ! The correct result is: Given that katom = iatom (i.e. ratom_3 = ratom_1),
        ! then xl(ix,jbeta) = -xl(ix,mbeta).
        if (katom .ne. iatom .or.          &  
        &   abs(xl(1,jbeta) + xl(1,mbeta)) .gt. 0.001d0 .or.       &
        &   abs(xl(2,jbeta) + xl(2,mbeta)) .gt. 0.001d0 .or.       &
        &   abs(xl(3,jbeta) + xl(3,mbeta)) .gt. 0.001d0) then
        else
          ! Found it
          icount = icount + 1
          neigh_back (iatom,ineigh) = jneigh
        end if
      end do
      if (icount .eq. 1) then
        ! Everything is fine
      else if (icount .eq. 0) then
        ! Some neighbors are right on the cusp and will only get counted going
        ! one directions do to slight numerical differences.  We remove the
        ! extra neighbor and START OVER.
        neighn(iatom) = neighn(iatom) - 1
        do iloop = ineigh, neighn(iatom)
          neigh_b(iloop,iatom) = neigh_b(iloop+1,iatom)
          neigh_j(iloop,iatom) = neigh_j(iloop+1,iatom)
        end do
      goto 727
      else if (icount .gt. 1) then   ! BAD-Lots of debug information
        write (*,*) ' icount =' , icount
        write (*,*) ' The variable icount MUST be ONE. NO EXCEPTIONS  ! '
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
end
