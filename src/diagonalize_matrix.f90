subroutine diagonalize_matrix

    use variables    
    integer :: n=0  
   
    open(unit=10, file=filename_in, status='old')
    do i = 1, norbitals
        read(10, *) (yyyy(i,j), j=1,norbitals)
        write(*, *) (yyyy(i,j), j=1,norbitals)
    end do
    close(unit=10)

    lrwork = 3*norbitals - 2
    allocate (rwork(lrwork))

    lwork = 1
    allocate (work(lwork))
    call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work,-1, rwork, info)
    lwork = work(1)
    deallocate (work)
    allocate (work(lwork))
    

    call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, lwork, rwork, info)

    if (info == 0) then
      print *, "Eigenvalues:"
      do i = 1, norbitals
        write(*,*) eigen(i)
      end do
      print*,'-------'
!      print *, "Diagonal Matrix:"
!      do i = 1, norbitals
!        write(*, '(<2*norbitals>F5.2)') (yyyy(i,j), j=1,norbitals)
!      end do
    else
      print *, "Error en la diagonalización (info =", info, ")"
    end if

end subroutine

