subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax, itype, in1, in2)
    use M_fdata
    implicit none
    integer, intent (in) :: in1, in2
    integer, intent (in) :: interaction
    integer, intent (in) :: iounit
    integer, intent (in) :: itype
    integer, intent (in) :: num_nonzero
    integer, intent (in) :: numz
    real, intent (in) :: zmax
    integer ipoint
    integer integral
    real, dimension (ME2c_max, nfofx) :: gstore
    if (interaction .ne. 8) then
        do ipoint = 1, numz
            read (iounit,*) (gstore(integral,ipoint), integral = 1, num_nonzero)
        end do
        do ipoint = 1, numz
            do integral = 1, num_nonzero
                xintegral_2c(integral,ipoint,itype,in1,in2) =  gstore(integral,ipoint)
            end do
        end do
    else
        do ipoint = 1, numz
            read (iounit,*) gstore(1,ipoint)
            xintegral_2c(1,ipoint,itype,in1,in2) = gstore(1,ipoint)
        end do
    end if

    return
end subroutine readdata_2c
