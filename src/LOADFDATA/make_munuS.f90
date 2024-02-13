subroutine make_munuS ()
    use M_fdata
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    allocate (index_maxS (nspecies, nspecies))

    MES_max = 0

    do in1 = 1, nspecies
        do in2 = 1, nspecies
            index = 0
            do issh1 = 1 , nssh(in1)
                do issh2 = 1, nssh(in2)
                    index = index + 1
                end do
            end do
            if (index .gt. MES_max) MES_max = index
        end do
    end do

    allocate (muS (MES_max, nspecies, nspecies))
    allocate (mvalueS (MES_max, nspecies, nspecies))
    allocate (nuS (MES_max, nspecies, nspecies))

    do in1 = 1, nspecies
        do in2 = 1, nspecies
            index = 0
            n1 = 0
            do issh1 = 1 , nssh(in1)
                n1 = n1 + 1
                n2 = 0
                do issh2 = 1, nssh(in2)
                    n2 = n2 + 1
                    index = index + 1
                    muS(index,in1,in2) = n1
                    nuS(index,in1,in2) = n2
                    mvalueS(index,in1,in2) = 0
                end do
            end do
            index_maxS(in1,in2) = index
        end do
    end do

    return
end subroutine make_munuS
