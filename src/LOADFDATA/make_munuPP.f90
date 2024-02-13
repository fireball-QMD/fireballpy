subroutine make_munuPP ()
    use M_fdata
    implicit none
    integer imu
    integer index
    integer in1, in2
    integer issh1, issh2
    integer l1, l2
    integer n1, n2

    allocate (index_maxPP (nspecies, nspecies))
    allocate (num_orbPP (nspecies))

    do in1 = 1, nspecies
        num_orbPP(in1) = 0
        do issh1 = 1 , nsshPP(in1)
            l1 = lsshPP(issh1,in1)
            num_orbPP(in1) = num_orbPP(in1) + 2*l1 + 1
        end do
    end do

    ME2cPP_max = 0
    do in1 = 1, nspecies
        do in2 = 1, nspecies
            index = 0
            do issh1 = 1, nssh(in1)
                l1 = lssh(issh1,in1)
                do issh2 = 1, nsshPP(in2) 
                    l2 = lsshPP(issh2,in2) 
                    do imu = -min(l1,l2), min(l1,l2) 
                        index = index + 1 
                    end do 
                end do 
            end do 
            if (index .gt. ME2cPP_max) ME2cPP_max = index 
        end do 
    end do

    allocate (muPP (ME2cPP_max, nspecies, nspecies)) 
    allocate (nuPP (ME2cPP_max, nspecies, nspecies)) 

    if (ME2cPP_max .gt. ME2c_max) ME2c_max = ME2cPP_max

    do in1 = 1, nspecies
        do in2 = 1, nspecies
            index = 0
            n1 = 0
            do issh1 = 1 , nssh(in1)
                l1 = lssh(issh1,in1)
                n1 = n1 + l1 + 1
                n2 = 0
                do issh2 = 1, nsshPP(in2)
                    l2 = lsshPP(issh2,in2)
                    n2 = n2 + l2 + 1
                    do imu = -min(l1,l2), min(l1,l2)
                        index = index + 1
                        muPP(index,in1,in2) = n1 + imu
                        nuPP(index,in1,in2) = n2 + imu
                    end do
                    n2 = n2 + l2
                end do
                n1 = n1 + l1
            end do
            index_maxPP(in1,in2) = index
        end do
    end do

    return
end subroutine make_munuPP
