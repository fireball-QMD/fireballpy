! These common blocks deal with DFT exchange-correlation
! You should not need to change any of this. But here it all is anyway.

        integer iderorb (nspec_max)           ! shell where charge is changed
        integer nsshxc (nspec_max)            ! number of shells
        integer lsshxc (nspec_max, nsh_max)   ! number of shells
 
        real(kind=long) xnocc (nsh_max, nspec_max)
        real(kind=long) dqorb (nspec_max)              ! dq used in changed shell
 
        common /xorbocc/ xnocc, dqorb, iderorb, nsshxc, lsshxc
! jel-dq
        real(kind=long) dqint (nsh_max, nspec_max)      ! intraatomic charge transfer
	common /dq/ dqint
