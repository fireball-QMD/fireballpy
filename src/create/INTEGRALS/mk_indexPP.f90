! mk_indexPP.f
! This program makes the indeces for the pseduopotential.
! Now there is a special "trick" with this.
! In create, we call mk_index. This sets everything up.
! The we call THE SAME PROGRAM twocenter to do the two center integrals,
! of which vnl is one of them. The problem is that the shell structure
! for the non-local pseudopotential is not the sam as the
!**********************************************************
!* JOM   MK_INDEX_3C (fortran77)  29/06/98
!* Fireball2000-version
!**********************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! We index (inm=1,inmmax) the non-zero 3-center matrix  elements H(i,j)
! between atom1 and atom2;
! For each non-zero matrix  element (inm) we keep the following
! information :
!    N-shell of left wavefunction : NLEFT (inm)
!    N-shell of right wavefunction : NRIGHT (inm)
!    L-value of left wavefunction : LLEFT (inm)
!    L-value of right wavefunction : LRIGHT (inm)
!    M-value of left wavefunction : MLEFT (inm)
!    M-value of right wavefunction : MRIGHT (inm)
!    number of non-zero elements (2c): index_max2c
!    number of non-zero elements (3c): index_max3c
!
! Atoms 1 and 2 (bondcharge) are along the Z-axis;
! The third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
       subroutine mk_indexPP ( IN1, IN2, nspecmax, nshmax, inter3max, &
     &         NSSH, nsshPP, LSSH, LsshPP, NLEFT, LLEFT, MLEFT, &
     &                        NRIGHT, LRIGHT, MRIGHT, &
     &                        index_max2cPP)
       use precision
       IMPLICIT NONE
 
       INTEGER in1,in2,inm, nshmax, nspecmax, inter3max
       INTEGER NSSH(nspecmax)
       integer nsshPP(nspecmax)
       INTEGER LSSH(nspecmax,nshmax)
       integer LsshPP(nspecmax,nshmax)
       INTEGER NLEFT(nspecmax,nspecmax,inter3max)
       INTEGER NRIGHT(nspecmax,nspecmax,inter3max)
       INTEGER LLEFT(nspecmax,nspecmax,inter3max)
       INTEGER LRIGHT(nspecmax,nspecmax,inter3max)
       INTEGER MLEFT(nspecmax,nspecmax,inter3max)
       INTEGER MRIGHT(nspecmax,nspecmax,inter3max)
 
 
       INTEGER index_max2cPP(nspecmax,nspecmax)
 
 
! Local variables
 
       INTEGER I1, I2, L1, L2 , K
 
! Let's begin with interactions with M1 = M2    (M=0 case)
! These are all cases we get for the two center integrals
! The left is the orbital, and the right is the potential.
!
       INM = 0
       DO I1 = 1, NSSH(IN1)
         L1 = LSSH(IN1,I1)
         DO I2 = 1 , nsshPP(IN2)
           L2 = lsshPP(IN2,I2)
           DO K = -MIN (L1,L2) , MIN (L1,L2)
             INM = INM + 1
             NLEFT (in1,in2,INM) = I1
             NRIGHT(in1,in2,INM) = I2
             LLEFT (in1,in2,INM) = L1
             LRIGHT(in1,in2,INM) = L2
             MLEFT (in1,in2,INM) = K
             MRIGHT(in1,in2,INM) = K
           END DO
         END DO
       END DO
!
       index_max2cPP(in1,in2) = inm
 
       RETURN
       END