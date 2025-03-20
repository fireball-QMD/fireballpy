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
!    N-shell of left wavefunction : NLEFT  (inm)
!    N-shell of right wavefunction : NRIGHT (inm)
!    L-value of left wavefunction : LLEFT  (inm)
!    L-value of right wavefunction : LRIGHT (inm)
!    M-value of left wavefunction : MLEFT  (inm)
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
 
        SUBROUTINE MK_INDEX ( IN1, IN2, nspecmax, nshmax, inter3max,
     1                        NSSH, LSSH, NLEFT3, LLEFT3, MLEFT3,
     2                        NRIGHT3, LRIGHT3, MRIGHT3,
     3                        index_max2c, index_max3c)
       IMPLICIT NONE
 
       INTEGER in1,in2,inm, nshmax, nspecmax, inter3max
       INTEGER NSSH(nspecmax)
       INTEGER LSSH(nspecmax,nshmax)
       INTEGER NLEFT3(nspecmax,nspecmax,inter3max)
       INTEGER NRIGHT3(nspecmax,nspecmax,inter3max)
       INTEGER LLEFT3(nspecmax,nspecmax,inter3max)
       INTEGER LRIGHT3(nspecmax,nspecmax,inter3max)
       INTEGER MLEFT3(nspecmax,nspecmax,inter3max)
       INTEGER MRIGHT3(nspecmax,nspecmax,inter3max)
 
       integer intermax
       parameter (intermax = 200)
 
       INTEGER NLEFT(intermax)
       INTEGER NRIGHT(intermax)
       INTEGER LLEFT(intermax)
       INTEGER LRIGHT(intermax)
       INTEGER MLEFT(intermax)
       INTEGER MRIGHT(intermax)
 
       INTEGER index_max2c(nspecmax,nspecmax)
       INTEGER index_max3c(nspecmax,nspecmax)
 
 
! Local variables
 
       INTEGER I1, I2, L1, L2 , K
 
! Let's begin with interactions with M1 = M2    (M=0 case)
! These are all cases we get for the two center integrals
 
       INM = 0
       DO I1 = 1, NSSH(IN1)
         L1 = LSSH(IN1,I1)
         DO I2 = 1 , NSSH(IN2)
           L2 = LSSH(IN2,I2)
           DO K = -MIN (L1,L2) , MIN (L1,L2)
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = K
             MRIGHT (INM) = K
           END DO
         END DO
       END DO
!
       index_max2c(in1,in2) = inm
!
!      All further interactions are what comes
!      in addition in the three center case.
!      So: in=1....index_max2c(in1,in2) gets us all
!          2-center stuff
!
!          in=1....index_max3c(in1,in2)  gets us
!          all 3-center stuff
!
 
! Now, the interactions with M1 = M2 +- 1      (M=1 case)
!
       DO I1 = 1 , NSSH(IN1)
         L1 = LSSH(IN1,I1)
         DO I2 = 1 , NSSH(IN2)
           L2 = LSSH(IN2,I2)
!
           IF ( L1 .eq. 0 .AND. L2 .ne. 0 ) THEN
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = 0
             MRIGHT (INM) = +1
           END IF
!
           IF ( L1 .eq. 1 ) THEN
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = +1
             MRIGHT (INM) = 0
!
             IF ( L2 .ne. 0 ) THEN
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = 0
               MRIGHT (INM) = +1
             END IF
!
             IF ( L2 .eq. 2 ) THEN
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = +1
               MRIGHT (INM) = +2
!
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = -1
               MRIGHT (INM) = -2
             END IF
           END IF
!
           IF ( L1 .eq. 2 ) THEN
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = +1
             MRIGHT (INM) = 0
!
             IF ( L2 .ne. 0 ) THEN
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = 0
               MRIGHT (INM) = +1
!
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = -2
               MRIGHT (INM) = -1
!
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = +2
               MRIGHT (INM) = +1
             END IF
!
             IF ( L2 .eq. 2 ) THEN
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = +1
               MRIGHT (INM) = +2
!
               INM = INM + 1
               NLEFT  (INM) = I1
               NRIGHT (INM) = I2
               LLEFT  (INM) = L1
               LRIGHT (INM) = L2
               MLEFT  (INM) = -1
               MRIGHT (INM) = -2
             END IF
           END IF
!
         END DO
       END DO
 
! Now, the interactions with M1 = M2 +- 2      (M=2 case)
!
       DO I1 = 1 , NSSH(IN1)
         L1 = LSSH(IN1,I1)
         DO I2 = 1 , NSSH(IN2)
           L2 = LSSH(IN2,I2)
!
           IF ( L1 .eq. 2 ) THEN
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = +2
             MRIGHT (INM) = 0
           END IF
!
           IF ( L2 .eq. 2 ) THEN
             INM = INM + 1
             NLEFT  (INM) = I1
             NRIGHT (INM) = I2
             LLEFT  (INM) = L1
             LRIGHT (INM) = L2
             MLEFT  (INM) = 0
             MRIGHT (INM) = +2
           END IF
!
         END DO
       END DO
 
        if (inm.gt.inter3max) then
        write(*,*)' inm .gt. inter3max (inter_max) in mk_index.f'
        write(*,*)' Make inter_max greater than or equal to ',inm
        write(*,*)' in parameters.inc!'
        stop
        end if

       index_max3c(in1,in2)=inm
       do i1 = 1, inm
         nleft3 (in1,in2,i1)=nleft(i1)
         nright3(in1,in2,i1)=nright(i1)
         lleft3 (in1,in2,i1)=lleft(i1)
         lright3(in1,in2,i1)=lright(i1)
         mleft3 (in1,in2,i1)=mleft(i1)
         mright3(in1,in2,i1)=mright(i1)
       enddo
 
       RETURN
       END
