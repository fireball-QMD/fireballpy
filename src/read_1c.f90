! read_1c.f90
! Program Description
! ===========================================================================
!       This routine reads in the 1-center (exchange-correlation)
! interactions. These 1-center interactions are contributions as described
! in the Horsefield approach.  This routine also reads in the variables which
! are needed to compute changes of the exchange correlation for the of charge
! transfer
!
! ===========================================================================
        subroutine read_1c ()
        !use charges
        !use dimensions
        !use integrals
        !use interactions
        use M_fdata
        implicit none
        integer :: nsup = 0

        integer iline
        integer Nlines_vdip1c_max
        integer trash
        integer in1, in2
        integer ins
        integer issh
        integer isorp
        integer itype
        integer jssh
        integer kssh
        integer kkssh
        integer numsh
 
        character (len=70) message
 
        logical skip_it

        real, dimension (nspecies) :: idshell

        integer, dimension (nsh_max) :: imask
        integer ideriv
        integer iissh, jjssh
        character (len = 200) extension
        character (len = 200) filename
        character (len = 200) root
          


! ***************************************************************************
!
!            M c W E D A   E X C H A N G E - C O R R E L A T I O N
!
! *************************************************************************** 
        if (itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then 
        
         allocate(exc1c0 (nspecies,nsh_max,nsh_max))
         allocate(nuxc1c (nspecies,nsh_max,nsh_max))
         allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
         allocate(d2exc1c (nspecies,nsh_max,nsh_max))
         allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
         allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))
        
         open (unit = 36, file = trim(fdataLocation)//'/xc1c_dqi.dat', status = 'unknown')
        
         do iline = 1, 4
          read (36,100) message
         end do
        
         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then
           read (36,*) itype, numsh

           do issh = 1, numsh
            read (36,*)
           end do
           read (36,*)
           do issh = 1, numsh
            read (36,*)
           end do

          else

           read (36,*) itype, numsh
           if (numsh .ne. nssh(in2)) then
            write (*,*) ' numsh .ne. nssh in read_1c.f90 '
            write (*,*) itype, numsh, in1, nssh(in2)
            stop
           end if

!           read (36,*) exc1c0(in2,iss)
           do issh = 1, numsh
            read (36,*) (exc1c0(in2,issh,jssh),jssh=1,numsh)
           end do
           read (36,*)
           do issh = 1, numsh
            read (36,*) (nuxc1c(in2,issh,jssh),jssh=1,numsh)
           end do
!           exc1c0(in2) = exc1c0(in2)
           do issh = 1, numsh
            do jssh = 1, numsh
             nuxc1c(in2,issh,jssh) = nuxc1c(in2,issh,jssh)
             exc1c0(in2,issh,jssh) = exc1c0(in2,issh,jssh)
            end do
           end do

! End loop over 2. derivative (only non-diagonal elements d^2 / dqi dqj)  

! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

! ==================================================================
!       READ FILE             nuxc1crho.dat
! ==================================================================
         open (unit = 36, file = trim(fdataLocation)//'/nuxc1crho.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (36,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
! skip unsed species
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then

           do 
            read (36,*) itype, numsh, kkssh
            do issh = 1, numsh
             read (36,*)
            end do
            if (numsh .eq. kkssh) exit
           end do ! do kssh 

          else

           do kssh = 1, nssh(in2)
            read (36,*) itype, numsh, kkssh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (36,*) (dnuxc1c(in2,issh,jssh,kssh),jssh=1,numsh)
            end do

            do issh = 1, numsh
             do jssh = 1, numsh
              dnuxc1c(in2,issh,jssh,kssh) = dnuxc1c(in2,issh,jssh,kssh)
             end do
            end do
           end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

! ==================================================================
!       READ FILE             exc1crho.dat
! ==================================================================
         open (unit = 36, file = trim(fdataLocation)//'/exc1crho.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (36,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
! skip unsed species
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then
           do
            read (36,*) itype, numsh, kkssh
            do issh = 1, numsh
             read (36,*)
            end do
            if (numsh .eq. kkssh) exit
           end do ! do kssh

          else

           do kssh = 1, nssh(in2)
            read (36,*) itype, numsh, kkssh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (36,*) (dexc1c(in2,issh,jssh,kssh),jssh=1,numsh)
            end do

            do issh = 1, numsh
             do jssh = 1, numsh
              dexc1c(in2,issh,jssh,kssh) = dexc1c(in2,issh,jssh,kssh)
             end do
            end do
           end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1
         if (itheory_xc .eq. 4) then
         end if !end if itheory_xc .eq. 4
         end if ! if(itheory_xc.eq.2 .or. itheory_xc .eq. 4) 


      !+++++++++++++++++++++++++++++++NEW JUNE 2019+++++++++++++++++++++++++++
      !.........................Vip 1c...........................................
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (V_intra_dip .eq. 1) then
      allocate(Nlines_vdip1c(nspecies))
      Nlines_vdip1c_max=0
      root = trim(fdataLocation)//'/vdip_onecenter'
      do in1 = 1,nspecies
       write (extension,'(''_'',i2.2)') in1
       filename = append_string (root,extension)
      open (unit = 36, file = filename, status = 'unknown')
      read(36,*) Nlines_vdip1c(in1)
      if (Nlines_vdip1c(in1) .gt. Nlines_vdip1c_max) then
      Nlines_vdip1c_max=Nlines_vdip1c(in1)
      end if
      close(36)
      end do !end do in1
        
       allocate(muR(Nlines_vdip1c_max,nspecies))
       allocate(nuR(Nlines_vdip1c_max,nspecies))
       allocate(alphaR(Nlines_vdip1c_max,nspecies))
       allocate(betaR(Nlines_vdip1c_max,nspecies))
       allocate(IR(Nlines_vdip1c_max,nspecies))

       muR    = 0.0d0
       nuR    = 0.0d0
       alphaR = 0.0d0
       betaR  = 0.0d0
       IR     = 0.0d0

         
      do in1 = 1,nspecies
       write (extension,'(''_'',i2.2)') in1
       filename = append_string (root,extension)
       open (unit = 36, file = filename, status = 'unknown')
       read(36,*) trash
         
       do iline = 1,Nlines_vdip1c(in1)
          
          
        read(36,*) muR(iline,in1), nuR(iline,in1), alphaR(iline,in1), betaR(iline,in1), IR(iline,in1)


      end do !end do iline = 1,Nlines_vdip1c

      close(36)
      write(*,*) 'Alles gut bisher' !Ankais
      end do !end do in1 = 1,nspecies

      end if ! if (V_intra_dip .eq. 1)
      !+++++++++++++++++++++++++++++++NEW JUNE 2019+++++++++++++++++++++++++++
      !.........................END OF Vip 1c...........................................
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (a70)
200     format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
!300     format (<numsh>f16.6)       
        return
        end subroutine read_1c
      
