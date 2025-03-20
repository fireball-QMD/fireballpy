! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jianjun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Brigham Young University - Hao Wang
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! onecentervdip.f90  
! Program Description
! ===========================================================================
!
! ===========================================================================
! Code written iy:
! Dani i JOM 
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine onecentervdip (nsh_max, nspec, nspec_max, itype, nsshxc,    &
     &                   lsshxc, drr_rho, rcutoffa_max, what, signature)
        use precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nsh_max
        integer, intent (in) :: nspec
        integer, intent (in) :: nspec_max
        integer, intent (in) :: itype

        integer, intent (in), dimension (nspec_max, nsh_max) :: lsshxc
        integer, intent (in), dimension (nspec_max) :: nsshxc

        real(kind=long), intent (in), dimension (nspec_max) :: drr_rho
        real(kind=long), intent (in), dimension (nspec_max) :: rcutoffa_max


        character (len=70), intent (in), dimension (nspec_max) :: what

        character (len=70), intent (in) :: signature


! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
        real(kind=long), parameter :: eq2 = 14.39975d0
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irho
        integer irhop
        integer issh
        integer jssh
        integer lalpha
        integer lmu
        integer lqn
        integer malpha
        integer mmu
        integer mqn
        integer nnrho
        integer counter
        integer counter_ini
        integer num_orb
        integer Nlines
        integer Max_Nlines
        integer iline

        real(kind=long) coefficient
        real(kind=long) cg1
        real(kind=long) cg2
        real(kind=long) cg3
        real(kind=long) cg4
        real(kind=long) drho
        real(kind=long) psi1
        real(kind=long) psi2
        real(kind=long) psi3
        real(kind=long) psi4
        real(kind=long) r
        real(kind=long) rp
        real(kind=long) rhomax
        real(kind=long) rhomin
        real(kind=long) sumr
        real(kind=long) sumrp
        real(kind=long) aux
        real(kind=long) Rint
        real(kind=long) radial_integral

        real(kind=long), dimension (:), allocatable :: factor
        real(kind=long), dimension (:), allocatable :: rpoint

        real(kind=long), external :: psiofr

        integer :: ialpha, ibeta, imu, inu, issh1, issh2, issh3, issh4
        integer :: l,l1,l2,l3,l4,m,m1,m2,m3,m4
        integer :: il1,il2,il3,il4
        real(kind=long) :: gauntReal, I

        integer, dimension(:,:), allocatable :: orb2lm
        integer, dimension(:), allocatable :: orb2shell
        integer , dimension(:,:), allocatable :: Orbitals
        real(kind=long), dimension(:), allocatable :: Total

        character (len = 200) extension
        character (len = 200) filename
        character (len = 200) root
        character (len = 200) append_string



! Procedure
! ===========================================================================
! Open the file to store the onecenter data.

         !root = trim(fdataLocation)//'/vdip_onecenter' 
         root = 'coutput/vdip_onecenter' 
         write (extension,'(''_'',i2.2)') itype
         !filename = trim(root)//trim(extension)
         filename = append_string (root,extension)

         open (unit = 36, file = filename, status = 'unknown')


! Set up the header for the output file.
!        write (36,100)
!        write (36,*) ' All one center matrix elements created by:'
!        write (36,200) signature

      !  do itype = 1, nspec
      !   write (36,300) what(itype)
      !  end do
      !  write (36,100)


             !----------        ORB2LM
  
          !Count number of orbitals
                num_orb = 0
             do issh = 1,nsshxc(itype)
                l = lsshxc(itype,issh)
                num_orb = num_orb+2*l+1
             end do

        

                 !integer, dimension(:,:,:), allocatable :: orb2lm
     
                 !allocate(orb2lm(nspecies,num_max_orb??,2))
                 allocate(orb2lm(num_orb,2))

                  !do itype = 1,nspecies
                    counter = 1
                    do issh = 1,nsshxc(itype)
                       l = lsshxc(itype,issh)
                       do m = -l,l
                          orb2lm(counter,1) = l
                          orb2lm(counter,2) = m
                          counter = counter+1
                       end do !end do imu = -l,l
                    end do !end do issh
                  !end do !end do itype

             !-----------


             !-------- ORB2SHELL
              allocate (orb2shell(num_orb))
              counter = 1
              do issh = 1,nsshxc(itype)
                 counter_ini = counter
                 l = lsshxc(itype,issh)
                 do imu = counter_ini,counter_ini+2*l
                    orb2shell(imu) = issh
                    counter = imu+1
                 end do !end imu
              end do ! end do issh = 1,nssh(in1) 

 
        !do itype = 1, nspec
 
! Initialize the limits of integration for the radial integral.
! Set up the grid points for the rho integration.
!         rhomin = 0.0d0
!         rhomax = rcutoffa_max(itype)
!         drho = drr_rho(itype)
!         nnrho = int((rhomax - rhomin)/drho) + 1
!         allocate (rpoint(nnrho))
!         allocate (factor(nnrho))


!test 
!         do il1 = 1 , nsshxc(itype)
!          l1=lsshxc(itype,il1)
!         write(36,*)'itype =',itype,'l =',l1 
!         do irho = 1, nnrho
!           r = rpoint(irho)
!           if (r .lt. 1.0d-04) r = 1.0d-04
!           psi1 = psiofr(itype,il1,r)
!           write(36,*)irho,r,psi1
!         end do
!         end do



 
!         do irho = 1, nnrho
!          rpoint(irho) = float(irho - 1)*drho
! Set up the Simpson rule factors:
!          factor(irho) = 2.0d0*drho/3.0d0
!          if (mod(irho,2) .eq. 0) factor(irho) = 4.0d0*drho/3.0d0
!          if (irho .eq. 1 .or. irho .eq. nnrho) factor(irho) = drho/3.0d0
!         end do !irho

              !new
               Max_Nlines = 100  !Upper bound : n*n*(n-1)*(n+1)/4, with n=num_orb
               allocate(Orbitals(Max_Nlines,4))
               Orbitals = 0
               allocate(Total(Max_Nlines))
               Total = 0.0d0
               Nlines = 0
               write(*,*) 'Hi, nice to see you!' !Ankais
               do imu = 1, num_orb  !(itype)
               do inu = imu, num_orb  !(itype)
               do ialpha = 1, num_orb  !(itype)
               do ibeta = ialpha+1, num_orb  !(itype)
               !if (imu .ne. ialpha .and. imu .ne. ibeta .and. inu .ne. ialpha .and. inu .ne. ibeta) then
               !if condition then  !Impose Symmetry
                            l1 = orb2lm(imu,1)
                            l2 = orb2lm(inu,1)
                            l3 = orb2lm(ialpha,1)
                            l4 = orb2lm(ibeta,1)
                            m1 = orb2lm(imu,2)
                            m2 = orb2lm(inu,2)
                            m3 = orb2lm(ialpha,2)
                            m4 = orb2lm(ibeta,2)
                            issh1 = orb2shell(imu)
                            issh2 = orb2shell(inu)
                            issh3 = orb2shell(ialpha)
                            issh4 = orb2shell(ibeta)
                          if (abs(l3-l4) .eq. 1) then !Dipolar approx.
                          I = 0.0d0
                          do l = 0,4    !This l comes from the Laplace expansion of 1/|r-r'|
                              !Compute the radial integral...
                              !R
                              Rint = radial_integral(itype, l, issh1, issh2,issh3, issh4, nspec_max, drr_rho, rcutoffa_max)
                              I=I+Rint*gauntReal(l,l1,l2,l3,l4,m1,m2,m3,m4)
                              !
                          end do ! end do l = 0,4
                         if (I .gt. 0.000001) then
                         Nlines = Nlines+1
                         Orbitals(Nlines,1) = imu
                         Orbitals(Nlines,2) = inu
                         Orbitals(Nlines,3) = ialpha
                         Orbitals(Nlines,4) = ibeta
                         Total(Nlines) = I
                         end if !end (if I .ne. 0)
                         end if !end if (abs(l3-l4) .eq. 1) !Dipolar approx.
                         !end if ! Impose symmetry
                         end do !ibeta
                         end do !ialpha
                         end do !inu
                         end do !imu
              ! end of new
 
                         !nsshxc instead of nssh??
                         !lsshxc instead of lssh??

      
     ! Write the Integrals to file 
                   write(36,501) Nlines
 do iline = 1,Nlines
write(36, 500) Orbitals(iline,1), Orbitals(iline,2), Orbitals(iline,3), Orbitals(iline,4), Total(iline)
 end do ! end do iline = 1,Nlines

     ! End of writing Integrals to file



               deallocate(orb2lm)
               deallocate(orb2shell)

 
 
! ****************************************************************************
! End loop over the species.
      
!        end do  ! end do itype

        !write (36,*) '  '
        write (*,*) '  '
        write (*,*) ' Writing output to: coutput/vdip_onecenter.dat '
        write (*,*) '  '
        close (unit = 36)

! Format itatements
! ===========================================================================
100     format (70('='))
200     format (2x, a45)
300     format (a70)
400     format (8d20.10)
500     format (4i3,1f12.8)
501     format (i3)
 
        return
        end




















