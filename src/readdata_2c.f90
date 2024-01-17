! readdata_2c.f90
! Program Description
! ===========================================================================
!       This routine reads the data from the 2-center integral files. When
! read, the information is stored in the array xintegral_2c.  This array
! is the field that stores all non-vanishing matrix elements for a general
! 2-center integral.  There are maximal ME2c_max non-vanishing matrix
! elements given on a grid of maximal nfofx data points.  The exact dimensions
! for a given interaction, and a given pair of atoms are numz and num_nonzero.
!
! ===========================================================================
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
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        

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
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end subroutine readdata_2c
