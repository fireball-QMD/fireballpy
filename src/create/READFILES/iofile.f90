! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
! fireball-qmd is a free (GPLv3) open project.

!
! This program is free software: you can redistribute it and/or modify
! by the Government is subject to restrictions as set forth in
! subdivsion { (b) (3) (ii) } of the Rights in Technical Data and
! Computer Software clause at 52.227-7013.
 
      subroutine iofile(root,suffix,index,iunit)
      use constants
      use precision
      implicit none
 
      integer index
      integer iunit
 
      integer i, j, lr, ls
      integer nplace
      integer nfile
 
      parameter (nplace=3,nfile=256)
      character(len=*) root,suffix
      character(len=80) form
      character(len=(nfile)) filen
 
!
! If root='junk/junker'
!    suffix='data'
!    index=12
!    iunit=25
! Then open file junk/junker012.data as unit 25 if open succeeds
! the parameter nplace gives the number of decimal places for the
! index number
! Warning: this file does internal i/o.  This will not work with OpenMP,
!          because all threads share the same internal i/o device.
!
      if (nplace.ge.10) then
         write (6,'('' nplace greater than 9 in iofile'')')
         stop 'error in iofile'
      endif
!
! find true length of root and suffix
!
      lr=len(root)
      do i=lr,1,-1
       if (root(i:i).ne.' ') go to 20
      end do
   20 lr=i
      ls=len(suffix)
      do i=ls,1,-1
       if (suffix(i:i).ne.' ') go to 40
      end do
   40 ls=i
      if (lr+ls+nplace+1.gt.len(filen)) then
         write (6,*) ' filen too short for input in iofile'
         write (6,*) root,suffix,index,iunit
         stop 'error in iofile'
      endif
!
! form file name
! fill in root here
!
      do i=1,lr
       filen(i:i)=root(i:i)
      end do
      if (index.lt.0.or.index.ge.10**nplace) then
         write (6,'('' index out of range in iofile'')')
         write (6,*) root,suffix,index,iunit
         stop 'error in iofile'
      endif
!
! write the format to form this will be (i3.3) for nplace=3
!
      write (form,'(''(i'',i1,''.'',i1,'')'')') nplace,nplace
!
! fill in number here
!
      write (filen(lr+1:lr+nplace),form) index
!
! fill in the dot here
!
      filen(lr+nplace+1:lr+nplace+1)='.'
!
! fill in the suffix here
!
      do i=1,ls
       j=lr+nplace+1+i
       filen(j:j)=suffix(i:i)
      end do
!
! fill rest with spaces in case there is garbage there now
!
      do i=lr+nplace+ls+2,len(filen)
       filen(i:i)=' '
      end do
!
! filename is formed, open the file
!
      open(unit=iunit,file=filen,status='unknown',err=80)
      if (verbose) write (*,*) filen
      return
   80 write (6,'('' error opening file in iofile '')')
      write (6,*) root,suffix,index,iunit
      stop 'error in iofile'
      end
