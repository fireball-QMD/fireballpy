! This routine calculates the Ewald sum for a crystal with a given basis. This is specially designed for molecules with a given dipole moment.
subroutine get_ewald (iauxforce)
  use iso_c_binding
  use M_constants, only: pi
  use M_system, only: icluster, natoms, ratom, imass, ewald, dewald, fewald, Qin, a1vec, a2vec, a3vec
  use M_fdata, only: nssh, Qneutral
  implicit none
  integer(c_long), intent (in) :: iauxforce
  integer(c_long) iatom
  integer(c_long) ig1
  integer(c_long) ig2
  integer(c_long) ig3
  integer(c_long) ig1mx
  integer(c_long) ig2mx
  integer(c_long) ig3mx
  integer(c_long) il1
  integer(c_long) il2
  integer(c_long) il3
  integer(c_long) il1mx
  integer(c_long) il2mx
  integer(c_long) il3mx
  integer(c_long) in1
  integer(c_long) issh
  integer(c_long) ix
  integer(c_long) jatom
  real(c_double) argument
  real(c_double) my_erfc
  real(c_double) derfcdr
  real(c_double) distance
  real(c_double) factor
  real(c_double) factorf
  real(c_double) g1mag2, g2mag2, g3mag2
  real(c_double) gdotb
  real(c_double) gmax
  real(c_double) gmin2
  real(c_double) gsq
  real(c_double) kappa
  real(c_double) QQ
  real(c_double) r1mag2, r2mag2, r3mag2
  real(c_double) rmax
  real(c_double) rmin2
  real(c_double) stuff
  real(c_double) volcel
  real(c_double), dimension (3) :: cvec
  real(c_double), dimension (3) :: eta
  real(c_double), dimension (3, natoms) :: fewald1, fewald2
  real(c_double), dimension (3) :: g
  real(c_double), dimension (3) :: g1, g2, g3
  real(c_double), dimension (natoms) :: Q, Q0
  real(c_double), dimension (3) :: vecl
  ewald = 0.0d0
  if (iauxforce .eq. 1) dewald = 0.0d0
  if (iauxforce .eq. 1) fewald = 0.0d0
  do iatom = 1, natoms
   Q(iatom) = 0.0d0
   Q0(iatom) = 0.0d0
   in1 = imass(iatom)
   do issh = 1, nssh(in1)
    Q(iatom) = Q(iatom) + Qin(issh,iatom)
    Q0(iatom) = Q0(iatom) + Qneutral(issh,in1)
   end do
  end do
  call cross (a2vec, a3vec, cvec)
  volcel = a1vec(1)*cvec(1) + a1vec(2)*cvec(2) + a1vec(3)*cvec(3)
  g1(:) = 2.0d0*pi*cvec(:)/volcel
  call cross (a3vec, a1vec, cvec)
  g2(:) = 2.0d0*pi*cvec(:)/volcel
  call cross (a1vec, a2vec, cvec)
  g3(:) = 2.0d0*pi*cvec(:)/volcel
  volcel = abs(volcel)
  gmax = 5.0d0
  rmax = 5.0d0
  g1mag2 = g1(1)**2 + g1(2)**2 + g1(3)**2
  g2mag2 = g2(1)**2 + g2(2)**2 + g2(3)**2
  g3mag2 = g3(1)**2 + g3(2)**2 + g3(3)**2
  r1mag2 = a1vec(1)**2 + a1vec(2)**2 + a1vec(3)**2
  r2mag2 = a2vec(1)**2 + a2vec(2)**2 + a2vec(3)**2
  r3mag2 = a3vec(1)**2 + a3vec(2)**2 + a3vec(3)**2
  rmin2 = r1mag2
  if (r2mag2 .lt. rmin2) rmin2 = r2mag2
  if (r3mag2 .lt. rmin2) rmin2 = r3mag2
  gmin2 = g1mag2
  if (g2mag2 .lt. gmin2) gmin2 = g2mag2
  if (g3mag2 .lt. gmin2) gmin2 = g3mag2
  kappa = sqrt(sqrt(gmin2/(4.0d0*rmin2)))
  ig1mx = nint(gmax * sqrt(4.0d0*kappa**2/g1mag2) + 1.0d0, c_long)
  ig2mx = nint(gmax * sqrt(4.0d0*kappa**2/g2mag2) + 1.0d0, c_long)
  ig3mx = nint(gmax * sqrt(4.0d0*kappa**2/g3mag2) + 1.0d0, c_long)
  if (ig1mx .le. 1) ig1mx = 2
  if (ig2mx .le. 1) ig2mx = 2
  if (ig3mx .le. 1) ig3mx = 2
  il1mx = nint(rmax * sqrt(1.0d0/(kappa**2*r1mag2)) + 1.0d0, c_long)
  il2mx = nint(rmax * sqrt(1.0d0/(kappa**2*r2mag2)) + 1.0d0, c_long)
  il3mx = nint(rmax * sqrt(1.0d0/(kappa**2*r3mag2)) + 1.0d0, c_long)
  if (il1mx .le. 1) il1mx = 2
  if (il2mx .le. 1) il2mx = 2
  if (il3mx .le. 1) il3mx = 2
  if (iauxforce .eq. 1) fewald1 = 0.0d0
  if (icluster .eq. 1) then
   ig1mx = 0
   ig2mx = 0
   ig3mx = 0
  end if !AQUI pensar quitar
  do ig1 = -ig1mx, ig1mx
   do ig2 = -ig2mx, ig2mx
    do ig3 = -ig3mx, ig3mx
     if (.not. (ig1 .eq. 0 .and. ig2 .eq. 0 .and. ig3 .eq. 0)) then
      g(:) = ig1*g1(:) + ig2*g2(:) + ig3*g3(:)
      gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
      argument = gsq/(4.0d0*kappa**2)
      stuff = 4.0d0*pi*exp(-argument)/(gsq*volcel)
      do iatom = 1, natoms
       do jatom = iatom, natoms
        factor = 1.0d0*stuff
        factorf = 2.0d0*stuff
        if (jatom .eq. iatom) factor = 0.5d0*stuff
        gdotb = g(1)*(ratom(1,iatom) - ratom(1,jatom)) + g(2)*(ratom(2,iatom) - ratom(2,jatom)) + g(3)*(ratom(3,iatom) - ratom(3,jatom))
        QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
        ewald(iatom,jatom) = ewald(iatom,jatom) + factor*cos(gdotb)
        ewald(jatom,iatom) = ewald(jatom,iatom) + factor*cos(gdotb)
        if (iauxforce .eq. 1) then
         do ix = 1, 3
           fewald1(ix,iatom) =  fewald1(ix,iatom) + QQ*factorf*sin(gdotb)*g(ix)
         end do
         do ix = 1, 3
           fewald1(ix,jatom) =  fewald1(ix,jatom) - QQ*factorf*sin(gdotb)*g(ix)
         end do
         do ix = 1, 3
           dewald(ix,iatom,jatom) =  dewald(ix,iatom,jatom) - factor*sin(gdotb)*g(ix)
         end do
         do ix = 1, 3
           dewald(ix,jatom,iatom) =  dewald(ix,jatom,iatom) + factor*sin(gdotb)*g(ix)
         end do
        end if
       end do
      end do
     end if
    end do
   end do
  end do
  if (iauxforce .eq. 1) fewald2 = 0.0d0
  if (icluster .eq. 1) then
   il1mx = 0
   il2mx = 0
   il3mx = 0
   kappa = 0.0d0
  end if
  do il1 = -il1mx, il1mx
   do il2 = -il2mx, il2mx
    do il3 = -il3mx, il3mx
     do iatom = 1, natoms
      do jatom = iatom, natoms
       factor = 1.0d0
       factorf = 2.0d0
       if (jatom .eq. iatom) factor = 0.5d0
       QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
       vecl(:) = il1*a1vec(:) + il2*a2vec(:) + il3*a3vec(:)
       eta(:) = vecl(:) + ratom(:,iatom) - ratom(:,jatom)
       distance = sqrt(eta(1)**2 + eta(2)**2 + eta(3)**2)
       if (distance .gt. 0.0001d0) then
        argument = kappa*distance
        ewald(iatom,jatom) =  ewald(iatom,jatom) + factor*my_erfc(argument)/distance
        ewald(jatom,iatom) =  ewald(jatom,iatom) + factor*my_erfc(argument)/distance
        derfcdr = (2.0d0*exp(-argument**2)*kappa/sqrt(pi) + my_erfc(argument)/distance)/distance**2
        if (iauxforce .eq. 1) then

         do ix = 1, 3
          fewald2(ix,iatom) =  fewald2(ix,iatom) + QQ*eta(ix)*factorf*derfcdr
         end do
         do ix = 1, 3
          fewald2(ix,jatom) =  fewald2(ix,jatom) - QQ*eta(ix)*factorf*derfcdr
         end do
         do ix = 1, 3
          dewald(ix,iatom,jatom) =  dewald(ix,iatom,jatom) - eta(ix)*factor*derfcdr
         end do
         do ix = 1, 3
          dewald(ix,jatom,iatom) =  dewald(ix,jatom,iatom) + eta(ix)*factor*derfcdr
         end do
        end if
       end if
      end do
     end do
    end do
   end do
  end do
  do iatom = 1, natoms
   ewald(iatom,iatom) = ewald(iatom,iatom) - 2.0d0*kappa/sqrt(pi)
  end do
  if (iauxforce .eq. 1) fewald = fewald1 + fewald2
  return
end subroutine get_ewald

subroutine cross (a, b, c)
  use iso_c_binding
  implicit none
  real(c_double), intent(in), dimension(3) :: a
  real(c_double), intent(in), dimension(3) :: b
  real(c_double), intent(out), dimension(3) :: c
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
  return
end subroutine cross

function my_erfc (x)
  use iso_c_binding
  implicit none
  real(c_double), intent(in) :: x
  real(c_double) :: s, ax, t, my_erfc

  ax = abs(x)
  if( ax .lt. 1.0d-10) then
    my_erfc = 1.0d0
    return
  end if
  s = x/ax
  t = 1.0d0/(1.0d0 + 0.3275911d0*ax)
  my_erfc = 1.0d0 - s*(1.0d0 - exp(-ax*ax)*t*(0.254829592d0 + t*(-0.284496736d0 + t*(1.421413741d0 + t*(-1.453152027d0 + t*1.061405429d0)))))
end function my_erfc
