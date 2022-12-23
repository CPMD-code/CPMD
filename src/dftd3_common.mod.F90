! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
!
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
! (for BJ-damping)
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
!
! Repackaging of the original code without any change in the functionality:
!
! Copyright (C) 2016, Bálint Aradi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

module dftd3_common
  implicit none

  ! Working precision (double precision)
  integer, parameter :: wp = kind(1.0d0)

  ! Large integers
  integer, parameter :: int64 = selected_int_kind(12)

  ! global ad hoc parameters
  real(wp), parameter :: k1 = 16.0
  real(wp), parameter :: k2 = 4./3.

  ! reasonable choices are between 3 and 5
  ! this gives smoth curves with maxima around the integer values
  ! k3=3 give for CN=0 a slightly smaller value than computed
  ! for the free atom. This also yields to larger CN for atoms
  ! in larger molecules but with the same chem. environment
  ! which is physically not right
  ! values >5 might lead to bumps in the potential
  real(wp), parameter :: k3 = -4.


  real(wp), parameter :: autoang = 0.52917726d0
  real(wp), parameter :: autokcal = 627.509541d0
  real(wp), parameter :: autoev = 27.21138505
  ! J/mol nm^6 - > au
  real(wp), parameter :: c6conv = 1.d-3/2625.4999d0/((autoang / 10.0d0)**6)


contains

  subroutine limit(iat,jat,iadr,jadr)
    integer, intent(inout) :: iat,jat
    integer, intent(out) :: iadr,jadr
    integer :: i
    iadr=1
    jadr=1
    i=100
    do while (iat .gt. 100)
      iat=iat-100
      iadr=iadr+1
    end do

    i=100
    do while (jat .gt.100)
      jat=jat-100
      jadr=jadr+1
    end do

  end subroutine limit

end module dftd3_common
