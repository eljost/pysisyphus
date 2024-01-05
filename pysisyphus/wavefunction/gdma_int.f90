! Compile with
!     gfortran -Wall -Wextra -shared -Ofast gdma_int.f90 -o gdma_int_fortran.so
!

module mod_gdma_int

   use iso_c_binding, only: dp => c_double, int32 => c_int32_t

   implicit none

   real(dp), parameter :: NORM2 = 0.5773502691896258d0
   real(dp), parameter :: NORM3 = 0.2581988897471611d0
   real(dp), parameter :: NORM4 = 0.09759000729485332d0
   real(dp), parameter :: NORM22 = 0.3333333333333333d0
   ! Cutoff value for exp-argument
   real(dp), parameter :: cutoff = -36d0

contains

   subroutine prefacts(La, RA, res)
      integer(int32), intent(in) :: La
      real(dp), intent(in) :: RA(3)
      real(dp), intent(out) :: res(:)

      real(dp) :: dx, dy, dz, dx2, dy2, dz2
      real(dp) :: dx3, dy3, dz3, dx4, dy4, dz4

      dx = RA(1)
      dy = RA(2)
      dz = RA(3)

      if (La > 1) then
         dx2 = dx*dx
         dy2 = dy*dy
         dz2 = dz*dz
      end if
      if (La > 2) then
         dx3 = dx2*dx
         dy3 = dy2*dy
         dz3 = dz2*dz
      end if

      if (La > 3) then
         dx4 = dx2*dx2
         dy4 = dy2*dy2
         dz4 = dz2*dz2
      end if

      ! s-Orbital
      if (La == 0) then
         res(1) = 1.0
         ! p-Orbital
      elseif (La == 1) then
         res(1) = dx
         res(2) = dy
         res(3) = dz
         ! d-Orbital
      elseif (La == 2) then
         res(1) = NORM2*dx2
         res(2) = dx*dy
         res(3) = dx*dz
         res(4) = NORM2*dy2
         res(5) = dy*dz
         res(6) = NORM2*dz2
         ! f-Orbital
      elseif (La == 3) then
         res(1) = NORM3*dx3
         res(2) = NORM2*dx2*dy
         res(3) = NORM2*dx2*dz
         res(4) = NORM2*dx*dy2
         res(5) = dx*dy*dz
         res(6) = NORM2*dx*dz2
         res(7) = NORM3*dy3
         res(8) = NORM2*dy2*dz
         res(9) = NORM2*dy*dz2
         res(10) = NORM3*dz3
         ! g-Orbital
      elseif (La == 4) then
         res(1) = NORM4*dx4
         res(2) = NORM3*dx3*dy
         res(3) = NORM3*dx3*dz
         res(4) = NORM22*dx2*dy2
         res(5) = NORM2*dx2*dy*dz
         res(6) = NORM22*dx2*dz2
         res(7) = NORM3*dx*dy3
         res(8) = NORM2*dx*dy2*dz
         res(9) = NORM2*dx*dy*dz2
         res(10) = NORM3*dx*dz3
         res(11) = NORM4*dy4
         res(12) = NORM3*dy3*dz
         res(13) = NORM22*dy2*dz2
         res(14) = NORM4*dz4
      end if

   end subroutine prefacts

   pure integer function cart_size(L)
      integer(int32), intent(in) :: L
      cart_size = (L + 2)*(L + 1)/2
   end function cart_size

   subroutine eval_prim_density(nprims, npoints, nbfs, Ls_inds, primdata, coords3d, &
                                P, switch, rho) bind(c, name="eval_prim_density")
      ! Number of primitives, grid points and (contracted) basis functions
      integer(int32), intent(in), value :: nprims, npoints, nbfs
      ! Integer array; one column per primitive, containing angular momentum and
      ! starting index of shell in the density matrix.
      integer(int32), intent(in) :: Ls_inds(2, nprims)
      ! Double array; one column per primitive, containing contraction coefficient,
      ! orbital exponent and three center coordinates.
      real(dp), intent(in) :: primdata(5, nprims)
      ! Double array; one column per grid point, containing X, Y and Z coordinate.
      real(dp), intent(in) :: coords3d(3, npoints)
      ! Double array; Cartesian density matrix.
      real(dp), intent(in) :: P(nbfs, nbfs)
      ! Switch value; numerical integration is carried out for total exponents below
      ! this value.
      real(dp), intent(in), value :: switch
      ! Double array; holds the numerically integrated density.
      real(dp), intent(out) :: rho(npoints)

      ! Contraction coefficients da, db and orbital exponents ax, bx
      real(dp) :: da, ax, db, bx
      ! Total exponent px, reduced exponent mux
      real(dp) :: px, mux
      ! Pre-exponential factor K, exp-term Pexp of overlap distribution
      real(dp) :: K, Pexp
      ! Primitive centers A and B as well as grid point R
      real(dp) :: A(3), B(3), R(3)
      ! Center of overlap distribution Povlp, distance of P to grid point R RP
      ! and its squared array RP2
      real(dp) :: Povlp(3), RP(3), RP2(3)
      ! Angular moment dependent prefactors
      real(dp), allocatable :: prefacts_a(:), prefacts_b(:)
      ! Primitive indices ai and bi
      integer(int32) :: ai, bi
      ! Total angular momenta of primitives La and Lb
      integer(int32) :: La, Lb
      ! (Cartesian) starting indices and sizes of Cartesian shells
      integer(int32) :: cart_index_a, cart_index_b, cart_size_a, cart_size_b
      ! Indices
      integer(int32) :: nu, mu, i
      ! Maximum angular momentum of primitives
      integer(int32) :: Lmax
      ! Symmetry factor
      real(dp) :: factor

      ! Initialize density array
      rho = 0d0

      ! Allocate prefactor arrays once using the maximum L values
      Lmax = maxval(Ls_inds(1, :))
      allocate (prefacts_a(cart_size(Lmax)))
      allocate (prefacts_b(cart_size(Lmax)))

      ! Loop over primitive pairs/first primitive
      do ai = 1, nprims
         da = primdata(1, ai)
         ax = primdata(2, ai)
         A = primdata(3:5, ai)
         La = Ls_inds(1, ai)
         cart_index_a = Ls_inds(2, ai)
         cart_size_a = cart_size(La)

         ! Loop over second primitive
         do bi = ai, nprims
            db = primdata(1, bi)
            bx = primdata(2, bi)

            ! Only carry out numerical integration for diffuse total exponents,
            ! i.e, small total exponents px below the 'switch'-threshold.
            if (ax + bx >= switch) then
               cycle
            end if

            B = primdata(3:5, bi)
            Lb = Ls_inds(1, bi)
            cart_index_b = Ls_inds(2, bi)
            cart_size_b = cart_size(Lb)

            ! Calculate various quantities for the overlap distribution of two
            ! (primitive) Gaussians.
            !
            ! Total exponent
            px = ax + bx
            ! Reduced exponent
            mux = ax*bx/px
            ! Skip primitive pair when exp-argument for pre-exponential factor
            ! is too small.
            if (-mux*sum((A - B))**2 <= cutoff) then
               cycle
            end if

            ! Take symmetry into account, as we only loop over unique primitive
            ! combinations.
            factor = 2.0_dp
            if (ai == bi) factor = 1.0_dp

            ! Pre-exponential factor w/ contraction coeffcients and symmetry factor
            K = da*db*factor*exp(-mux*sum((A - B)**2))
            ! Center-of-charge coordinate
            Povlp = (ax*A + bx*B)/px

            ! Loop over grid points
            do i = 1, npoints
               R = coords3d(:, i)
               RP = R - Povlp
               RP2 = RP**2
               ! Check if exp-argument is below the threshold, if so we skip
               ! this grid point.
               if (-px*sum(RP2) <= cutoff) then
                  cycle
               end if
               Pexp = exp(-px*sum(RP2))

               ! Determine angular momentum dependent prefactors of Gaussian
               ! overlap distribution for both (primitive) shells.
               call prefacts(La, R - A, prefacts_a)
               call prefacts(Lb, R - B, prefacts_b)

               ! Contract everything with the appropriate density matrix
               ! entries.
               do nu = 1, cart_size_a
                  do mu = 1, cart_size_b
                     rho(i) = rho(i) + &
                              K*Pexp &
                              *P(cart_index_a + nu, cart_index_b + mu) &
                              *prefacts_a(nu) &
                              *prefacts_b(mu)
                  end do
               end do

            end do  ! End of loop over grid points R
         end do  ! End of loop over prims b
      end do  ! End of loop over prims a

      deallocate (prefacts_a)
      deallocate (prefacts_b)
   end subroutine eval_prim_density

end module mod_gdma_int
