module mod_smkm7i907a0p


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine assemble_vector_ex02(ne1, ne2, ne3, p1, p2, p3, spans_1, &
        spans_2, spans_3, basis_1, basis_2, basis_3, weights_1, &
        weights_2, weights_3, points_1, points_2, points_3, knots_1, &
        knots_2, knots_3, vector_d, ovlp_value, rhs)

    implicit none

    integer(i64), value :: ne1
    integer(i64), value :: ne2
    integer(i64), value :: ne3
    integer(i64), value :: p1
    integer(i64), value :: p2
    integer(i64), value :: p3
    integer(i64), intent(in) :: spans_1(0:)
    integer(i64), intent(in) :: spans_2(0:)
    integer(i64), intent(in) :: spans_3(0:)
    real(f64), intent(in) :: basis_1(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_2(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_3(0:,0:,0:,0:)
    real(f64), intent(in) :: weights_1(0:,0:)
    real(f64), intent(in) :: weights_2(0:,0:)
    real(f64), intent(in) :: weights_3(0:,0:)
    real(f64), intent(in) :: points_1(0:,0:)
    real(f64), intent(in) :: points_2(0:,0:)
    real(f64), intent(in) :: points_3(0:,0:)
    real(f64), intent(in) :: knots_1(0:)
    real(f64), intent(in) :: knots_2(0:)
    real(f64), intent(in) :: knots_3(0:)
    real(f64), intent(in) :: vector_d(0:,0:)
    real(f64), value :: ovlp_value
    real(f64), intent(inout) :: rhs(0:)
    integer(i64) :: k1
    integer(i64) :: k2
    real(f64), allocatable :: lcoeffs_d(:,:)
    real(f64), allocatable :: lvalues_u(:)
    integer(i64) :: nders
    integer(i64) :: degree
    real(f64), allocatable :: left(:)
    real(f64), allocatable :: right(:)
    real(f64), allocatable :: a(:,:)
    real(f64), allocatable :: ndu(:,:)
    real(f64), allocatable :: ders(:,:)
    real(f64), allocatable :: basis2(:)
    integer(i64) :: i
    real(f64) :: xq
    integer(i64) :: low
    integer(i64) :: high
    integer(i64) :: span
    integer(i64) :: j
    real(f64) :: saved
    integer(i64) :: r
    real(f64) :: temp
    integer(i64) :: s1
    integer(i64) :: s2
    integer(i64) :: k
    real(f64) :: d
    integer(i64) :: rk
    integer(i64) :: pk
    integer(i64) :: j1
    integer(i64) :: j2
    integer(i64) :: ij
    integer(i64) :: ie1
    integer(i64) :: i_span_1
    integer(i64) :: g1
    integer(i64) :: span_2
    real(f64) :: u
    integer(i64) :: il_1
    integer(i64) :: il_2
    real(f64) :: bi_0
    real(f64) :: coeff_d
    integer(i64) :: i1
    real(f64) :: v
    real(f64) :: wvol

    !... sizes
    k1 = size(weights_1, 1_i64, i64)
    k2 = size(weights_2, 1_i64, i64)
    !...
    !...
    allocate(lcoeffs_d(0:p3 + 1_i64 - 1_i64, 0:p2 + 1_i64 - 1_i64))
    lcoeffs_d = 0.0_f64
    !..
    allocate(lvalues_u(0:k1 - 1_i64))
    lvalues_u = 0.0_f64
    !---Computes All basis in a new points
    nders = 0_i64
    degree = p2
    !..
    allocate(left(0:degree - 1_i64))
    allocate(right(0:degree - 1_i64))
    allocate(a(0:degree + 1_i64 - 1_i64, 0:1_i64))
    allocate(ndu(0:degree + 1_i64 - 1_i64, 0:degree + 1_i64 - 1_i64))
    allocate(ders(0:degree + 1_i64 - 1_i64, 0:nders + 1_i64 - 1_i64))
    ders = 0.0_f64
    !...
    allocate(basis2(0:degree + 1_i64 - 1_i64))
    basis2 = 0.0_f64
    do i = 0_i64, 0_i64
      !span = find_span( knots, degree, xq )
      xq = ovlp_value
      !~~~~~~~~~~~~~~~
      !Knot index at left/right boundary
      low = degree
      high = size(knots_2, kind=i64) - 1_i64 - degree
      !Check if point is exactly on left/right boundary, or outside domain
      if (xq <= knots_2(low)) then
        span = low
      end if
      if (xq >= knots_2(high)) then
        span = high - 1_i64
      else
        !Perform binary search
        span = FLOOR((low + high)/2.0_f64,i64)
        do while (xq < knots_2(span) .or. xq >= knots_2(span + 1_i64))
          if (xq < knots_2(span)) then
            high = span
          else
            low = span
          end if
          span = FLOOR((low + high)/2.0_f64,i64)
          !compute inverse of knot differences and save them into lower triangular part of ndu
          !compute basis functions and save them into upper triangular part of ndu
        end do
      end if
      ndu(0_i64, 0_i64) = 1.0_f64
      do j = 0_i64, degree - 1_i64
        left(j) = xq - knots_2(span - j)
        right(j) = knots_2(span + 1_i64 + j) - xq
        saved = 0.0_f64
        do r = 0_i64, j + 1_i64 - 1_i64
          ndu(r, j + 1_i64) = 1.0_f64 / (right(r) + left(j - r))
          temp = ndu(j, r) * ndu(r, j + 1_i64)
          ndu(j + 1_i64, r) = saved + right(r) * temp
          saved = left(j - r) * temp
        end do
        ndu(j + 1_i64, j + 1_i64) = saved
      end do
      !Compute derivatives in 2D output array 'ders'
      ders(:, 0_i64) = ndu(degree, :)
      do r = 0_i64, degree + 1_i64 - 1_i64
        s1 = 0_i64
        s2 = 1_i64
        a(0_i64, 0_i64) = 1.0_f64
        do k = 1_i64, nders + 1_i64 - 1_i64
          d = 0.0_f64
          rk = r - k
          pk = degree - k
          if (r >= k) then
            a(0_i64, s2) = a(0_i64, s1) * ndu(rk, pk + 1_i64)
            d = a(0_i64, s2) * ndu(pk, rk)
          end if
          j1 = merge(1_i64, -rk, rk > -1_i64)
          j2 = merge(k - 1_i64, degree - r, r - 1_i64 <= pk)
          do ij = j1, j2 + 1_i64 - 1_i64
            a(ij, s2) = (a(ij, s1) - a(ij - 1_i64, s1)) * ndu(rk + ij, &
                  pk + 1_i64)
          end do
          do ij = j1, j2 + 1_i64 - 1_i64
            d = d + a(ij, s2) * ndu(pk, rk + ij)
          end do
          if (r <= pk) then
            a(k, s2) = (-a(k - 1_i64, s1)) * ndu(r, pk + 1_i64)
            d = d + a(k, s2) * ndu(pk, r)
          end if
          ders(r, k) = d
          j = s1
          s1 = s2
          s2 = j
        end do
      end do
      basis2(:) = ders(:, 0_i64)
    end do
    !... build rhs
    do ie1 = 0_i64, ne1 - 1_i64
      i_span_1 = spans_1(ie1)
      do g1 = 0_i64, k1 - 1_i64
        !... We compute firstly the span in new adapted points
        !...
        xq = ovlp_value
        degree = p2
        low = degree
        high = size(knots_2, kind=i64) - 1_i64 - degree
        !Check if point is exactly on left/right boundary, or outside domain
        if (xq <= knots_2(low)) then
          span = low
        end if
        if (xq >= knots_2(high)) then
          span = high - 1_i64
        else
          !Perform binary search
          span = FLOOR((low + high)/2.0_f64,i64)
          do while (xq < knots_2(span) .or. xq >= knots_2(span + 1_i64))
            if (xq < knots_2(span)) then
              high = span
            else
              low = span
            end if
            span = FLOOR((low + high)/2.0_f64,i64)
          end do
        end if
        span_2 = span
        !...
        lcoeffs_d(:, :) = vector_d(i_span_1:i_span_1 + p3 + 1_i64 - &
              1_i64, span_2:span_2 + p2 + 1_i64 - 1_i64)
        !...
        u = 0.0_f64
        do il_1 = 0_i64, p2 + 1_i64 - 1_i64
          do il_2 = 0_i64, p3 + 1_i64 - 1_i64
            bi_0 = basis2(il_1) * basis_3(g1, 0_i64, il_2, ie1)
            !...
            coeff_d = lcoeffs_d(il_2, il_1)
            !...
            u = u + coeff_d * bi_0
          end do
        end do
        lvalues_u(g1) = u
      end do
      do il_1 = 0_i64, p1 + 1_i64 - 1_i64
        i1 = i_span_1 - p1 + il_1
        v = 0.0_f64
        do g1 = 0_i64, k1 - 1_i64
          bi_0 = basis_1(g1, 0_i64, il_1, ie1)
          !...
          wvol = weights_1(g1, ie1)
          !...
          u = lvalues_u(g1)
          !...
          v = v + bi_0 * u * wvol
        end do
        rhs(i1 + p1) = rhs(i1 + p1) + v
      end do
    end do
    !...
    if (allocated(lvalues_u)) then
      deallocate(lvalues_u)
    end if
    if (allocated(left)) then
      deallocate(left)
    end if
    if (allocated(ders)) then
      deallocate(ders)
    end if
    if (allocated(right)) then
      deallocate(right)
    end if
    if (allocated(lcoeffs_d)) then
      deallocate(lcoeffs_d)
    end if
    if (allocated(a)) then
      deallocate(a)
    end if
    if (allocated(basis2)) then
      deallocate(basis2)
    end if
    if (allocated(ndu)) then
      deallocate(ndu)
    end if

  end subroutine assemble_vector_ex02
  !........................................

end module mod_smkm7i907a0p
