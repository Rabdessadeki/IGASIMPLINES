module mod_coktgbjp3gub


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine assemble_vector_ex01(ne1, ne2, ne3, ne4, p1, p2, p3, p4, &
        spans_1, spans_2, spans_3, spans_4, basis_1, basis_2, basis_3, &
        basis_4, weights_1, weights_2, weights_3, weights_4, points_1, &
        points_2, points_3, points_4, knots_1, knots_2, knots_3, &
        knots_4, vector_d, ovlp_value, S_DDM, domain_nb, rhs)

    implicit none

    integer(i64), value :: ne1
    integer(i64), value :: ne2
    integer(i64), value :: ne3
    integer(i64), value :: ne4
    integer(i64), value :: p1
    integer(i64), value :: p2
    integer(i64), value :: p3
    integer(i64), value :: p4
    integer(i64), intent(in) :: spans_1(0:)
    integer(i64), intent(in) :: spans_2(0:)
    integer(i64), intent(in) :: spans_3(0:)
    integer(i64), intent(in) :: spans_4(0:)
    real(f64), intent(in) :: basis_1(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_2(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_3(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_4(0:,0:,0:,0:)
    real(f64), intent(in) :: weights_1(0:,0:)
    real(f64), intent(in) :: weights_2(0:,0:)
    real(f64), intent(in) :: weights_3(0:,0:)
    real(f64), intent(in) :: weights_4(0:,0:)
    real(f64), intent(in) :: points_1(0:,0:)
    real(f64), intent(in) :: points_2(0:,0:)
    real(f64), intent(in) :: points_3(0:,0:)
    real(f64), intent(in) :: points_4(0:,0:)
    real(f64), intent(in) :: knots_1(0:)
    real(f64), intent(in) :: knots_2(0:)
    real(f64), intent(in) :: knots_3(0:)
    real(f64), intent(in) :: knots_4(0:)
    real(f64), intent(in) :: vector_d(0:,0:)
    real(f64), value :: ovlp_value
    real(f64), value :: S_DDM
    integer(i64), value :: domain_nb
    real(f64), intent(inout) :: rhs(0:,0:)
    integer(i64) :: k1
    integer(i64) :: k2
    integer(i64) :: ie1
    integer(i64) :: ie2
    integer(i64) :: i_spans1
    integer(i64) :: i_spans2
    integer(i64) :: il1
    integer(i64) :: i1
    integer(i64) :: il2
    integer(i64) :: i2
    real(f64) :: v
    integer(i64) :: ik1
    integer(i64) :: ik2
    real(f64) :: bx_0
    real(f64) :: bx_1
    real(f64) :: by_0
    real(f64) :: by_1
    real(f64) :: wval
    real(f64) :: x
    real(f64) :: y
    real(f64) :: f
    real(f64), allocatable :: lcoeffs_d(:,:)
    real(f64), allocatable :: lvalues_u(:)
    real(f64), allocatable :: lvalues_ux(:)
    integer(i64) :: nders
    integer(i64) :: degree
    real(f64), allocatable :: left(:)
    real(f64), allocatable :: right(:)
    real(f64), allocatable :: a(:,:)
    real(f64), allocatable :: ndu(:,:)
    real(f64), allocatable :: ders(:,:)
    real(f64), allocatable :: basis3(:,:)
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
    integer(i64) :: i1_ovrlp
    real(f64) :: neum_sign
    integer(i64) :: i_span_2
    integer(i64) :: g2
    integer(i64) :: span_3
    real(f64) :: u
    real(f64) :: ux
    integer(i64) :: il_1
    integer(i64) :: il_2
    real(f64) :: bi_0
    real(f64) :: bi_x
    real(f64) :: coeff_d
    real(f64) :: wsurf

    k1 = size(weights_1, 1_i64, i64)
    k2 = size(weights_2, 1_i64, i64)
    do ie1 = 0_i64, ne1 - 1_i64
      do ie2 = 0_i64, ne2 - 1_i64
        i_spans1 = spans_1(ie1)
        i_spans2 = spans_2(ie2)
        do il1 = 0_i64, p1 + 1_i64 - 1_i64
          i1 = i_spans1 - p1 + il1
          do il2 = 0_i64, p2 + 1_i64 - 1_i64
            i2 = i_spans2 - p2 + il2
            v = 0.0_f64
            do ik1 = 0_i64, k1 - 1_i64
              do ik2 = 0_i64, k2 - 1_i64
                bx_0 = basis_1(ik1, 0_i64, il1, ie1)
                bx_1 = basis_1(ik1, 1_i64, il1, ie1)
                by_0 = basis_2(ik2, 0_i64, il2, ie2)
                by_1 = basis_2(ik2, 1_i64, il2, ie2)
                wval = weights_1(ik1, ie1) * weights_2(ik2, ie2)
                x = points_1(ik1, ie1)
                y = points_2(ik2, ie2)
                f = 2_i64 * (3.141592653589793_f64 * &
                      3.141592653589793_f64) * sin( &
                      3.141592653589793_f64 * x) * sin( &
                      3.141592653589793_f64 * y)
                v = v + f * by_0 * bx_0 * wval
              end do
            end do
            rhs(i2 + p2, i1 + p1) = rhs(i2 + p2, i1 + p1) + v
          end do
        end do
      end do
    end do
    allocate(lcoeffs_d(0:p2 + 1_i64 - 1_i64, 0:p1 + 1_i64 - 1_i64))
    lcoeffs_d = 0.0_f64
    allocate(lvalues_u(0:k2 - 1_i64))
    lvalues_u = 0.0_f64
    allocate(lvalues_ux(0:k2 - 1_i64))
    lvalues_ux = 0.0_f64
    !---Computes All basis in a new points
    nders = 1_i64
    degree = p3
    !..
    allocate(left(0:degree - 1_i64))
    allocate(right(0:degree - 1_i64))
    allocate(a(0:degree + 1_i64 - 1_i64, 0:1_i64))
    allocate(ndu(0:degree + 1_i64 - 1_i64, 0:degree + 1_i64 - 1_i64))
    allocate(ders(0:degree + 1_i64 - 1_i64, 0:nders + 1_i64 - 1_i64))
    ders = 0.0_f64
    !...
    allocate(basis3(0:degree + 1_i64 - 1_i64, 0:nders + 1_i64 - 1_i64))
    basis3 = 0.0_f64
    do i = 0_i64, 0_i64
      !span = find_span( knots, degree, xq )
      xq = ovlp_value
      !~~~~~~~~~~~~~~~
      !span = find_span( knots, degree, xq )
      !~~~~~~~~~~~~~~~
      !Knot index at left/right boundary
      low = degree
      high = size(knots_3, kind=i64) - 1_i64 - degree
      !Check if point is exactly on left/right boundary, or outside domain
      if (xq <= knots_3(low)) then
        span = low
      end if
      if (xq >= knots_3(high)) then
        span = high - 1_i64
      else
        !Perform binary search
        span = FLOOR((low + high)/2.0_f64,i64)
        do while (xq < knots_3(span) .or. xq >= knots_3(span + 1_i64))
          if (xq < knots_3(span)) then
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
        left(j) = xq - knots_3(span - j)
        right(j) = knots_3(span + 1_i64 + j) - xq
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
      !Multiply derivatives by correct factors
      r = degree
      ders(:, 1_i64) = ders(:, 1_i64) * r
      basis3(:, 0_i64) = ders(:, 0_i64)
      basis3(:, 1_i64) = ders(:, 1_i64)
    end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Assembles Neumann Condition
    if (domain_nb == 0_i64) then
      i1_ovrlp = ne1 + 2_i64 * p1 - 1_i64
      neum_sign = 1.0_f64
    else
      i1_ovrlp = p1
      neum_sign = -1.0_f64
    end if
    do ie2 = 0_i64, ne2 - 1_i64
      i_span_2 = spans_2(ie2)
      do g2 = 0_i64, k2 - 1_i64
        !... We compute firstly the span in new adapted points
        !...
        xq = ovlp_value
        degree = p3
        low = degree
        high = size(knots_3, kind=i64) - 1_i64 - degree
        !Check if point is exactly on left/right boundary, or outside domain
        if (xq <= knots_3(low)) then
          span = low
        end if
        if (xq >= knots_3(high)) then
          span = high - 1_i64
        else
          !Perform binary search
          span = FLOOR((low + high)/2.0_f64,i64)
          do while (xq < knots_3(span) .or. xq >= knots_3(span + 1_i64))
            if (xq < knots_3(span)) then
              high = span
            else
              low = span
            end if
            span = FLOOR((low + high)/2.0_f64,i64)
          end do
        end if
        span_3 = span
        !...
        lcoeffs_d(:, :) = vector_d(i_span_2:i_span_2 + p4 + 1_i64 - &
              1_i64, span_3:span_3 + p3 + 1_i64 - 1_i64)
        !...
        u = 0.0_f64
        ux = 0.0_f64
        do il_1 = 0_i64, p3 + 1_i64 - 1_i64
          do il_2 = 0_i64, p4 + 1_i64 - 1_i64
            bi_0 = basis3(il_1, 0_i64) * basis_4(g2, 0_i64, il_2, ie2)
            bi_x = basis3(il_1, 1_i64) * basis_4(g2, 0_i64, il_2, ie2)
            !...
            coeff_d = lcoeffs_d(il_2, il_1)
            !...
            u = u + coeff_d * bi_0
            ux = ux + coeff_d * bi_x
          end do
        end do
        lvalues_u(g2) = u
        lvalues_ux(g2) = ux
      end do
      do il_2 = 0_i64, p2 + 1_i64 - 1_i64
        i2 = i_span_2 - p2 + il_2
        v = 0.0_f64
        do g2 = 0_i64, k2 - 1_i64
          bi_0 = basis_2(g2, 0_i64, il_2, ie2)
          wsurf = weights_2(g2, ie2)
          ux = lvalues_ux(g2)
          u = lvalues_u(g2)
          !..
          v = v + bi_0 * (neum_sign * ux + S_DDM * u) * wsurf
        end do
        rhs(i2 + p2, i1_ovrlp) = rhs(i2 + p2, i1_ovrlp) + v
      end do
    end do
    if (allocated(basis3)) then
      deallocate(basis3)
    end if
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
    if (allocated(lvalues_ux)) then
      deallocate(lvalues_ux)
    end if
    if (allocated(ndu)) then
      deallocate(ndu)
    end if

  end subroutine assemble_vector_ex01
  !........................................

end module mod_coktgbjp3gub
