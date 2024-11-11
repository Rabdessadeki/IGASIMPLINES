module mod_7pbjpmwnmb1r


  use, intrinsic :: ISO_C_Binding, only : f64 => C_DOUBLE , i64 => &
        C_INT64_T
  implicit none

  contains

  !........................................
  subroutine assemble_vector_ex01(ne1, ne2, p1, p2, spans_1, spans_2, &
        basis_1, basis_2, weights_1, weights_2, points_1, points_2, rhs &
        )

    implicit none

    integer(i64), value :: ne1
    integer(i64), value :: ne2
    integer(i64), value :: p1
    integer(i64), value :: p2
    integer(i64), intent(in) :: spans_1(0:)
    integer(i64), intent(in) :: spans_2(0:)
    real(f64), intent(in) :: basis_1(0:,0:,0:,0:)
    real(f64), intent(in) :: basis_2(0:,0:,0:,0:)
    real(f64), intent(in) :: weights_1(0:,0:)
    real(f64), intent(in) :: weights_2(0:,0:)
    real(f64), intent(in) :: points_1(0:,0:)
    real(f64), intent(in) :: points_2(0:,0:)
    real(f64), intent(inout) :: rhs(0:,0:)
    integer(i64) :: k1
    integer(i64) :: k2
    real(f64), allocatable :: J_mat(:,:)
    integer(i64) :: ie1
    integer(i64) :: ie2
    integer(i64) :: i_spans1
    integer(i64) :: i_spans2
    integer(i64) :: g1
    integer(i64) :: g2
    real(f64) :: x
    real(f64) :: y
    real(f64) :: F1x
    real(f64) :: F1y
    real(f64) :: F2y
    real(f64) :: F2x
    real(f64) :: det_jac
    real(f64) :: x1
    real(f64) :: x2
    real(f64) :: f
    integer(i64) :: il1
    integer(i64) :: i1
    integer(i64) :: il2
    integer(i64) :: i2
    real(f64) :: v
    real(f64) :: bx_0
    real(f64) :: wvol

    k1 = size(weights_1, 1_i64, i64)
    k2 = size(weights_2, 1_i64, i64)
    allocate(J_mat(0:k2 - 1_i64, 0:k1 - 1_i64))
    J_mat = 0.0_f64
    do ie1 = 0_i64, ne1 - 1_i64
      do ie2 = 0_i64, ne2 - 1_i64
        i_spans1 = spans_1(ie1)
        i_spans2 = spans_2(ie2)
        do g1 = 0_i64, k1 - 1_i64
          do g2 = 0_i64, k2 - 1_i64
            x = points_1(g1, ie1)
            y = points_2(g2, ie2)
            F1x = 2.0_f64 * sqrt(1.0_f64 - 0.5_f64 * (2.0_f64 * y - &
                  1.0_f64) ** 2_i64)
            F1y = (-(2.0_f64 * x - 1.0_f64)) * (2.0_f64 * y - 1.0_f64) / &
                  sqrt(1.0_f64 - 0.5_f64 * (2.0_f64 * y - 1.0_f64) ** &
                  2_i64)
            F2y = 2.0_f64 * sqrt(1.0_f64 - 0.5_f64 * (2.0_f64 * x - &
                  1.0_f64) ** 2_i64)
            F2x = (-(2.0_f64 * x - 1.0_f64)) * (2.0_f64 * y - 1.0_f64) / &
                  sqrt(1.0_f64 - 0.5_f64 * (2.0_f64 * x - 1.0_f64) ** &
                  2_i64)
            det_jac = F1x * F2y - F2x * F1y
            x1 = (2.0_f64 * x - 1.0_f64) * sqrt(1.0_f64 - 0.5_f64 * ( &
                  2.0_f64 * y - 1.0_f64) ** 2_i64)
            x2 = (2.0_f64 * y - 1.0_f64) * sqrt(1.0_f64 - 0.5_f64 * ( &
                  2.0_f64 * x - 1.0_f64) ** 2_i64)
            f = 4.0_f64
            J_mat(g2, g1) = det_jac * f
          end do
        end do
        do il1 = 0_i64, p1 + 1_i64 - 1_i64
          i1 = i_spans1 - p1 + il1
          do il2 = 0_i64, p2 + 1_i64 - 1_i64
            i2 = i_spans2 - p2 + il2
            v = 0.0_f64
            do g1 = 0_i64, k1 - 1_i64
              do g2 = 0_i64, k2 - 1_i64
                bx_0 = basis_1(g1, 0_i64, il1, ie1) * basis_2(g2, 0_i64, &
                      il2, ie2)
                wvol = weights_1(g1, ie1) * weights_2(g2, ie2)
                v = v + bx_0 * wvol * J_mat(g2, g1)
              end do
            end do
            rhs(i2 + p2, i1 + p1) = rhs(i2 + p2, i1 + p1) + v
          end do
        end do
      end do
    end do
    if (allocated(J_mat)) then
      deallocate(J_mat)
    end if

  end subroutine assemble_vector_ex01
  !........................................

end module mod_7pbjpmwnmb1r
