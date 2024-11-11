module mod_6o9jbnk6atzw


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine assemble_vector_ex001(ne1, ne2, p1, p2, spans_1, spans_2, &
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

  end subroutine assemble_vector_ex001
  !........................................

end module mod_6o9jbnk6atzw
