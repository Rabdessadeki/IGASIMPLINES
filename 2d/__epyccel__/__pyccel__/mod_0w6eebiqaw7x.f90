module mod_0w6eebiqaw7x


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine assemble_mass_1D(ne, degree, spans, basis, weights, points, &
        matrix)

    implicit none

    integer(i64), value :: ne
    integer(i64), value :: degree
    integer(i64), intent(in) :: spans(0:)
    real(f64), intent(in) :: basis(0:,0:,0:,0:)
    real(f64), intent(in) :: weights(0:,0:)
    real(f64), intent(in) :: points(0:,0:)
    real(f64), intent(inout) :: matrix(0:,0:)
    integer(i64) :: k1
    integer(i64) :: ie1
    integer(i64) :: i_span_1
    integer(i64) :: il_1
    integer(i64) :: i1
    integer(i64) :: il_2
    integer(i64) :: i2
    real(f64) :: v
    integer(i64) :: g1
    real(f64) :: bi_0
    real(f64) :: bj_0
    real(f64) :: bi_x
    real(f64) :: bj_x
    real(f64) :: wvol

    !... sizes
    k1 = size(weights, 1_i64, i64)
    !... build matrices
    do ie1 = 0_i64, ne - 1_i64
      i_span_1 = spans(ie1)
      !evaluation dependant uniquement de l'element
      do il_1 = 0_i64, degree + 1_i64 - 1_i64
        i1 = i_span_1 - degree + il_1
        do il_2 = 0_i64, degree + 1_i64 - 1_i64
          i2 = i_span_1 - degree + il_2
          v = 0.0_f64
          do g1 = 0_i64, k1 - 1_i64
            bi_0 = basis(g1, 0_i64, il_1, ie1)
            bj_0 = basis(g1, 0_i64, il_2, ie1)
            bi_x = basis(g1, 1_i64, il_1, ie1)
            bj_x = basis(g1, 1_i64, il_2, ie1)
            wvol = weights(g1, ie1)
            v = v + bi_0 * bj_0 * wvol
          end do
          matrix(degree + i2 - i1, degree + i1) = matrix(degree + i2 - &
                i1, degree + i1) + v
        end do
      end do
    end do
    !... end of mass assembling

  end subroutine assemble_mass_1D
  !........................................

end module mod_0w6eebiqaw7x
