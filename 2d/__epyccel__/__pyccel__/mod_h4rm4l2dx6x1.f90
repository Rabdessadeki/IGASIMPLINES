module mod_h4rm4l2dx6x1


  use, intrinsic :: ISO_C_Binding, only : i64 => C_INT64_T , f64 => &
        C_DOUBLE
  implicit none

  contains

  !........................................
  subroutine assemble_norm_ex01(ne1, ne2, p1, p2, spans_1, spans_2, &
        basis_1, basis_2, weights_1, weights_2, points_1, points_2, &
        vector_u, rhs)

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
    real(f64), intent(in) :: vector_u(0:,0:)
    real(f64), intent(inout) :: rhs(0:,0:)
    integer(i64) :: k1
    integer(i64) :: k2
    real(f64), allocatable :: lcoeffs_u(:,:)
    real(f64), allocatable :: lvalues_u(:,:)
    real(f64), allocatable :: lvalues_ux(:,:)
    real(f64), allocatable :: lvalues_uy(:,:)
    real(f64) :: norm_l2
    real(f64) :: norm_H1
    integer(i64) :: ie1
    integer(i64) :: ie2
    integer(i64) :: i_span_1
    integer(i64) :: i_span_2
    integer(i64) :: il_1
    integer(i64) :: il_2
    real(f64) :: coef
    integer(i64) :: ik1
    integer(i64) :: ik2
    real(f64) :: b0
    real(f64) :: dbx
    real(f64) :: dby
    real(f64) :: v
    real(f64) :: w
    real(f64) :: wvol
    real(f64) :: x
    real(f64) :: y
    real(f64) :: u
    real(f64) :: ux
    real(f64) :: uy
    real(f64) :: uh
    real(f64) :: uhx
    real(f64) :: uhy

    !... sizes
    k1 = size(weights_1, 1_i64, i64)
    k2 = size(weights_2, 1_i64, i64)
    !...
    allocate(lcoeffs_u(0:p2 + 1_i64 - 1_i64, 0:p1 + 1_i64 - 1_i64))
    lcoeffs_u = 0.0_f64
    allocate(lvalues_u(0:k2 - 1_i64, 0:k1 - 1_i64))
    lvalues_u = 0.0_f64
    allocate(lvalues_ux(0:k2 - 1_i64, 0:k1 - 1_i64))
    lvalues_ux = 0.0_f64
    allocate(lvalues_uy(0:k2 - 1_i64, 0:k1 - 1_i64))
    lvalues_uy = 0.0_f64
    norm_l2 = 0.0_f64
    norm_H1 = 0.0_f64
    do ie1 = 0_i64, ne1 - 1_i64
      do ie2 = 0_i64, ne2 - 1_i64
        i_span_1 = spans_1(ie1)
        i_span_2 = spans_1(ie2)
        lvalues_ux(:, :) = 0.0_f64
        lvalues_u(:, :) = 0.0_f64
        lvalues_uy(:, :) = 0.0_f64
        lcoeffs_u(:, :) = vector_u(i_span_2:i_span_2 + p2 + 1_i64 - &
              1_i64, i_span_1:i_span_1 + p1 + 1_i64 - 1_i64)
        do il_1 = 0_i64, p1 + 1_i64 - 1_i64
          do il_2 = 0_i64, p2 + 1_i64 - 1_i64
            coef = lcoeffs_u(il_2, il_1)
            do ik1 = 0_i64, k1 - 1_i64
              do ik2 = 0_i64, k1 - 1_i64
                b0 = basis_1(ik1, 0_i64, il_1, ie1) * basis_2(ik2, 0_i64 &
                      , il_2, ie2)
                dbx = basis_1(ik1, 1_i64, il_1, ie1) * basis_2(ik2, &
                      0_i64, il_2, ie2)
                dby = basis_1(ik1, 0_i64, il_1, ie1) * basis_2(ik2, &
                      1_i64, il_2, ie2)
                lvalues_u(ik2, ik1) = lvalues_u(ik2, ik1) + coef * b0
                lvalues_ux(ik2, ik1) = lvalues_ux(ik2, ik1) + coef * dbx
                lvalues_uy(ik2, ik1) = lvalues_uy(ik2, ik1) + coef * dby
              end do
            end do
          end do
        end do
        v = 0.0_f64
        w = 0.0_f64
        do ik1 = 0_i64, k1 - 1_i64
          do ik2 = 0_i64, k2 - 1_i64
            wvol = weights_1(ik1, ie1) * weights_2(ik2, ie2)
            x = points_1(ik1, ie1)
            y = points_2(ik2, ie2)
            !... test 0
            !_______________________CommentBlock_______________________!
            !                                                          !
            !                    u    = x*y*(x-1)*(y-1)                !
            !                    ux   =  x*y*(y - 1) + y*(x - 1)*(y - 1)!
            !                    uy   = x*y*(x - 1) + x*(x - 1)*(y - 1)!
            !                                                          !
            !__________________________________________________________!
            u = sin(3.141592653589793_f64 * x) * sin( &
                  3.141592653589793_f64 * y)
            ux = 3.141592653589793_f64 * cos(3.141592653589793_f64 * x) &
                  * sin(3.141592653589793_f64 * y)
            uy = 3.141592653589793_f64 * sin(3.141592653589793_f64 * x) &
                  * cos(3.141592653589793_f64 * y)
            uh = lvalues_u(ik2, ik1)
            uhx = lvalues_ux(ik2, ik1)
            uhy = lvalues_uy(ik2, ik1)
            v = v + (u - uh) ** 2_i64 * wvol
            w = w + ((ux - uhx) ** 2_i64 + (uy - uhy) ** 2_i64) * wvol
          end do
        end do
        norm_l2 = norm_l2 + v
        norm_H1 = norm_H1 + w
      end do
    end do
    norm_l2 = sqrt(norm_l2)
    norm_H1 = sqrt(norm_H1)
    rhs(p2, p1) = norm_l2
    rhs(p2 + 1_i64, p1) = norm_H1
    !...
    if (allocated(lcoeffs_u)) then
      deallocate(lcoeffs_u)
    end if
    if (allocated(lvalues_uy)) then
      deallocate(lvalues_uy)
    end if
    if (allocated(lvalues_ux)) then
      deallocate(lvalues_ux)
    end if
    if (allocated(lvalues_u)) then
      deallocate(lvalues_u)
    end if

  end subroutine assemble_norm_ex01
  !........................................

end module mod_h4rm4l2dx6x1
