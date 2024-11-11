module bind_c_mod_pvrvdwly0zgb

  use mod_pvrvdwly0zgb

  use, intrinsic :: ISO_C_Binding, only : C_F_Pointer , f64 => C_DOUBLE &
        , c_ptr , i64 => C_INT64_T
  implicit none

  contains

  !........................................
  subroutine bind_c_assemble_mass_1d(ne, degree, bound_spans, &
        bound_spans_shape_1, bound_spans_stride_1, bound_basis, &
        bound_basis_shape_1, bound_basis_shape_2, bound_basis_shape_3, &
        bound_basis_shape_4, bound_basis_stride_1, bound_basis_stride_2 &
        , bound_basis_stride_3, bound_basis_stride_4, bound_weights, &
        bound_weights_shape_1, bound_weights_shape_2, &
        bound_weights_stride_1, bound_weights_stride_2, bound_points, &
        bound_points_shape_1, bound_points_shape_2, &
        bound_points_stride_1, bound_points_stride_2, bound_matrix, &
        bound_matrix_shape_1, bound_matrix_shape_2, &
        bound_matrix_stride_1, bound_matrix_stride_2) bind(c)

    implicit none

    integer(i64), value :: ne
    integer(i64), value :: degree
    type(c_ptr), value :: bound_spans
    integer(i64), value :: bound_spans_shape_1
    integer(i64), value :: bound_spans_stride_1
    type(c_ptr), value :: bound_basis
    integer(i64), value :: bound_basis_shape_1
    integer(i64), value :: bound_basis_shape_2
    integer(i64), value :: bound_basis_shape_3
    integer(i64), value :: bound_basis_shape_4
    integer(i64), value :: bound_basis_stride_1
    integer(i64), value :: bound_basis_stride_2
    integer(i64), value :: bound_basis_stride_3
    integer(i64), value :: bound_basis_stride_4
    type(c_ptr), value :: bound_weights
    integer(i64), value :: bound_weights_shape_1
    integer(i64), value :: bound_weights_shape_2
    integer(i64), value :: bound_weights_stride_1
    integer(i64), value :: bound_weights_stride_2
    type(c_ptr), value :: bound_points
    integer(i64), value :: bound_points_shape_1
    integer(i64), value :: bound_points_shape_2
    integer(i64), value :: bound_points_stride_1
    integer(i64), value :: bound_points_stride_2
    type(c_ptr), value :: bound_matrix
    integer(i64), value :: bound_matrix_shape_1
    integer(i64), value :: bound_matrix_shape_2
    integer(i64), value :: bound_matrix_stride_1
    integer(i64), value :: bound_matrix_stride_2
    integer(i64), pointer :: spans(:)
    real(f64), pointer :: basis(:,:,:,:)
    real(f64), pointer :: weights(:,:)
    real(f64), pointer :: points(:,:)
    real(f64), pointer :: matrix(:,:)

    call C_F_Pointer(bound_spans, spans, [bound_spans_shape_1 * &
          bound_spans_stride_1])
    call C_F_Pointer(bound_basis, basis, [bound_basis_shape_4 * &
          bound_basis_stride_4,bound_basis_shape_3 * &
          bound_basis_stride_3,bound_basis_shape_2 * &
          bound_basis_stride_2,bound_basis_shape_1 * &
          bound_basis_stride_1])
    call C_F_Pointer(bound_weights, weights, [bound_weights_shape_2 * &
          bound_weights_stride_2,bound_weights_shape_1 * &
          bound_weights_stride_1])
    call C_F_Pointer(bound_points, points, [bound_points_shape_2 * &
          bound_points_stride_2,bound_points_shape_1 * &
          bound_points_stride_1])
    call C_F_Pointer(bound_matrix, matrix, [bound_matrix_shape_2 * &
          bound_matrix_stride_2,bound_matrix_shape_1 * &
          bound_matrix_stride_1])
    call assemble_mass_1D(ne = ne, degree = degree, spans = spans(1_i64 &
          ::bound_spans_stride_1), basis = basis(1_i64:: &
          bound_basis_stride_4, 1_i64::bound_basis_stride_3, 1_i64:: &
          bound_basis_stride_2, 1_i64::bound_basis_stride_1), weights = &
          weights(1_i64::bound_weights_stride_2, 1_i64:: &
          bound_weights_stride_1), points = points(1_i64:: &
          bound_points_stride_2, 1_i64::bound_points_stride_1), matrix &
          = matrix(1_i64::bound_matrix_stride_2, 1_i64:: &
          bound_matrix_stride_1))

  end subroutine bind_c_assemble_mass_1d
  !........................................

end module bind_c_mod_pvrvdwly0zgb
