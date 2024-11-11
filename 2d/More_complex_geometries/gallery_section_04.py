__all__ = ['assemble_matrix_ex01',
           'assemble_vector_ex01',
           'assemble_norm_ex01'
]

from pyccel.decorators import types


# assembles matrix 

@types('int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:,:,:]')
def assemble_matrix_ex01(ne1, ne2, p1, p2, spans_1, spans_2, basis_1, basis_2, weights_1, weights_2, points_1, points_2, matrix):
	from numpy import exp
	from numpy import cos
	from numpy import sin
	from numpy import pi
	from numpy import sqrt
	from numpy import zeros

	k1 = weights_1.shape[1]
	k2 = weights_2.shape[1]

	#.. circle

	# ...
	lcoeffs_u  = zeros((p1+1,p2+1))
	arr_J_mat0 = zeros((k1,k2))
	arr_J_mat1 = zeros((k1,k2))
	arr_J_mat2 = zeros((k1,k2))
	arr_J_mat3 = zeros((k1,k2))
	J_mat      = zeros((k1,k2))

    # ... build matrices
	for ie1 in range(0, ne1):
		i_span_1 = spans_1[ie1]
		for ie2 in range(0, ne2):
			i_span_2 = spans_2[ie2]

			for g1 in range(0, k1):
				for g2 in range(0, k2):

					x1       = points_1[ie1,g1]
					x2       = points_2[ie2,g2]


					F1x      = 2.0*sqrt(1.0-0.5*(2.0*x2-1.0)**2)
					F1y      = -(2.0*x1-1.0)*(2.0*x2-1.0)/sqrt(1.0-0.5*(2.0*x2-1.0)**2)
					F2y      = 2.0*sqrt(1.0-0.5*(2.0*x1-1.0)**2)
					F2x      = -(2.0*x1-1.0)*(2.0*x2-1.0)/sqrt(1.0-0.5*(2.0*x1-1.0)**2)
					det_Hess = abs(F1x*F2y-F1y*F2x)

					#F=(F1, F2)
					arr_J_mat0[g1,g2] = F2y
					arr_J_mat1[g1,g2] = F1x
					arr_J_mat2[g1,g2] = F1y
					arr_J_mat3[g1,g2] = F2x
			
					J_mat[g1,g2]      = det_Hess
		                    
			for il_1 in range(0, p1+1):
				for il_2 in range(0, p2+1):
					for jl_1 in range(0, p1+1):
				   		for jl_2 in range(0, p2+1):
				   			i1 = i_span_1 - p1 + il_1
				   			j1 = i_span_1 - p1 + jl_1
				   			i2 = i_span_2 - p2 + il_2
				   			j2 = i_span_2 - p2 + jl_2
				   			v  = 0.0
				   			for g1 in range(0, k1):
				   				for g2 in range(0, k2):
				   					
				   					bi_x1 = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2]
				   					bi_x2 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2]
				   					
				   					bj_x1 = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2]
				   					bj_x2 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2]
				   					
				   					bi_x = arr_J_mat0[g1,g2] * bi_x1 - arr_J_mat3[g1,g2] * bi_x2
				   					bi_y = arr_J_mat1[g1,g2] * bi_x2 - arr_J_mat2[g1,g2] * bi_x1
				   					
				   					bj_x = arr_J_mat0[g1,g2] * bj_x1 - arr_J_mat3[g1,g2] * bj_x2
				   					bj_y = arr_J_mat1[g1,g2] * bj_x2 - arr_J_mat2[g1,g2] * bj_x1
				   					
				   					wvol = weights_1[ie1, g1] * weights_2[ie2, g2]
				   					
				   					v   += (bi_x * bj_x + bi_y * bj_y ) * wvol / J_mat[g1,g2]
				   			matrix[p1+i1, p2+i2, p1+j1-i1, p2+j2-i2]  += v
								
								    
							
						
@types('int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')
def assemble_vector_ex01(ne1, ne2, p1, p2, spans_1, spans_2, basis_1, basis_2, weights_1, weights_2, points_1, points_2, rhs):

	from numpy import exp
	from numpy import cos
	from numpy import sin
	from numpy import pi
	from numpy import sqrt
	from numpy import zeros, empty

	k1 = weights_1.shape[1]
	k2 = weights_2.shape[1]
	J_mat      = zeros((k1,k2))

	for ie1 in range(ne1):
		for ie2 in range(ne2):
			i_spans1 = spans_1[ ie1 ]
			i_spans2 = spans_2[ ie2 ]
			for g1 in range(0, k1):
				for g2 in range(0, k2):
				
					x              =  points_1[ie1, g1]
					y              =  points_2[ie2, g2]
					F1x      = 2.0*sqrt(1.0-0.5*(2.0*y-1.0)**2)
					F1y      = -(2.0*x-1.0)*(2.0*y-1.0)/sqrt(1.0-0.5*(2.0*y-1.0)**2)
					F2y      = 2.0*sqrt(1.0-0.5*(2.0*x-1.0)**2)
					F2x      = -(2.0*x-1.0)*(2.0*y-1.0)/sqrt(1.0-0.5*(2.0*x-1.0)**2)
					det_jac        = (F1x * F2y  - F2x * F1y)
					x1             = (2.0*x-1.0)*sqrt(1.0-0.5*(2.0*y-1.0)**2)
					x2             = (2.0*y-1.0)*sqrt(1.0-0.5*(2.0*x-1.0)**2)
					f              = 4.
					J_mat [g1 , g2] = det_jac * f
			for il1 in range( p1 + 1 ):
				i1 = i_spans1 - p1 + il1
				for il2 in range( p2 + 1):
					i2 = i_spans2 - p2 + il2
					v = 0.0
					for g1 in range(k1):
						for g2 in range(k2):
							bx_0           = basis_1[ie1, il1, 0, g1]* basis_2[ie2, il2, 0, g2]	
							wvol           = weights_1[ ie1 , g1] * weights_2[ ie2 , g2]
					
							v+= bx_0 * wvol * J_mat[g1, g2]
					rhs[i1 + p1 , i2 + p2 ]+=v

#================================ Assemble rhs  ==============================================


@types('int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]',  'double[:,:]','double[:,:]')
def assemble_norm_ex01(ne1, ne2, p1, p2, spans_1, spans_2,  basis_1, basis_2,  weights_1, weights_2, points_1, points_2, vector_u, rhs):

	from numpy import exp
	from numpy import cos
	from numpy import sin
	from numpy import pi
	from numpy import sqrt
	from numpy import zeros

	# ... sizes
	k1 = weights_1.shape[1]
	k2 = weights_2.shape[1]
	# ...

	lcoeffs_u = zeros((p1+1,p2+1))
	lvalues_u = zeros((k1, k2))
	lvalues_ux = zeros((k1, k2))
	lvalues_uy = zeros((k1, k2))
	norm_l2 = 0.
	norm_H1 = 0.
	for ie1 in range(0, ne1):
		for ie2 in range(ne2):
			i_span_1 = spans_1[ie1]
			i_span_2 = spans_1[ie2]
			lcoeffs_u[ :, : ]   = vector_u[i_span_1 : i_span_1 + p1 + 1 ,  i_span_2 : i_span_2 + p2 +1 ]
			lvalues_u[:,:]   = 0.0
			lvalues_ux[:,:]   = 0.0
			lvalues_uy[:,:]   = 0.0
			for il1 in range(p1 + 1):
				for il2 in range(p2+ 1):
					coef = lcoeffs_u[il1, il2]
					for ik1 in range(0, k1):
						for ik2 in range(0, k2):
							bi00 = basis_1[ie1 , il1, 0, ik1] * basis_2[ie2 , il2, 0, ik2]
							bi10 = basis_1[ie1 , il1, 1, ik1] * basis_2[ie2 , il2, 0, ik2]
							bi01 = basis_1[ie1 , il1, 0, ik1] * basis_2[ie2 , il2, 1, ik2]

							lvalues_u[ik1 , ik2]  += coef * bi00
							lvalues_ux[ik1 , ik2] += coef * bi10
							lvalues_uy[ik1 , ik2] += coef * bi01
													
						
			for ik1 in range(0, k1):
				for ik2 in range(0, k2):
					wvol = weights_1[ie1, ik1] *weights_2[ie2, ik2]
					x              =  points_1[ie1, ik1]	
					y              =  points_2[ie2, ik2]

					F1x      = 2.0*sqrt(1.0-0.5*(2.0*y-1.0)**2)
					F1y      = -(2.0*x-1.0)*(2.0*y-1.0)/sqrt(1.0-0.5*(2.0*y-1.0)**2)
					F2y      = 2.0*sqrt(1.0-0.5*(2.0*x-1.0)**2)
					F2x      = -(2.0*x-1.0)*(2.0*y-1.0)/sqrt(1.0-0.5*(2.0*x-1.0)**2)

					x1             = (2.0*x-1.0)*sqrt(1.0-0.5*(2.0*y-1.0)**2)
					x2             = (2.0*y-1.0)*sqrt(1.0-0.5*(2.0*x-1.0)**2)

					det_jac        = (F1x * F2y  - F2x * F1y)
					u              = -x1**2 - x2**2 + 1
					ux             = -2 * x1
					uy             = -2 * x2
					
					uh             = lvalues_u[ik1, ik2]
					uhx0            = lvalues_ux[ik1, ik2]
					uhy0            = lvalues_uy[ik1, ik2]
					
					uhx            = (uhx0 * F2y - uhy0 * F2x) / det_jac
					uhy            = (-uhx0 * F1y + uhy0 * F1x) / det_jac

					norm_l2  += ( u - uh ) ** 2 * wvol * det_jac
					norm_H1  += ( ( ux - uhx) **2  + (uy - uhy  ) ** 2 ) * wvol * det_jac
			

	norm_l2 = sqrt(norm_l2)
	norm_H1 = sqrt(norm_H1)

	rhs[p1 , p2]   = norm_l2
	rhs[p1, p2+1] = norm_H1
	# ...
