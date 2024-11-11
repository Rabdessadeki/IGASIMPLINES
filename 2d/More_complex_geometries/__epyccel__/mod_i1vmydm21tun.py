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
