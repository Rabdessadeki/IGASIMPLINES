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
