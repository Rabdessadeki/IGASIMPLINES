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


            lvalues_ux[ : , : ]  = 0.0
            lvalues_u[ : , : ]   = 0.0
            lvalues_uy[ : , : ]  = 0.0
            lcoeffs_u[ :, : ]   = vector_u[i_span_1 : i_span_1 + p1 + 1 ,  i_span_2 : i_span_2 + p2 +1 ]
         
            for il_1 in range(0, p1+1):
            	for il_2 in range(0, p2+1):
            	
                    coef = lcoeffs_u[ il_1 , il_2 ]
                    
                    for ik1 in range(0, k1):
                    	for ik2 in range(0, k1):
                    	
	                        b0  = basis_1[ie1, il_1, 0, ik1] * basis_2[ie2, il_2, 0, ik2]
	                        dbx = basis_1[ie1, il_1, 1, ik1] * basis_2[ie2, il_2, 0, ik2]
	                        dby = basis_1[ie1, il_1, 0, ik1] * basis_2[ie2, il_2, 1, ik2]

                        	lvalues_u [ik1, ik2] += coef * b0
                        	lvalues_ux[ik1, ik2] += coef * dbx
                        	lvalues_uy[ik1, ik2] += coef * dby

            v = 0.0
            w = 0.0
            for ik1 in range(0, k1):
            	for ik2 in range(0, k2):
                    wvol = weights_1[ie1, ik1] *weights_2[ie2, ik2] 
                    x    = points_1[ie1, ik1]
                    y    = points_2[ie2, ik2]

                    # ... test 0
                    '''
                    u    = x*y*(x-1)*(y-1)
                    ux   =  x*y*(y - 1) + y*(x - 1)*(y - 1)
                    uy   = x*y*(x - 1) + x*(x - 1)*(y - 1)
                    '''
                    u     = sin( pi * x ) * sin( pi * y)
                    ux    = pi  * cos(pi * x) * sin(pi * y)
                    uy    = pi  * sin(pi * x) * cos(pi * y)
                    
                    
                    uh  = lvalues_u [ik1, ik2]
                    uhx = lvalues_ux[ik1, ik2]
                    uhy = lvalues_uy[ik1, ik2]

                    v  += ( u - uh ) ** 2 * wvol
                    w  += ( ( ux - uhx ) **2 + (uy -uhy ) ** 2 ) * wvol
            norm_l2 += v
            norm_H1 += w

    norm_l2 = sqrt(norm_l2)
    norm_H1 = sqrt(norm_H1)

    rhs[p1 , p2]   = norm_l2
    rhs[p1, p2+1] = norm_H1
    # ...
