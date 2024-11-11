@types('int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')
def assemble_vector_ex001(ne1, ne2, p1, p2, spans_1, spans_2, basis_1, basis_2, weights_1, weights_2, points_1, points_2, rhs):

    from numpy import exp
    from numpy import cos
    from numpy import sin
    from numpy import pi
    from numpy import sqrt
    from numpy import zeros, empty
    
    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
   
    
    for ie1 in range(ne1):
	    for ie2 in range(ne2):
	    
		    i_spans1 = spans_1[ ie1 ]
		    i_spans2 = spans_2[ ie2 ]	
		        
		    for il1 in range( p1 + 1 ):
		    	i1 = i_spans1 - p1 + il1
		    	for il2 in range( p2 + 1):
		    		i2 = i_spans2 - p2 + il2
		    		v = 0.0
		    		for ik1 in range(k1):
		    			for ik2 in range(k2):
		    				bx_0 = basis_1[ie1, il1, 0, ik1]
		    				bx_1 = basis_1[ie1, il1, 1, ik1]
		    				by_0 = basis_2[ie2, il2, 0, ik2]
		    				by_1 = basis_2[ie2, il2, 1, ik2]
		    				wval = weights_1[ ie1 , ik1] * weights_2[ ie2 , ik2]
		    				x = points_1[ie1, ik1]
		    				y = points_2[ie2, ik2]
		    				f = 2*pi**2*sin(pi*x )*sin(pi*y)
		    				v+= (  f * by_0 * bx_0  )  * wval
		    		rhs[i1 + p1 , i2 + p2 ]+=v
