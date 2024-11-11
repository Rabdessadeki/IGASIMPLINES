__all__ = ['assemble_mass_1D',
           'assemble_stiffness_1D',
           'assemble_vector_ex01',
           'assemble_norm_ex01'
]

from pyccel.decorators import types


# assembles mass matrix 1D
#==============================================================================
@types('int', 'int', 'int[:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')

def assemble_mass_1D(ne, degree, spans, basis, weights, points, matrix):

    # ... sizes
    k1 = weights.shape[1]
    # ... build matrices
    for ie1 in range(0, ne):
            i_span_1 = spans[ie1]        
            # evaluation dependant uniquement de l'element

            for il_1 in range(0, degree+1):
            
                i1 = i_span_1 - degree + il_1
                
                for il_2 in range(0, degree+1):
                
                            i2 = i_span_1 - degree + il_2
                            v = 0.0
                            
                            for g1 in range(0, k1):
                                
                                    bi_0 = basis[ie1, il_1, 0, g1]
                                    bj_0 = basis[ie1, il_2, 0, g1]
                                    
                                    bi_x = basis[ie1, il_1, 1, g1]
                                    bj_x = basis[ie1, il_2, 1, g1]
                                    
                                    wvol = weights[ie1, g1]
                                    
                                    v +=  bi_0 * bj_0 * wvol

                            matrix[ degree + i1 , degree + i2 - i1 ]  += v
    # ... end of mass assembling

    
# assemble stiffness matrix 1D
@types('int', 'int', 'int[:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')
def assemble_stiffness_1D(ne, degree, spans, basis, weights, points, matrix):

    # ... sizes
    k1 = weights.shape[1]
    # ... build matrices
    for ie1 in range(0, ne):
            i_span_1 = spans[ie1]        
            # evaluation dependant uniquement de l'element

            for il_1 in range(0, degree+1):
            
                i1 = i_span_1 - degree + il_1
                
                for il_2 in range(0, degree+1):
                
                            i2 = i_span_1 - degree + il_2
                            v = 0.0
                            
                            for g1 in range(0, k1):                                
                               
                                    bi_x = basis[ie1, il_1, 1, g1]
                                    bj_x = basis[ie1, il_2, 1, g1]
                                    
                                    wvol = weights[ie1, g1]
                                    
                                    v +=  bi_x * bj_x * wvol

                            matrix[ degree + i1 , degree + i2 - i1 ]  += v
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

#================================ Assemble rhs  ==============================================

#---1 : In uniform mesh
@types('int', 'int', 'int', 'int','int', 'int', 'int', 'int', 'int[:]', 'int[:]','int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'real[:]', 'real[:]', 'real[:]', 'real[:]', 'double[:,:]', 'real', 'real', 'int', 'double[:,:]')
def assemble_vector_ex01(ne1, ne2, ne3, ne4, p1, p2, p3, p4, spans_1, spans_2,  spans_3, spans_4, basis_1, basis_2, basis_3, basis_4, weights_1, weights_2, weights_3, weights_4, points_1, points_2, points_3, points_4, knots_1, knots_2, knots_3, knots_4, vector_d, ovlp_value, S_DDM, domain_nb, rhs):

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

			
    lcoeffs_d   = zeros((p1+1,p2+1))
    lvalues_u      = zeros(k2)
    lvalues_ux     = zeros(k2)
    #---Computes All basis in a new points
    nders          = 1
    degree         = p3
    #..
    left           = empty( degree )
    right          = empty( degree )
    a              = empty( (       2, degree+1) )
    ndu            = empty( (degree+1, degree+1) )
    ders           = zeros( (     nders+1, degree+1) ) # output array
    #...
    basis3         = zeros((nders+1, degree+1))
    for i in range(1):
            #span = find_span( knots, degree, xq )
            xq = ovlp_value
            #~~~~~~~~~~~~~~~
            #span = find_span( knots, degree, xq )
            #~~~~~~~~~~~~~~~
            # Knot index at left/right boundary
            low  = degree
            high = len(knots_3)-1-degree
            # Check if point is exactly on left/right boundary, or outside domain
            if xq <= knots_3[low ]: 
                 span = low
            if xq >= knots_3[high]: 
                 span = high-1
            else :
              # Perform binary search
              span = (low+high)//2
              while xq < knots_3[span] or xq >= knots_3[span+1]:
                 if xq < knots_3[span]:
                     high = span
                 else:
                     low  = span
                 span = (low+high)//2
            ndu[0,0] = 1.0
            for j in range(0,degree):
                left [j] = xq - knots_3[span-j]
                right[j] = knots_3[span+1+j] - xq
                saved    = 0.0
                for r in range(0,j+1):
                    # compute inverse of knot differences and save them into lower triangular part of ndu
                    ndu[j+1,r] = 1.0 / (right[r] + left[j-r])
                    # compute basis functions and save them into upper triangular part of ndu
                    temp       = ndu[r,j] * ndu[j+1,r]
                    ndu[r,j+1] = saved + right[r] * temp
                    saved      = left[j-r] * temp
                ndu[j+1,j+1] = saved	

            # Compute derivatives in 2D output array 'ders'
            ders[0,:] = ndu[:,degree]
            for r in range(0,degree+1):
                s1 = 0
                s2 = 1
                a[0,0] = 1.0
                for k in range(1,nders+1):
                    d  = 0.0
                    rk = r-k
                    pk = degree-k
                    if r >= k:
                       a[s2,0] = a[s1,0] * ndu[pk+1,rk]
                       d = a[s2,0] * ndu[rk,pk]
                    j1 = 1   if (rk  > -1 ) else -rk
                    j2 = k-1 if (r-1 <= pk) else degree-r
                    for ij in range(j1,j2+1):
                        a[s2,ij] = (a[s1,ij] - a[s1,ij-1]) * ndu[pk+1,rk+ij]
                    for ij in range(j1,j2+1):
                        d += a[s2,ij]* ndu[rk+ij,pk]
                    if r <= pk:
                       a[s2,k] = - a[s1,k-1] * ndu[pk+1,r]
                       d += a[s2,k] * ndu[r,pk]
                    ders[k,r] = d
                    j  = s1
                    s1 = s2
                    s2 = j
            # Multiply derivatives by correct factors
            r = degree
            ders[1,:] = ders[1,:] * r
            basis3[0,:] = ders[0,:]
            basis3[1,:] = ders[1,:]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Assembles Neumann Condition
    if domain_nb == 0 :
      i1_ovrlp  = ne1+2*p1-1
      neum_sign = 1.
    else :
      i1_ovrlp  = p1
      neum_sign = -1.
                
    for ie2 in range(0, ne2):
           i_span_2 = spans_2[ie2]
           
           for g2 in range(0, k2):
                    #... We compute firstly the span in new adapted points                             
                    #...                    
                    xq        = ovlp_value
                    degree    = p3
                    low       = degree
                    high      = len(knots_3)-1-degree
                    # Check if point is exactly on left/right boundary, or outside domain
                    if xq <= knots_3[low ]: 
                         span = low
                    if xq >= knots_3[high]: 
                         span = high-1
                    else :
                      # Perform binary search
                      span = (low+high)//2
                      while xq < knots_3[span] or xq >= knots_3[span+1]:
                         if xq < knots_3[span]:
                             high = span
                         else:
                             low  = span
                         span = (low+high)//2
                    span_3    = span 
                    # ...
                    lcoeffs_d[:, :] = vector_d[span_3 : span_3+p3+1, i_span_2 : i_span_2+p4+1]
                    # ...
                    u  = 0.
                    ux  = 0.
                    for il_1 in range(0, p3+1):
                       for il_2 in range(0, p4+1):
                             bi_0      = basis3[0,il_1] * basis_4[ie2, il_2, 0, g2]
                             bi_x      = basis3[1,il_1] * basis_4[ie2, il_2, 0, g2]
                             # ...
                             coeff_d   = lcoeffs_d[il_1, il_2]
                             # ...
                             u        +=  coeff_d*bi_0
                             ux       +=  coeff_d*bi_x
                    lvalues_u[g2]  = u
                    lvalues_ux[g2] = ux
                                                
           for il_2 in range(0, p2+1):
                 i2 = i_span_2 - p2 + il_2

                 v = 0.0
                 for g2 in range(0, k2):
                       bi_0     =  basis_2[ie2, il_2, 0, g2]
                       wsurf    =  weights_2[ie2, g2]
                       ux = lvalues_ux[g2]
                       u = lvalues_u[g2]
                       #.. 
                       v       += bi_0 * (neum_sign*ux+S_DDM*u) * wsurf
                      
                 rhs[i1_ovrlp,i2+p2] += v
###################################################################################
#---1 : In uniform mesh
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]','int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]','double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'real[:]', 'real[:]', 'real[:]', 'double[:,:]', 'real', 'double[:]')
def assemble_vector_ex02(ne1, ne2, ne3, p1, p2, p3, spans_1, spans_2,  spans_3, basis_1, basis_2, basis_3, weights_1, weights_2, weights_3, points_1, points_2, points_3,  knots_1, knots_2, knots_3, vector_d, ovlp_value, rhs):

    from numpy import exp
    from numpy import pi
    from numpy import sin
    from numpy import arctan2
    from numpy import cos, cosh
    from numpy import sqrt
    from numpy import zeros
    from numpy import empty

    # ... sizes
    k1        = weights_1.shape[1]
    k2        = weights_2.shape[1]
    # ...
    # ...
    lcoeffs_d   = zeros((p2+1,p3+1))
    # ..
    lvalues_u   = zeros(k1)
    #---Computes All basis in a new points
    nders          = 0
    degree         = p2
    #..
    left           = empty( degree )
    right          = empty( degree )
    a              = empty( (       2, degree+1) )
    ndu            = empty( (degree+1, degree+1) )
    ders           = zeros( (     nders+1, degree+1) ) # output array
    #...
    basis2         = zeros(degree+1)
    for i in range(1):
            #span = find_span( knots, degree, xq )
            xq = ovlp_value
            #~~~~~~~~~~~~~~~
            # Knot index at left/right boundary
            low  = degree
            high = len(knots_2)-1-degree
            # Check if point is exactly on left/right boundary, or outside domain
            if xq <= knots_2[low ]: 
                 span = low
            if xq >= knots_2[high]: 
                 span = high-1
            else :
              # Perform binary search
              span = (low+high)//2
              while xq < knots_2[span] or xq >= knots_2[span+1]:
                 if xq < knots_2[span]:
                     high = span
                 else:
                     low  = span
                 span = (low+high)//2
            ndu[0,0] = 1.0
            for j in range(0,degree):
                left [j] = xq - knots_2[span-j]
                right[j] = knots_2[span+1+j] - xq
                saved    = 0.0
                for r in range(0,j+1):
                    # compute inverse of knot differences and save them into lower triangular part of ndu
                    ndu[j+1,r] = 1.0 / (right[r] + left[j-r])
                    # compute basis functions and save them into upper triangular part of ndu
                    temp       = ndu[r,j] * ndu[j+1,r]
                    ndu[r,j+1] = saved + right[r] * temp
                    saved      = left[j-r] * temp
                ndu[j+1,j+1] = saved	

            # Compute derivatives in 2D output array 'ders'
            ders[0,:] = ndu[:,degree]
            for r in range(0,degree+1):
                s1 = 0
                s2 = 1
                a[0,0] = 1.0
                for k in range(1,nders+1):
                    d  = 0.0
                    rk = r-k
                    pk = degree-k
                    if r >= k:
                       a[s2,0] = a[s1,0] * ndu[pk+1,rk]
                       d = a[s2,0] * ndu[rk,pk]
                    j1 = 1   if (rk  > -1 ) else -rk
                    j2 = k-1 if (r-1 <= pk) else degree-r
                    for ij in range(j1,j2+1):
                        a[s2,ij] = (a[s1,ij] - a[s1,ij-1]) * ndu[pk+1,rk+ij]
                    for ij in range(j1,j2+1):
                        d += a[s2,ij]* ndu[rk+ij,pk]
                    if r <= pk:
                       a[s2,k] = - a[s1,k-1] * ndu[pk+1,r]
                       d += a[s2,k] * ndu[r,pk]
                    ders[k,r] = d
                    j  = s1
                    s1 = s2
                    s2 = j
            basis2[:] = ders[0,:]
    # ... build rhs
    for ie1 in range(0, ne1):
            i_span_1 = spans_1[ie1]

            for g1 in range(0, k1):

                    #... We compute firstly the span in new adapted points                             
                    #...                    
                    xq        = ovlp_value
                    degree    = p2
                    low       = degree
                    high      = len(knots_2)-1-degree
                    # Check if point is exactly on left/right boundary, or outside domain
                    if xq <= knots_2[low ]: 
                         span = low
                    if xq >= knots_2[high]: 
                         span = high-1
                    else :
                      # Perform binary search
                      span = (low+high)//2
                      while xq < knots_2[span] or xq >= knots_2[span+1]:
                         if xq < knots_2[span]:
                             high = span
                         else:
                             low  = span
                         span = (low+high)//2
                    span_2    = span 
                    # ...
                    lcoeffs_d[:, : ] = vector_d[span_2 : span_2+p2+1, i_span_1 : i_span_1+p3+1]
                    # ...
                    u  = 0.
                    for il_1 in range(0, p2+1):
                       for il_2 in range(0, p3+1):
                             bi_0    = basis2[il_1] * basis_3[ie1, il_2, 0, g1]
                             # ...
                             coeff_d = lcoeffs_d[il_1, il_2]
                             # ...
                             u      +=  coeff_d*bi_0
                    lvalues_u[g1] = u

            for il_1 in range(0, p1+1):

                      i1 = i_span_1 - p1 + il_1

                      v = 0.0
                      for g1 in range(0, k1):

                              bi_0 = basis_1[ie1, il_1, 0, g1]
                              # ...
                              wvol = weights_1[ie1, g1]
                              # ...
                              u    = lvalues_u[g1]
                              # ...        
                              v   += bi_0 * u * wvol

                      rhs[i1+p1] += v 
                             
    # ...    
#==============================================================================Assemble l2 and H1 error norm

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
