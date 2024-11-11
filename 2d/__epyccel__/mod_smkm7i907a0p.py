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
