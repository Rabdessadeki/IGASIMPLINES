




'''
from numpy import random
n = 4
A = random.randint(10,size = (n, n, n, n))
A[0, : , :, :] =0. 
A[:, 0, : ,:] = 0.

A[n-1, :, : ,:] = 0.
A[:, n-1, : ,:] = 0.
A[0,:,0,0]      = 1.
A[:,0,0,0]      = 1.
A[n-1,:,0,0]      = 1.
A[:,n-1,0,0]      = 1.


print(A)
'''


from simplines import compile_kernel, apply_dirichlet

from simplines import SplineSpace
from simplines import TensorSpace
from simplines import StencilMatrix
from simplines import StencilVector
from simplines import pyccel_sol_field_2d
from simplines import StencilVectorSpace

# .. Matrices in 1D ..



#---In Poisson equation
from gallery_section_04 import assemble_matrix_ex01 
from gallery_section_04 import assemble_vector_ex01  

from gallery_section_04 import assemble_norm_ex01      




assemble_stiffness_2d   = compile_kernel(assemble_matrix_ex01, arity=2)

assemble_rhs_2d        = compile_kernel(assemble_vector_ex01, arity=1)

assemble_norm_l2     = compile_kernel(assemble_norm_ex01, arity=1)

#from matplotlib.pyplot import plot, show
import matplotlib.pyplot            as     plt
from   mpl_toolkits.axes_grid1      import make_axes_locatable
from   mpl_toolkits.mplot3d         import axes3d
from   matplotlib                   import cm
from   mpl_toolkits.mplot3d.axes3d  import get_test_data
from   matplotlib.ticker            import LinearLocator, FormatStrFormatter
#..
from   scipy.sparse                 import kron
from   scipy.sparse                 import csr_matrix
from   scipy.sparse                 import csc_matrix, linalg as sla
from   numpy                        import zeros, linalg, asarray, linspace
from   numpy                        import cos, sin, pi, exp, sqrt, arctan2
from   tabulate                     import tabulate
import numpy                        as     np
import timeit
import time

from tabulate import tabulate

#==============================================================================
#.......Poisson ALGORITHM

class Poisson_2d(object):
   def __init__(self, V1, V2, Vt):
       
       stiffness           = assemble_stiffness_2d(Vt )
       stiffness           = apply_dirichlet(Vt, stiffness, dirichlet = [[True, True],[True, True]])
       
       M                   = stiffness.tosparse()

       self.lu             = sla.splu(csc_matrix(M))
       self.Vt             = Vt


       
   def solve(self, u_d = None ):
       u                   = StencilVector(self.Vt.vector_space)
       # ++++
       #--Assembles a right hand side of Poisson equation
       rhs                 = StencilVector(self.Vt.vector_space)
       rhs                 = assemble_rhs_2d( self.Vt, out = rhs )

      
       rhs                 = apply_dirichlet(self.Vt, rhs, dirichlet = [[True, True],[True, True]])
     
       b                   = rhs.toarray()
       # ...
       x                   = self.lu.solve(b)       
       # ...
       x                   = x.reshape(self.Vt.nbasis)
       #...
       u.from_array(self.Vt, x)
       # ...
       Norm                = assemble_norm_l2(self.Vt, fields=[u]) 
       norm                = Norm.toarray()
       l2_norm             = norm[0]
       H1_norm             = norm[1]  
       return u, x, l2_norm, H1_norm 
       
		


degree      = 2
nelements   = 128
iter_max    = 10
quad_degree = degree + 1 
alpha       = 0.5
beta        = 0.5
left_v      = 0.
right_v     = 1.
h           = 1./nelements
kmax        = np.pi/h
kmin        = np.pi
R_Const     = 1./beta#(kmin**2*kmax**2)**(1/4)

# create the spline space for each direction
V_1             = SplineSpace(degree = degree, nelements = nelements, nderiv = 2, quad_degree = quad_degree)
V_2             = SplineSpace(degree = degree, nelements = nelements, nderiv = 2, quad_degree = quad_degree)
V_12            = TensorSpace(V_1, V_2)

P = Poisson_2d(V_1, V_2, V_12)	

u_00    = StencilVector(V_12.vector_space)
u_11,   xuh0, l2_norm0, H1_norm0     = P.solve( )
print('The error without domain decompsition in the circulare domain')
print('-----> L^2-error ={} -----> H^1-error = {}'.format(l2_norm0, H1_norm0))

fig = plt.figure(figsize=(15, 10))

nbpts = 200
u, ux, uy, sX, sY               = pyccel_sol_field_2d((nbpts, nbpts),  xuh0,   V_12.knots, V_12.degree)
# ...
X = (2.0*sX-1.0) * sqrt(1.-0.5*(2.0*sY-1.0)**2)
Y = (2.0*sY-1.0) * sqrt(1.-0.5*(2.0*sX-1.0)**2)
	# ...
sol = lambda x , y  : 1 - x**2 - y**2

ax = fig.add_subplot(121, projection='3d')
surf0 =		ax.plot_wireframe(X[:,:], Y[:,:], u[:,:], color ='black') 


ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_title('Approximate Solution ')
ax.set_xlabel('F1',  fontweight ='bold')
ax.set_ylabel('F2',  fontweight ='bold')
#fig.plot(surf0, shrink=0.5, aspect=25)

plt.grid()


ax = fig.add_subplot(122, projection='3d')
surf0 =		ax.plot_wireframe(X[:,:], Y[:,:], sol(X[:,:], Y[:,:]), color ='black') 


ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_title('Exact Solution ')
ax.set_xlabel('X',  fontweight ='bold')
ax.set_ylabel('Y',  fontweight ='bold')
#fig.plot(surf0, shrink=0.5, aspect=25)
plt.savefig('Poisson3D_Cicule.png')
plt.grid()
plt.show()
 
	
'''	
	
	u_0,   xuh0, l2_norm0, H1_norm0       = DDM_0.solve( u_d = u_1 )
	
	u_1, xuh_1, l2_norm1, H1_norm1        = DDM_1.solve( u_d =  u_00 )
	xuh_0.append(xuh0)
	xuh_01.append(xuh_1)
	u_00 = u_0 
	l2_err = sqrt(l2_norm0**2 + l2_norm1**2)
	H1_err = sqrt(H1_norm0**2 + H1_norm1**2)
	L2.append(l2_err)
	H1.append(H1_err)
	


#---Compute a solution
u, ux, uy, X, Y               = pyccel_sol_field_2d((nbpts, nbpts),  xuh0,   V_01.knots, V_01.degree)
# ...
u_1, ux_1, uy_1, X_1, Y_1     = pyccel_sol_field_2d((nbpts, nbpts),  xuh_1,   V_02.knots, V_02.degree)
u_00  = []
u_01 = []
	for i in range(iter_max):
	    u_00.append(pyccel_sol_field_2d((nbpts, nbpts),  xuh_0[i],   V_01.knots, V_01.degree)[0][:,50])
	    u_01.append(pyccel_sol_field_2d((nbpts, nbpts),  xuh_01[i],   V_02.knots, V_02.degree)[0][:,50])
	   
	# solut = lambda x , y : x*y*(x-1)*(y-1)       
	solut = lambda  x, y :  sin( pi*x)* sin( pi*y)

	plt.figure() 
	#plt.axes().set_aspect('equal')
	plt.subplot(221)
	for i in range(iter_max-1):
	     plt.plot(X[:,50], u_00[i], '-k', linewidth = 1.)
	     plt.plot(X_1[:,50], u_01[i], '-k', linewidth = 1.)
	plt.plot(X_1[:,50], u_01[i+1], '-k', linewidth = 1., label='$\mathbf{Un_1-iter(i)}$')
	plt.plot(X[:,50], u_00[i+1], '-k', linewidth = 1., label='$\mathbf{Un_0-iter(i)}$')
	plt.grid(True)
	
	plt.legend()
	plt.subplot(222)
	plt.plot(X[:,50], u[:,50],  '--g', label = '$\mathbf{Un_0-iter-max}$' )
	plt.plot(X_1[:,50], u_1[:,50],  '--b', label = '$\mathbf{Un_1-iter-max}$')
	plt.grid(True)
	plt.legend()
	plt.subplot(223)
	
	plt.plot(L2, 'X', color = 'green')

	plt.plot(H1, 's', color = 'magenta')
	plt.yscale('log')
	plt.legend(['L2 norm', 'H1 norm'])
	plt.grid()
	plt.savefig('Behvoir_onXaxis_withl2normParallel.png')
	plt.show()
	
	
	nbpts = 20
	u, ux, uy, X, Y               = pyccel_sol_field_2d((nbpts, nbpts),  xuh0,   V_01.knots, V_01.degree)
	# ...
	u_1, ux_1, uy_1, X_1, Y_1     = pyccel_sol_field_2d((nbpts, nbpts),  xuh_1,   V_02.knots, V_02.degree)
	# set up a figure twice as wide as it is tall
	fig = plt.figure(figsize=plt.figaspect(0.5))
	#===============
	# First subplot
	# set up the axes for the first plot
	ax = fig.add_subplot(1, 2, 1, projection='3d')
	# plot a 3D surface like in the example mplot3d/surface3d_demo
	surf0 =		ax.plot_wireframe(X[:,:], Y[:,:], u[:,:], color ='gray') 

	surf0 = ax.plot_wireframe(X_1[:,:], Y_1[:,:], u_1[:,:], color ='green')
	ax.set_xlim(0.0, 1.0)
	ax.set_ylim(0.0, 1.0)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	#ax.set_title('Approximate solution in uniform mesh')
	ax.set_xlabel('X',  fontweight ='bold')
	ax.set_ylabel('Y',  fontweight ='bold')
	# Add a color bar which maps values to colors.
	#fig.colorbar(surf0, shrink=0.5, aspect=25)

	#===============
	# Second subplot
	ax = fig.add_subplot(1, 2, 2, projection='3d')
	
	surf = ax.plot_surface(X_1[:,:], Y_1[:,:], solut(X_1[:,:], Y_1[:,:]), cmap=cm.coolwarm, linewidth=0, antialiased=False)
	surf = ax.plot_surface(X[:,:], Y[:,:], solut(X[:,:], Y[:,:]), cmap='viridis',linewidth=0, antialiased=False)
	
    
	surf = ax.plot_wireframe(X[:,:], Y[:,:], solut(X[:,:], Y[:,:]), color ='magenta')
	surf =		ax.plot_wireframe(X_1[:,:], Y_1[:,:], solut(X_1[:,:], Y_1[:,:]), color ='magenta')
	ax.set_xlim(0.0, 1.0)
	ax.set_ylim(0.0, 1.0)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	#ax.set_title('Approximate Solution in adaptive meshes')
	ax.set_xlabel('F1',  fontweight ='bold')
	ax.set_ylabel('F2',  fontweight ='bold')
	#fig.colorbar(surf, shrink=0.5, aspect=25)
	plt.grid()

	plt.show()
'''
