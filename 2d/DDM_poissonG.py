from simplines import compile_kernel, apply_dirichlet

from simplines import SplineSpace
from simplines import TensorSpace
from simplines import StencilMatrix
from simplines import StencilVector
from simplines import pyccel_sol_field_2d
from simplines import StencilVectorSpace

# .. Matrices in 1D ..



#---In Poisson equation
from gallery_section_04 import assemble_mass_1D 
from gallery_section_04 import assemble_vector_ex001  
from gallery_section_04 import assemble_stiffness_1D
from gallery_section_04 import assemble_vector_ex01    
from gallery_section_04 import assemble_norm_ex01      
from gallery_section_04 import assemble_vector_ex02     


assemble_mass          = compile_kernel(assemble_mass_1D, arity=2)
assemble_stiffness     = compile_kernel(assemble_stiffness_1D, arity=2)
Proj_sol               = compile_kernel(assemble_vector_ex02, arity=1)
assemble_rhs           = compile_kernel(assemble_vector_ex01, arity=1)
assemble_rhs01         = compile_kernel(assemble_vector_ex001, arity=1)
assemble_norm_l2       = compile_kernel(assemble_norm_ex01, arity=1)

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

       
		
class Poisson_2d(object):
	def __init__(self, V1, V2, Vt):
		K1                   =  assemble_stiffness(V1)
		K1                   =  K1.tosparse()
		K1                   =  K1.toarray()[1:-1,1:-1]
		K1                   =  csr_matrix(K1)

		M1                   =  assemble_mass(V1)
		M1                   =  M1.tosparse()
		M1                   =  M1.toarray()[1:-1,1:-1]
		M1                   =  csr_matrix(M1)

		K2                   =  assemble_stiffness(V2)
		K2                   =  K2.tosparse()
		K2                   =  K2.toarray()[1:-1,1:-1]
		K2                   =  csr_matrix(K2)
			

		M2                   =  assemble_mass(V2)
		M2                   =  M2.tosparse()
		M2                   =  M2.toarray()[1:-1,1:-1]
		M2                   =  csr_matrix(M2) 
		
       	
		M                    =  kron(K1, M2) + kron(M1 , K2)
		M                    =  csr_matrix(M)
	
		self.lu              =  sla.splu(M)
		self.V2              =  V2
		self.V1              =  V1
		self.Vt              =  Vt	
		self.nbasis_1        =  (V1.nbasis-2)*(V2.nbasis-2)
		self.nbasis_2        =  [V1.nbasis-2, V2.nbasis-2]
	def solve(self, u_d):
		rhs                  = StencilVector(self.Vt.vector_space)
		rhs                  = assemble_rhs01( self.Vt, out = rhs) 
		b                    = rhs.toarray() 	 
		b                    = b.reshape(self.Vt.nbasis)
		b                    = b[ 1 : -1 , 1 : -1]
		b                    = b.reshape(self.nbasis_1)

		xkron                = self.lu.solve(b)       
		xkron                = xkron.reshape(self.nbasis_2)
		x                    = np.zeros(self.Vt.nbasis)
		x[ 1 : -1 , 1 : -1]  = xkron
	

		#... Dirichlet nboundary
		u                    = StencilVector(self.Vt.vector_space)
		u.from_array(self.Vt, x)
		# ...
		Norm                 = assemble_norm_l2(self.Vt, fields=[u]) 
		norm                 = Norm.toarray()
		
		l2_norm              = norm[0]
		H1_norm              = norm[1]         
		return u, x, l2_norm, H1_norm
	   		

degree      = 6
nelements   = 64
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
V_1          = SplineSpace(degree = degree, nelements = nelements, nderiv = 2, quad_degree = quad_degree)
V_2            = SplineSpace(degree = degree, nelements = nelements, nderiv = 2, quad_degree = quad_degree)
V_12         = TensorSpace(V_1, V_2)
P = Poisson_2d(V_1, V_2, V_12)	
u_00    = StencilVector(V_12.vector_space)
u_11,   xuh0, l2_norm0, H1_norm0     = P.solve( u_d = u_00)
print('the solution without DDM ')
print('iteration  -----> L^2-error ={} -----> H^1-error = {}'.format( l2_norm0, H1_norm0))

#

nbpts = 120
fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(121, projection='3d')
u, ux, uy, X, Y               = pyccel_sol_field_2d((nbpts, nbpts),  xuh0,   V_12.knots, V_12.degree)
# ...
# plot a 3D surface like in the example mplot3d/surface3d_demo

surf0 =		ax.plot_wireframe(X[:,:], Y[:,:], u[:,:], color ='black') 

#ax.set_xlim(0.0, 1.0)
#ax.set_ylim(0.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.set_title('Approximate Solution in adaptive meshes')
ax.set_xlabel('F1',  fontweight ='bold')
ax.set_ylabel('F2',  fontweight ='bold')
#fig.plot(surf0, shrink=0.5, aspect=25)
plt.grid()
Sol = lambda x,y : sin(pi*x)* sin(pi*y)
ax = fig.add_subplot(122, projection='3d')
surf0 =		ax.plot_wireframe(X[:,:], Y[:,:], Sol(X, Y), color ='black') 

#ax.set_xlim(0.0, 1.0)
#ax.set_ylim(0.0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.set_title('Approximate Solution in adaptive meshes')
ax.set_xlabel('X',  fontweight ='bold')
ax.set_ylabel('Y',  fontweight ='bold')
plt.savefig('Poisson3D.png')
plt.show()
 

