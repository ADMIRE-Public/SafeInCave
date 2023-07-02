import fenics as fe
import sympy as sy
import json
import numpy as np

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
month = 30*day
kPa = 1e3
MPa = 1e6
GPa = 1e9

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def strain2voigt(e):
	x = 1
	return fe.as_vector([e[0,0], e[1,1], e[2,2], x*e[0,1], x*e[0,2], x*e[1,2]])
	# return fe.as_vector([e[0,0], e[1,1], e[2,2], e[0,1], e[0,2], e[1,2]])

def voigt2stress(s):
    return fe.as_matrix([[s[0], s[3], s[4]],
		    		     [s[3], s[1], s[5]],
		    		     [s[4], s[5], s[2]]])

def epsilon(u):
	# return 0.5*(fe.nabla_grad(u) + fe.nabla_grad(u).T)
	# return 0.5*(fe.grad(u) + fe.grad(u).T)
	return fe.sym(fe.grad(u))

def sigma(C, eps):
	return voigt2stress(fe.dot(C, strain2voigt(eps)))

def local_projection(tensor, V):
	dv = fe.TrialFunction(V)
	v_ = fe.TestFunction(V)
	a_proj = fe.inner(dv, v_)*fe.dx
	b_proj = fe.inner(tensor, v_)*fe.dx
	solver = fe.LocalSolver(a_proj, b_proj)
	solver.factorize()
	u = fe.Function(V)
	solver.solve_local_rhs(u)
	return u

def constitutive_matrix_sy(E, nu):
	lame = E*nu/((1+nu)*(1-2*nu))
	# lame = 0.0
	G = E/(2 + 2*nu)
	x = 2
	M = sy.Matrix(6, 6, [ (2*G+lame),	lame,			lame,			0.,		0.,		0.,
						  lame,			(2*G+lame),		lame,			0.,		0.,		0.,
						  lame,			lame,			(2*G+lame),		0.,		0.,		0.,
							0.,			0.,				0.,				x*G,	0.,		0.,
							0.,			0.,				0.,				0.,		x*G,	0.,
							0.,			0.,				0.,				0.,		0.,		x*G])
	return M

def double_dot(A, B):
	# Performs the operation A:B, which returns a scalar. 
	# A and B are second order tensors (2d numpy arrays)
	return np.tensordot(A, B.T, axes=2)
	# return float(np.tensordot(A, B, axes=([0, 1], [0, 1])))

ppos = lambda x: (x+abs(x))/2.