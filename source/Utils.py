import torch as to
# from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai
# from Elements import *
from dolfin import *
import json

MPa = 1e6
minute = 60
hour = 60*minute
day = 24*hour
year = 365*day

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def local_projection(tensor, V):
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(tensor, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    u = Function(V)
    solver.solve_local_rhs(u)
    return u

def epsilon(u):
	return sym(grad(u))

def dotdot(C, eps):
	return voigt2stress(dot(C, strain2voigt(eps)))

def strain2voigt(e):
	x = 1
	return as_vector([e[0,0], e[1,1], e[2,2], x*e[0,1], x*e[0,2], x*e[1,2]])

def voigt2stress(s):
    return as_matrix([	[s[0], s[3], s[4]],
						[s[3], s[1], s[5]],
						[s[4], s[5], s[2]]])

def to_tensor(numpy_array):
	return to.tensor(numpy_array, dtype=to.float64)

def dotdot2(C0_torch, eps_tot_torch):
	n_elems = C0_torch.shape[0]
	eps_tot_voigt = to.zeros((n_elems, 6), dtype=to.float64)
	eps_tot_voigt[:,0] = eps_tot_torch[:,0,0]
	eps_tot_voigt[:,1] = eps_tot_torch[:,1,1]
	eps_tot_voigt[:,2] = eps_tot_torch[:,2,2]
	eps_tot_voigt[:,3] = eps_tot_torch[:,0,1]
	eps_tot_voigt[:,4] = eps_tot_torch[:,0,2]
	eps_tot_voigt[:,5] = eps_tot_torch[:,1,2]
	stress_voigt = to.bmm(C0_torch, eps_tot_voigt.unsqueeze(2)).squeeze(2)
	stress_torch = to.zeros_like(eps_tot_torch, dtype=to.float64)
	stress_torch[:,0,0] = stress_voigt[:,0]
	stress_torch[:,1,1] = stress_voigt[:,1]
	stress_torch[:,2,2] = stress_voigt[:,2]
	stress_torch[:,0,1] = stress_torch[:,1,0] = stress_voigt[:,3]
	stress_torch[:,0,2] = stress_torch[:,2,0] = stress_voigt[:,4]
	stress_torch[:,1,2] = stress_torch[:,2,1] = stress_voigt[:,5]
	return stress_torch


def get_list_of_elements(input_model, n_elems, element_class="Elastic"):
	from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai
	ELEMENT_DICT = {
		"Spring": Spring,
		"KelvinVoigt": Viscoelastic,
		"DislocationCreep": DislocationCreep,
		"ViscoplasticDesai": ViscoplasticDesai
	}
	list_of_elements = []
	props = input_model[element_class]
	for elem_name in props.keys():
		if props[elem_name]["active"] == True:
			element_parameters = props[elem_name]["parameters"]
			for param in element_parameters:
				element_parameters[param] = element_parameters[param]*to.ones(n_elems, dtype=to.float64)
			elem = ELEMENT_DICT[props[elem_name]["type"]](element_parameters)
			list_of_elements.append(elem)
	if element_class == "Elastic" and len(list_of_elements) == 0:
		raise Exception("Model must have at least 1 elastic element (Spring). None was given.")
	return list_of_elements


def get_list_of_elements_new(input_model, n_elems, element_class="Elastic"):
	from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai
	ELEMENT_DICT = {
		"Spring": Spring,
		"KelvinVoigt": Viscoelastic,
		"DislocationCreep": DislocationCreep,
		"ViscoplasticDesai": ViscoplasticDesai
	}
	list_of_elements = []
	props = input_model[element_class]
	for elem_name in props.keys():
		if props[elem_name]["active"] == True:
			element_parameters = props[elem_name]["parameters"]
			for param in element_parameters:
				element_parameters[param] = to.tensor(element_parameters[param])
				# element_parameters[param] = element_parameters[param]*to.ones(n_elems, dtype=to.float64)
			elem = ELEMENT_DICT[props[elem_name]["type"]](element_parameters)
			list_of_elements.append(elem)
	if element_class == "Elastic" and len(list_of_elements) == 0:
		raise Exception("Model must have at least 1 elastic element (Spring). None was given.")
	return list_of_elements