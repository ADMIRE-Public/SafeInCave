import safeincave as sf
import torch as to

def main():
	n_elems = 1
	mat = sf.Material(n_elems)

	# Elastic element
	E = 102e9*to.ones(n_elems)
	nu = 0.3*to.ones(n_elems)
	spring_0 = sf.Spring(E, nu, "spring")

	# Viscoelastic element
	eta = 105e11*to.ones(n_elems)
	E = 10e9*to.ones(n_elems)
	nu = 0.32*to.ones(n_elems)
	kelvin = sf.Viscoelastic(eta, E, nu, "kelvin")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)


if __name__ == '__main__':
	main()