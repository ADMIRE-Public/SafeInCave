# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .MassBC import BcHandler

import dolfinx as do
import ufl
from petsc4py import PETSc
import torch as to
from .MaterialProps import Material
from .Grid import GridHandlerGMSH
from abc import ABC, abstractmethod



class CGA(object):
	def __init__(self, x_initial, step, a=1.1, b=1.5, c=2.5):
		self.x = x_initial
		self.xs = [self.x]
		self.f_xs = []
		self.initialStep = step
		self.step = self.initialStep
		self.a = a
		self.b = b
		self.c = c
		self.direction = 1.
		self.n_reverses = 0


	def compute(self, f_x):
		self.anxious_goat_point = 0
		self.__saveFunctionValue(f_x)
		if len(self.f_xs) < 2:
			self.__walk()
			return self.x, self.anxious_goat_point
		else:
			self.__peacefulGoat()
		return self.x, self.anxious_goat_point


	def __peacefulGoat(self):
		if self.f_xs[-1] > self.f_xs[-2]:
			self.__walk()
			self.n_reverses = 0
			self.__increaseRushedStep()
		else:
			self.__reverse()
			self.__angryGoat()
			self.__walk()


	def __angryGoat(self):
		if self.n_reverses >= 2:
			self.__increaseStep()
			self.anxious_goat_point = 1


	def __walk(self):
		self.x += self.step*self.direction
		self.__saveDelta()


	def __reverse(self):
		self.direction = -self.direction
		self.__decreaseStep()
		self.n_reverses += 1


	def __increaseRushedStep(self):
		self.step = min(self.initialStep, self.step*self.a)


	def __decreaseStep(self):
		self.step = self.step/self.b


	def __increaseStep(self):
		self.step = min(self.initialStep, self.step*self.b**self.c)


	def __saveDelta(self):
		self.xs.append(self.x)


	def __saveFunctionValue(self, f_x):
		if len(self.f_xs) > 3:
			self.f_xs.pop(0)
		self.f_xs.append(f_x)



class MassEquationBase(ABC):
    def __init__(self, grid: GridHandlerGMSH, is_linear: bool=True):
        self.grid = grid
        self.is_linear = is_linear

        self.create_function_spaces()

        self.n_elems = self.DG0_1.dofmap.index_map.size_local + len(self.DG0_1.dofmap.index_map.ghosts)
        self.n_nodes = self.V.dofmap.index_map.size_local + len(self.V.dofmap.index_map.ghosts)
        self.dt = do.fem.Constant(self.grid.mesh, 1.0)

        self.create_trial_test_functions()
        self.create_ds_dx()
        self.create_fenicsx_fields()

    



    def create_fenicsx_fields(self) -> None:
        self.perm = do.fem.Function(self.DG0_1)
        self.M = do.fem.Function(self.DG0_1)
        self.P_old = do.fem.Function(self.V)
        self.P = do.fem.Function(self.V)
        self.X = do.fem.Function(self.V)


    def set_solver(self, solver: PETSc.KSP) -> None:
        self.solver = solver


    def set_boundary_conditions(self, bc: BcHandler) -> None:
        self.bc = bc


    def create_trial_test_functions(self) -> None:
        self.dp = ufl.TrialFunction(self.V)
        self.p_ = ufl.TestFunction(self.V)


    def create_function_spaces(self) -> None:
        self.DG0_1 = do.fem.functionspace(self.grid.mesh, ("DG", 0))
        self.V = do.fem.functionspace(self.grid.mesh, ("Lagrange", 1))


    def create_ds_dx(self) -> None:
        self.ds = ufl.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
        self.dx = ufl.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())


    def split_solution(self) -> None:
        self.P.x.array[:] = self.X.x.array


    def update_P_old(self) -> None:
        self.P_old.x.array[:] = self.P.x.array


    def set_initial_P(self, P_field: to.Tensor) -> None:
        self.P_old.x.array[:] = P_field
        self.P.x.array[:] = P_field


    @abstractmethod
    def solve(self, t: float, dt: float):
        pass


    @abstractmethod
    def set_material(self, material: Material) -> None:
        pass




class MassPorousMedia(MassEquationBase):
    def __init__(self, grid: GridHandlerGMSH, is_linear: bool=True, omega: float=1.0):
        super().__init__(grid, is_linear)
        self.omega = omega


    def set_material(self, material: Material) -> None:
        self.mat = material
        self.M.x.array[:] = self.mat.cf*self.mat.porosity + self.mat.cs*(self.mat.biot - self.mat.porosity)
        density = do.fem.Function(self.DG0_1)
        density.x.array[:] = self.mat.fluid_density
        self.rho_g = density*do.fem.Constant(self.grid.mesh, do.default_scalar_type(tuple(self.mat.gravity)))


    def solve(self, t: float, dt: float) -> None:
        # Update time step
        self.dt.value = dt

        # Update boundary conditions
        self.bc.update_dirichlet(t)
        self.bc.update_neumann(t)
    
        # Update boundary conditions
        self.bc.update_bcs(t)

        tol = 1e-9
        error = 2*tol
        ite = 0
        max_ite = 250

        while error > tol and ite < max_ite:
            # Hold pressure field from previous iteration
            p_k = self.P.x.array.copy()

            # Calculate permeability field
            self.perm.x.array[:] = self.mat.kappa.compute()

            # Build bilinear form
            a = (self.M*self.dp*self.p_ + self.dt*ufl.dot(self.perm*ufl.grad(self.dp), ufl.grad(self.p_)))*self.dx
            bilinear_form = do.fem.form(a)
            A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
            A.assemble()

            # Build linear form
            L = (self.M*self.P_old*self.p_)*self.dx + self.dt*sum(self.bc.neumann_bcs)
            L += self.dt*ufl.dot(self.perm*self.rho_g, ufl.grad(self.p_))*self.dx
            linear_form = do.fem.form(L)
            b = do.fem.petsc.assemble_vector(linear_form)
            do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
            b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            do.fem.petsc.set_bc(b, self.bc.dirichlet_bcs)
            b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            # Solve linear system
            self.solver.setOperators(A)
            self.solver.solve(b, self.X.x.petsc_vec)
            self.X.x.scatter_forward()
            self.split_solution()

            # Apply relaxation
            self.P.x.array[:] = (1 - self.omega)*p_k + self.omega*self.P.x.array

            # Increment iteration counter
            ite += 1

            # Compute error
            error = to.norm(to.from_numpy(self.P.x.array - p_k), p=to.inf).item()

            if self.is_linear:
                error = 0.0
                break

        # Update old temperature field
        self.update_P_old()

        return ite, error


class MassDeformablePorousMedia(MassEquationBase):
    def __init__(self, 
                 grid: GridHandlerGMSH, 
                 is_linear: bool=True, 
                 fixed_stress: bool=True,
                 omega: float=1.0,
                 delta: float=1.0):
        super().__init__(grid, is_linear)

        self.fixed_stress = fixed_stress
        self.omega = omega
        self.delta = delta


    def create_fenicsx_fields(self) -> None:
        self.biot = do.fem.Function(self.DG0_1)
        self.Q_bar = do.fem.Function(self.DG0_1)
        self.perm = do.fem.Function(self.DG0_1)
        self.M = do.fem.Function(self.DG0_1)
        self.G = do.fem.Function(self.DG0_1)
        self.P_old = do.fem.Function(self.V)
        self.P = do.fem.Function(self.V)
        self.X = do.fem.Function(self.V)
        self.h_cell_2 = None


    def set_stabilization_h(self, h : to.Tensor) -> None:
        self.h_cell_2 = do.fem.Function(self.DG0_1)
        self.h_cell_2.x.array[:] = h**2


    def set_material(self, material: Material) -> None:
        self.mat = material
        try:
            biot_to = 1 - self.mat.cs/self.mat.K
        except:
            raise ValueError("Material property 'K' has not been set. Elastic element was probably not defined.")
        self.biot.x.array[:] = biot_to
        self.M.x.array[:] = self.mat.cf*self.mat.porosity + self.mat.cs*(biot_to - self.mat.porosity)
        self.Q_bar.x.array[:] = (self.mat.biot**2)/(self.delta*self.mat.K)
        self.G.x.array[:] = self.mat.G


    def set_u(self, u: do.fem.function.Function) -> None:
        self.u = u


    def update_u_old(self, u: do.fem.function.Function) -> None:
        self.u_old = u.copy()


    def solve(self, t: float, dt: float) -> None:
        # Update time step
        self.dt.value = dt

        # Update boundary conditions
        self.bc.update_dirichlet(t)
        self.bc.update_neumann(t)
    
        # Update boundary conditions
        self.bc.update_bcs(t)

        # Hold pressure from previous fixed-stress iteration
        P_last = self.P.copy()

        # Define tolerances for non-linear loop
        tol = 1e-5
        error = 2*tol
        ite = 0
        max_ite = 50

        while error > tol and ite < max_ite:
            # Hold pressure field from previous iteration
            p_k = self.P.x.array.copy()

            # Calculate permeability field
            self.perm.x.array[:] = self.mat.kappa.compute()

            if self.h_cell_2 is not None:
                stab = self.h_cell_2*self.biot**2/self.G

            # Build bilinear form
            a = (self.M*self.dp*self.p_)*self.dx
            if self.fixed_stress:
                a += (self.Q_bar*self.dp*self.p_)*self.dx
            if self.h_cell_2 is not None:
                a += ufl.dot(stab*ufl.grad(self.dp), ufl.grad(self.p_))*self.dx
            a += self.dt*ufl.dot(self.perm*ufl.grad(self.dp), ufl.grad(self.p_))*self.dx
            bilinear_form = do.fem.form(a)
            A = do.fem.petsc.assemble_matrix(bilinear_form, bcs=self.bc.dirichlet_bcs)
            A.assemble()

            # Build linear form
            L = (self.M*self.P_old*self.p_)*self.dx
            if self.fixed_stress:
                L += (self.Q_bar*P_last*self.p_)*self.dx
            if self.h_cell_2 is not None:
                L += ufl.dot(stab*ufl.grad(self.P_old), ufl.grad(self.p_))*self.dx
            L -= (self.biot*self.p_*(ufl.div(self.u - self.u_old)))*self.dx
            L += self.dt*ufl.dot(self.perm*self.rho_g, ufl.grad(self.p_))*self.dx
            L += self.dt*sum(self.bc.neumann_bcs)
            linear_form = do.fem.form(L)
            b = do.fem.petsc.assemble_vector(linear_form)
            do.fem.petsc.apply_lifting(b, [bilinear_form], [self.bc.dirichlet_bcs])
            b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            do.fem.petsc.set_bc(b, self.bc.dirichlet_bcs)
            b.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)

            # Solve linear system
            self.solver.setOperators(A)
            self.solver.solve(b, self.X.x.petsc_vec)
            self.X.x.scatter_forward()
            self.split_solution()

            self.P.x.array[:] = (1 - self.omega)*p_k + self.omega*self.P.x.array

            # Increment iteration counter
            ite += 1

            # Compute error
            error = to.norm(to.from_numpy(self.P.x.array - p_k), p=to.inf).item()

            if self.is_linear:
                error = 0.0
                break

        return ite, error, P_last


