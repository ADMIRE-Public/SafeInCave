# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from typing import TYPE_CHECKING, Sequence, Union
import torch as to

from .IO import create_field_elems

if TYPE_CHECKING:
    from ..Mesh.Grid import GridHandlerGMSH
    from ..Materials.Material import Material

StressLike = Union[to.Tensor, Sequence[Sequence[float]]]


class ElementSet:
    """A set of element indices bound to the grid they belong to.

    Carries its own grid reference so downstream helpers (e.g.
    `build_uniform_stress`, `build_geostatic_gradient_stress`) only need this
    object, not a separate `grid` argument. Optionally also carries the
    `Material` assigned to those elements, so density can be read from it
    automatically instead of being passed in separately.
    """

    def __init__(
        self, grid: "GridHandlerGMSH", indices: Sequence[int],
        material: "Material" = None,
    ):
        self.grid = grid
        self.indices = list(indices)
        self.material = material

    @classmethod
    def region(
        cls, grid: "GridHandlerGMSH", name: str, material: "Material" = None,
    ) -> "ElementSet":
        """Build an ElementSet from a named mesh region (`grid.region_indices[name]`)."""
        if name not in grid.region_indices:
            raise ValueError(
                f"Unknown region {name!r}. Available regions: "
                f"{sorted(grid.region_indices.keys())}."
            )
        return cls(grid, grid.region_indices[name], material=material)


def build_uniform_stress(element_set: ElementSet, sigma: StressLike) -> to.Tensor:
    """
    Build a per-element initial-stress tensor from a single uniform value.

    Parameters
    ----------
    element_set : ElementSet
        Elements to assign `sigma` to; elements outside it are left at zero.
        To cover several regions with different values, sum multiple calls,
        e.g. `build_uniform_stress(set_a, sig_a) + build_uniform_stress(set_b, sig_b)`.
    sigma : (3, 3) tensor/array/list
        Stress tensor applied to every element in `element_set`.

    Returns
    -------
    torch.Tensor
        Tensor of shape `(element_set.grid.n_elems, 3, 3)` and dtype
        `torch.float64`, ready to pass to `LinearMomentumBase.apply_initial_stress`.
    """
    grid = element_set.grid
    sigma0 = to.zeros((grid.n_elems, 3, 3), dtype=to.float64)
    sigma0[element_set.indices] = to.as_tensor(sigma, dtype=to.float64)
    return sigma0


def build_geostatic_gradient_stress(
    element_set: ElementSet,
    minimum_horizontal_stress_ratio=1.0,
    maximum_horizontal_stress_ratio=1.0,
    density=None,
    gravity: float = 9.81,
    ref_pos: float = None,
    vertical_axis: int = 2,
) -> to.Tensor:
    """
    Build a depth-dependent geostatic initial-stress tensor.

    The vertical stress grows linearly with depth below `ref_pos`,
    `vertical_stress(z) = density * gravity * (ref_pos - z)`, and the two
    horizontal components are `minimum_horizontal_stress_ratio` and
    `maximum_horizontal_stress_ratio` times that vertical stress (no shear).

    Parameters
    ----------
    element_set : ElementSet
        Elements to assign the computed gradient stress to; elements outside
        it are left at zero. Also provides the grid used for `get_parameter`
        and for evaluating element-centroid coordinates.
    density : float, per-region sequence, or per-element tensor, optional
        Passed through `element_set.grid.get_parameter`. If omitted, taken
        from `element_set.material.density` (i.e. the material already
        assigned to these elements); if neither is available, raises.
    gravity : float
        Gravitational acceleration magnitude.
    ref_pos : float, optional
        Coordinate along `vertical_axis` where the vertical stress is zero
        (e.g. ground surface elevation). Defaults to the top of the model,
        i.e. the maximum mesh-node coordinate along `vertical_axis`.
    minimum_horizontal_stress_ratio, maximum_horizontal_stress_ratio :
        float, per-region sequence, or per-element tensor
        Ratio of minimum/maximum horizontal stress to vertical stress.
        Passed through `element_set.grid.get_parameter`. Assigned
        respectively to the lower- and higher-index axis among the two
        horizontal axes.
    vertical_axis : int
        Index (0=x, 1=y, 2=z) of the vertical/depth axis. Defaults to z.

    Returns
    -------
    torch.Tensor
        Tensor of shape `(element_set.grid.n_elems, 3, 3)` and dtype
        `torch.float64`, ready to pass to `LinearMomentumBase.apply_initial_stress`.
    """
    grid = element_set.grid

    if density is None:
        if element_set.material is None:
            raise ValueError(
                "density is missing: pass it explicitly, or build the ElementSet "
                "with material=... (e.g. ElementSet.region(grid, name, material=mat))."
            )
        density = element_set.material.density

    density_t = grid.get_parameter(density)
    k0_min_t = grid.get_parameter(minimum_horizontal_stress_ratio)
    k0_max_t = grid.get_parameter(maximum_horizontal_stress_ratio)

    coord_fns = [lambda x, y, z: x, lambda x, y, z: y, lambda x, y, z: z]
    coord = create_field_elems(grid, coord_fns[vertical_axis])

    if ref_pos is None:
        ref_pos = grid.mesh.geometry.x[:, vertical_axis].max()

    vertical_stress = density_t * gravity * (ref_pos - coord)

    horizontal_axes = [axis for axis in range(3) if axis != vertical_axis]
    diag = to.zeros((grid.n_elems, 3), dtype=to.float64)
    diag[:, vertical_axis] = vertical_stress
    diag[:, horizontal_axes[0]] = k0_min_t * vertical_stress
    diag[:, horizontal_axes[1]] = k0_max_t * vertical_stress

    sigma0 = to.diag_embed(diag)
    mask = to.zeros(grid.n_elems, dtype=to.bool)
    mask[element_set.indices] = True
    sigma0[~mask] = 0.0

    return sigma0
