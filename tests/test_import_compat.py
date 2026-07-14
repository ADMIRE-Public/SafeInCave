"""
Import compatibility test for SafeInCave refactoring.

This test ensures all current import idioms continue to work during the
package restructuring. It should be run before and after each refactoring
PR to verify zero breakage.
"""
import safeincave as sf
from safeincave import (
    MomentumBC,
    HeatBC,
    CavernBC,
    PostProcessingTools,
    Utils,
)


class TestImportCompatibility:
    """Verify all current import patterns work correctly."""

    def test_top_level_class_imports(self):
        """Test that top-level class imports work (sf.X pattern)."""
        # Core classes
        assert hasattr(sf, 'GridHandlerGMSH')
        assert hasattr(sf, 'HeatDiffusion')
        assert hasattr(sf, 'LinearMomentumBase')
        assert hasattr(sf, 'LinearMomentum')
        assert hasattr(sf, 'LinearMomentumMixed')
        assert hasattr(sf, 'Material')
        assert hasattr(sf, 'SaveFields')
        assert hasattr(sf, 'SimulationLogging')
        assert hasattr(sf, 'CavernThermodynamics')
        assert hasattr(sf, 'ScreenPrinter')
        
        # Simulators
        assert hasattr(sf, 'Simulator_TM')
        assert hasattr(sf, 'Simulator_T')
        assert hasattr(sf, 'Simulator_M')
        
        # Time controllers
        assert hasattr(sf, 'TimeControllerBase')
        assert hasattr(sf, 'TimeController')
        assert hasattr(sf, 'TimeControllerParabolic')
        assert hasattr(sf, 'TimeControllerAdaptive')
        
        # Convergence criteria
        assert hasattr(sf, 'ConvergenceCriterion')
        assert hasattr(sf, 'StrainBasedCriterion')
        assert hasattr(sf, 'ForceResidualCriterion')
        assert hasattr(sf, 'DisplacementIncrementCriterion')
        assert hasattr(sf, 'ForceDisplacementCriterion')
        
        # Utility functions
        assert hasattr(sf, 'register_variable')
        assert hasattr(sf, 'get_variable')
        assert hasattr(sf, 'list_registered_variables')

    def test_module_imports(self):
        """Test that module imports work (sf.Module pattern)."""
        assert hasattr(sf, 'MomentumBC')
        assert hasattr(sf, 'HeatBC')
        assert hasattr(sf, 'CavernBC')
        assert hasattr(sf, 'PostProcessingTools')
        assert hasattr(sf, 'Utils')

    def test_submodule_class_imports(self):
        """Test that submodule class imports work."""
        # MomentumBC classes
        assert hasattr(sf.MomentumBC, 'GeneralBC')
        assert hasattr(sf.MomentumBC, 'DirichletBC')
        assert hasattr(sf.MomentumBC, 'NeumannBC')
        assert hasattr(sf.MomentumBC, 'BcHandler')
        
        # HeatBC classes
        assert hasattr(sf.HeatBC, 'GeneralBC')
        assert hasattr(sf.HeatBC, 'DirichletBC')
        assert hasattr(sf.HeatBC, 'NeumannBC')
        assert hasattr(sf.HeatBC, 'RobinBC')
        assert hasattr(sf.HeatBC, 'BcHandler')
        
        # CavernBC classes
        assert hasattr(sf.CavernBC, 'CavernHandler')
        assert hasattr(sf.CavernBC, 'Cavern_T')
        assert hasattr(sf.CavernBC, 'Cavern_PT')
        assert hasattr(sf.CavernBC, 'Cavern_MassFlux')
        assert hasattr(sf.CavernBC, 'CavernVolumeComputer')
        assert hasattr(sf.CavernBC, 'HeatFluxComputer')

    def test_from_submodule_imports(self):
        """Test 'from safeincave.X import Y' pattern."""
        from safeincave.MomentumBC import GeneralBC, DirichletBC, NeumannBC, BcHandler
        from safeincave.HeatBC import GeneralBC as HeatGeneralBC
        from safeincave.CavernBC import Cavern_PT, CavernHandler
        from safeincave.Utils import create_field_nodes
        
        # Verify they're the right classes
        assert GeneralBC.__name__ == 'GeneralBC'
        assert HeatGeneralBC.__name__ == 'GeneralBC'
        assert Cavern_PT.__name__ == 'Cavern_PT'
        assert callable(create_field_nodes)

    def test_constitutive_models_imported(self):
        """Test that constitutive models are available."""
        # Check that some constitutive models are exported
        assert hasattr(sf, 'Spring')
        assert hasattr(sf, 'Thermoelastic')
        assert hasattr(sf, 'Viscoelastic')
        assert hasattr(sf, 'LinearDashpot')
        assert hasattr(sf, 'DislocationCreep')

    def test_extension_models_available(self):
        """Test that extension models are loaded if available."""
        # These may or may not be present depending on extensions
        # Just check the mechanism works
        if hasattr(sf, 'PlasticRankine'):
            assert callable(sf.PlasticRankine)
        if hasattr(sf, 'PlasticDruckerPrager'):
            assert callable(sf.PlasticDruckerPrager)

    def test_utils_functions_available(self):
        """Test that utility functions are accessible."""
        assert hasattr(sf.Utils, 'create_field_nodes')
        assert hasattr(sf.Utils, 'create_field_elems')
        assert hasattr(sf.Utils, 'numpy2torch')
        assert hasattr(sf.Utils, 'project')
        assert hasattr(sf.Utils, 'epsilon')

    def test_postprocessing_functions_available(self):
        """Test that postprocessing functions are accessible."""
        assert hasattr(sf.PostProcessingTools, 'read_cell_tensor')
        assert hasattr(sf.PostProcessingTools, 'read_cell_scalar')
        assert hasattr(sf.PostProcessingTools, 'read_node_scalar')
        assert hasattr(sf.PostProcessingTools, 'read_node_vector')
        assert hasattr(sf.PostProcessingTools, 'build_smoother')
        assert hasattr(sf.PostProcessingTools, 'build_mapping')

    def test_simulator_mout_exists(self):
        """Test that Simulator_Mout exists (even if not exported at top level)."""
        # Simulator_Mout is in Simulators.py but not exported from __init__.py
        # It should still be accessible via the module
        from safeincave.Simulators import Simulator_Mout
        assert Simulator_Mout.__name__ == 'Simulator_Mout'
