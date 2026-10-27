"""
spatial-sanctuary: Biophysical Reaction-Diffusion Transport Solver &
Evolutionary Sanctuary Dynamics in Compact EGFR-Mutant Non-Small Cell Lung Cancer.

Author: Luca Celebrano (Department of Clinical Medicine and Surgery, School of Medicine, University of Naples Federico II)
ORCID: 0009-0007-0931-3555
License: GNU General Public License v3.0 (GPL-3.0)
"""

from .biophysical_diffusion_model import (
    BiophysicalTransportModel,
    get_standard_modalities,
)
from .replicator_dynamics_game import (
    SpatialGameAssay,
)
from .clinical_validation import (
    MSKCCClinicalValidator,
    run_freireich_benchmark_validation,
)
from .master_pipeline import (
    MasterSpatialOncologyPipeline,
)

__version__ = "2.2.0"
__author__ = "Luca Celebrano"
__email__ = "l.celebrano@studenti.unina.it"

__all__ = [
    "BiophysicalTransportModel",
    "get_standard_modalities",
    "SpatialGameAssay",
    "MSKCCClinicalValidator",
    "run_freireich_benchmark_validation",
    "MasterSpatialOncologyPipeline",
    "__version__",
]
