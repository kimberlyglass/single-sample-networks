"""
Central configuration for the scalability/runtime benchmark.

Edit the values below to change the experimental grid. All other scripts
import from this file so there is a single source of truth for paths,
parameter grids, and repetition counts.
"""

from pathlib import Path

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
CODE_DIR = Path(__file__).resolve().parent
SIM_DIR = CODE_DIR.parent
DATA_DIR = SIM_DIR / "data"
RESULTS_DIR = SIM_DIR / "results"

MANIFEST_PATH = DATA_DIR / "manifest.csv"
RAW_TIMINGS_PATH = RESULTS_DIR / "raw_timings.csv"
SUMMARY_PATH = RESULTS_DIR / "summary_timings.csv"

# -----------------------------------------------------------------------
# Reproducibility
# -----------------------------------------------------------------------
BASE_SEED = 42

# -----------------------------------------------------------------------
# Number of repetitions per (n, p) setting
# -----------------------------------------------------------------------
N_REPS = 3

# -----------------------------------------------------------------------
# Methods to benchmark. method_ID follows run_SSN_method.py's convention:
#   1 = LIONESS, 2 = SSN, 3 = SWEET, 4 = BONOBO
# Methods are run individually (not method_ID=5/"ALL") so that each
# method's runtime and peak memory are measured in isolation.
# -----------------------------------------------------------------------
METHODS = {
    "LIONESS": 1,
    "SSN": 2,
    "SWEET": 3,
    "BONOBO": 4,
}

# -----------------------------------------------------------------------
# Experimental grid
# -----------------------------------------------------------------------

# Sweep 1: vary number of samples (n), fixed number of genes (p)
FIXED_P_FOR_N_SWEEP = 300
N_SAMPLES_GRID = [20, 50, 100, 200, 400, 800]

# Sweep 2: vary number of genes (p), fixed number of samples (n)
FIXED_N_FOR_P_SWEEP = 100
N_GENES_GRID = [50, 100, 200, 400, 800, 1200]

# -----------------------------------------------------------------------
# Synthetic data generation parameters
# -----------------------------------------------------------------------

N_BLOCKS = 10          # number of co-expression modules
WITHIN_BLOCK_RHO = 0.9  # target within-block correlation; matches the ~0.91
                        # implied correlation in the paper's MakeInSilicoData.m
                        # (unit-variance signal + N(0, 0.1) noise -> rho = 1/(1+0.1))
