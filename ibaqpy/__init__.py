"""
ibaqpy - A Python package for iBAQ (intensity-based absolute quantification) analysis.

This package provides tools for processing and analyzing proteomics data using
the iBAQ method, which allows for absolute quantification of proteins.
"""

import warnings

# Suppress numpy matrix deprecation warning
warnings.filterwarnings("ignore", category=PendingDeprecationWarning,
                       module="numpy.matrixlib.defmatrix")

__version__ = "0.0.5"

# Import logging configuration
from ibaqpy.ibaq.logging_config import initialize_logging

# Initialize logging with default settings
# Users can override these settings by calling initialize_logging with their own settings
initialize_logging()