"""
Protocol definitions for the Drug Response Prediction Subnet.

This module defines the communication protocols between validators and miners
for drug response prediction tasks using synthetic genomic data.
"""

from typing import Dict, Any, Optional

import bittensor as bt

from niome_subnet.genomics.model import Task

class GenomicsTaskSynapse(bt.Synapse):
    """Protocol for genomics simulation tasks."""

    # Input fields
    task: Optional[Task] = None
    encryption_key: str = ""
    timeout: Optional[float] = None  # Timeout window for submission

    # Output fields (encrypted via hybrid RSA+AES-GCM when encryption_key is set)
    encrypted_vcf: Optional[str] = None         # JSON payload from encryption.encrypt()
    encrypted_annotations: Optional[str] = None # JSON payload from encryption.encrypt()
    elapsed_time: Optional[float] = None
    signature: Optional[str] = None  # Cryptographic signature
    error: Optional[str] = None  # Set when variant calling fails

    def deserialize(self) -> bt.Synapse:
        """Deserialize the GenomicsTaskSynapse Object."""
        return self