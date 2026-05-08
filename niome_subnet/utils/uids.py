import random
from typing import List

import bittensor as bt
import numpy as np


def check_uid_availability(
    metagraph: "bt.metagraph.Metagraph", uid: int, vpermit_tao_limit: int
) -> bool:
    """Check if uid is available. The UID should be available if it is serving and has less than vpermit_tao_limit stake
    Args:
        metagraph (:obj: bt.metagraph.Metagraph): Metagraph object
        uid (int): uid to be checked
        vpermit_tao_limit (int): Validator permit tao limit
    Returns:
        bool: True if uid is available, False otherwise
    """
    # Filter non serving axons.
    if not metagraph.axons[uid].is_serving:
        return False
    # Filter validator permit > 1024 stake.
    if metagraph.validator_permit[uid]:
        if metagraph.S[uid] > vpermit_tao_limit:
            return False
    # Available otherwise.
    return True


def get_miner_uids(self) -> np.ndarray:
    """
    Filter out uids that are validators in the metagraph.
    """
    uids = []
    for uid in range(self.metagraph.n):
        # skip validators
        if self.metagraph.validator_trust[uid] > 0:
            continue

        # get current block (try property, fallback to attribute)
        try:
            current_block = int(self.block)
        except Exception:
            current_block = int(getattr(self, "current_block", 0))

        # epoch length from config or attribute
        epoch_length = getattr(self.config.neuron, "epoch_length", None)
        if epoch_length is None:
            epoch_length = getattr(self, "epoch_length", 0)

#        if (current_block - int(self.metagraph.last_update[uid])) <= epoch_length:
#            continue

        uids.append(uid)
    
    uids = np.array(uids)
    return uids


def get_random_uids(
    self, k: int, available_uids: List[int] = None
) -> np.ndarray:
    """Returns k available random uids from the metagraph.
    Args:
        k (int): Number of uids to return.
        exclude (List[int]): List of uids to exclude from the random sampling.
    Returns:
        uids (np.ndarray): Randomly sampled available uids.
    Notes:
        If `k` is larger than the number of available `uids`, set `k` to the number of available `uids`.
    """
    # Avoid truth-testing numpy arrays (which raises an error when they contain
    # multiple elements). Treat None or empty sequence as missing.
    if available_uids is None or (hasattr(available_uids, "__len__") and len(available_uids) == 0):
        available_uids = get_miner_uids(self)

    # If k is larger than the number of available uids, set k to the number of available uids.
    # Normalize available_uids to a Python list to make random.sample robust
    if isinstance(available_uids, np.ndarray):
        available_list = available_uids.tolist()
    else:
        try:
            available_list = list(available_uids)
        except Exception:
            # Fallback: wrap scalar into list
            available_list = [available_uids]

    # If k is larger than the number of available uids, set k to the number of available uids.
    k = min(k, len(available_list))

    # If nothing to sample, return empty array.
    if k == 0 or len(available_list) == 0:
        return np.array([], dtype=int)

    # Sample and return as numpy array
    uids = np.array(random.sample(available_list, k))
    return uids
