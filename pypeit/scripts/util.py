"""
Script utilities
"""
import sys
import numpy as np

def specified_args(args, options=None):
    """
    Determine which of the arguments defined in an ArgumentParser
    instance were defined by the user. 

    Args:
        args (`argparse.ArgumentParser`):
            Object used to parse the command-line arguments.
    """
    spec_keys = list(args._option_string_actions.keys())
    specified = []
    if options is None:
        options = sys.argv[1:]
    for key in options:
        if key in spec_keys:
            specified += [args._option_string_actions[key].dest]
    return np.unique(specified)
