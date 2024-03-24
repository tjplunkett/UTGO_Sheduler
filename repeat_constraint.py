"""
repeat_constraint.py

Author: Tom Plunkett (UTAS)

Purpose: Defines the RepeatConstraints class for use in the manual scheduler.

Date: 24/03/2024
"""
# Import necessary packages 
from astroplan import Constraint, is_observable, min_best_rescale
from astropy import units as u
from astropy.time import TimeDelta
import numpy as np

class RepeatConstraint(Constraint):
    """
    Constraint for making sure repeats behave properly (i.e don't image the same object again if just observed).
    
    params:
    min - The minimum number of hours (dimensionless) between observations of a target - defaults to 0.5 
    max - The maximum numbers of hours (dimensionless) between observations of a target - defaults to 4
    t_end - The end time of the previous observation of a target
    boolean_constraint - We want this constraint to be strict (i.e dont rank, just apply this mask)
    
    returns:
    mask - The times where this constraint is adhered to.
    """
    # Initialise object
    def __init__(self, min=None, max=None, t_end=None, boolean_constraint=True):
        self.min = min if min is not None else 0.5
        self.max = max if max is not None else 4
        self.t_end = t_end
        self.boolean_constraint = boolean_constraint
    
    # Do the actual work of applying this constraint. Check if observing times are available from 30 mins
    # up to 4 hours after the last observation.
    def compute_constraint(self, times, observer, targets):
        mask = np.logical_and(times > self.t_end + TimeDelta(self.min*u.hour), times < self.t_end + TimeDelta(self.max*u.hour))
        return mask
         