from astroplan import Constraint, is_observable, min_best_rescale

class TimeConstraint(Constraint):
    """
    Constraint for making sure repeats behave properly (i.e don't image the same object again if just observed)
    """
    def __init__(self, min=None, max=None, boolean_constraint=True):
        """
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable time between repeats. `None`
            indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable time between repeats. `None`
            indicates no limit.
        """
        self.min = min if min is not None else 0*u.hours
        self.max = max if max is not None else 2*u.hours
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times, observer, targets, filters):
        if self.boolean_constraint:
            if len(targets) > 1:
                for i in range(0, len(targets)):
                    if targets[i] = targets[i+1] and filters[i] = filters[i+1:
                        mask = np.logical_and(times > times + self.min, times < times + max_time))
            return mask

        # Need to figure this out...
        else:
            # rescale the values so that they become
            # scores between zero and one
            rescale = min_best_rescale(times, times+self.min,
                                       times+self.max, less_than_min=0)
            return rescale
 