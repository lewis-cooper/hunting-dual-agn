"""Cache repeated scalar population-prior evaluations without changing their values."""

from functools import lru_cache

import numpy as np


class CachedMagnitudePrior:
    def __init__(self, prior, kind, maxsize=100000,mag_step=.01):
        self.prior=prior
        self.kind=kind
        self.mag_edges=prior.mag_edges
        if kind=='quasar':
            self.evaluate=lru_cache(maxsize=maxsize)(
                lambda redshift,magnitude:float(prior(np.array([redshift]),magnitude)[0]))
            if type(prior).__name__=='GridQSOPrior':
                grid=np.arange(prior.mag_edges[0],prior.mag_edges[-1]+1e-9,mag_step)
                indices={float(magnitude):index for index,magnitude in enumerate(grid)}
                # Use the same magnitude-then-redshift interpolation order as
                # GridQSOPrior; only the repeated magnitude interpolation is batched.
                columns=np.array([np.interp(grid,prior.mag_centres,row) for row in prior.sigma])
                inside=(grid>=prior.mag_edges[0])&(grid<=prior.mag_edges[-1])
                @lru_cache(maxsize=64)
                def curve(redshift):
                    if not prior.z_edges[0]<=redshift<=prior.z_edges[-1]:
                        return np.zeros(len(grid))
                    values=np.array([np.interp(np.array([redshift]),prior.z_centres,column)[0]
                                     for column in columns.T])
                    return np.where(inside,values,0.)
                scalar=self.evaluate
                def evaluate(redshift,magnitude):
                    index=indices.get(magnitude)
                    return float(curve(redshift)[index]) if index is not None else scalar(redshift,magnitude)
                self.evaluate=evaluate
        else:
            self.evaluate=lru_cache(maxsize=maxsize)(
                lambda magnitude,longitude,latitude:float(prior(np.array([magnitude]),longitude,latitude)[0]))

    def __call__(self, first, *arguments):
        values=np.asarray(first)
        if values.size!=1 or any(np.asarray(value).size!=1 for value in arguments):
            return self.prior(first,*arguments)
        scalars=[float(values.reshape(-1)[0]),*[float(np.asarray(value).reshape(-1)[0]) for value in arguments]]
        return np.array([self.evaluate(*scalars)])


def cache_population_priors(model):
    """Only memoize repeated calls made by BlendModel._log_mag_density."""
    model.base.priors={label:(CachedMagnitudePrior(quasar,'quasar',mag_step=model.mag_step),
                              CachedMagnitudePrior(background,'background',mag_step=model.mag_step))
                       for label,(quasar,background) in model.base.priors.items()}
    return model
