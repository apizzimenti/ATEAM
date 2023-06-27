
import numpy as np


def constant(temperature):
    """
    A constant annealing schedule.

    Args:
        temperature (float): The temperature to be returned.

    Returns:
        A function passed to a model constructor.
    """
    def _(t):
        return temperature
    
    return _


def critical(field):
    """
    A constant annealing schedule which calculates the critical temperature
    of the model.

    Args:
        field (int): The order of the field we're over.

    Returns:
        A function passed to a Model constructor that returns the critical
        temperature of the Potts model.
    """
    p = -np.log(1-(np.sqrt(field)/(1+np.sqrt(field))))

    def _(t):
        return p
    
    return _
