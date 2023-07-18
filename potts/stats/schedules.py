
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
    critical = np.log(1-(np.sqrt(field)/(1+np.sqrt(field))))

    def _(t):
        return critical
    
    return _


def randomizedToConstant(steps, field, hold=1/2, distribution=np.random.normal):
    """
    A temperature schedule which samples spin values (according to `distribution`)
    *less than* the critical temperature of the model, then fixes the temperature
    at the critical temperature for some desired proportion of the runtime.

    Args:
        steps (int): Total number of steps the experiment takes.
        field (int): Order of the field we're over; used to compute critical
            temperature.
        hold (float): Value in [0,1] which represents the proportion of the time
            the temperature is fixed at the critical temperature. Defaults to
            half random, half fixed.
        distribution (Callable): Distribution from which we sample.

    Returns:
        A function that consumes a step number and returns a temperature.
    """
    holdAtStep = int(steps * (1-hold))
    critical = -np.log(1-np.sqrt(field)/(1+np.sqrt(field)))

    def _(t):
        if t < holdAtStep: return distribution()-1
        else: return critical

    return _


def linear(steps, low=-10, high=10):
    """
    A linear temperature schedule.

    Args:
        steps (int): Total number of steps in the experiment.
        low (float): Lowest temperature we assign.
        high (float): Highest temperature we assign.

    Returns:
        A function that consumes a step number and returns a temperature.
    """
    def _(t):
        q = t * ((high-low)/steps) + low
        return q if q != 0 else 0.0000000000001

    return _
