
def constant(temperature):
    """
    A constant annealing schedule.

    Args:
        temperature (float): The temperature to be returned.

    Returns:
        A function passed to a model constructor.
    """
    def _(t):
        return -temperature
    
    return _
