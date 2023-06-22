
def always():
    """
    An acceptance function which always accepts the proposed state; we do nothing
    with the arguments to this function.

    Args:
        chain (Chain): Chain object containing the model and current state.

    Returns:
        A function which always returns True.
    """
    def _(chain): return True
    return _
