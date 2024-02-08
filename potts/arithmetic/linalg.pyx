
import numpy as np

def sampleFromKernel(A, field, includes=None):
    """
    Uniformly randomly samples a cochain given the coboundary matrix A.
    """
    # If we're including all the columns, just set our submatrix equal to the
    # whole coboundary; otherwise, exclude some of the columns.
    if includes: B = A.take(includes, axis=0)
    else: B = A

    print(B.shape)
    
    # Find a basis for the kernel of B, then take a uniformly random linear
    # combination of the basis elements; this is equivalent to uniformly
    # randomly sampling a cocycle.
    K = B.null_space()
    M, _ = K.shape
    Q = field.Random(M, dtype=np.int32)
    return (Q@K).astype(np.int32)
