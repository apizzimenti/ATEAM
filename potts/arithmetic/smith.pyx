
import numpy as np

def SNF(M, z=2, S=None):
    """
    Computes the Smith normal(/canonical) form of the matrix M.
    """
    nonzero = np.transpose((M == 1).nonzero())
    m, n = M.shape
    x = 0

    if m == 0 or n == 0: return M, S

    if S is None: S = np.zeros(M.shape, dtype=int)

    # Column and row swaps.
    for k, l in nonzero:
        if x <= k and x <= l:
            M[[x,k]] = M[[k,x]]
            M[:,[x,l]] = M[:,[l,x]]
            S[x,k] = 1
            S[l,x] = 1
            break

    # Adding rows and columns.
    if len(nonzero):
        for i in range(x+1, m):
            if M[i,x]:
                M[i] = (M[i] + M[x]) % z

        for j in range(x+1, n):
            if M[x,j]:
                M[:,j] = (M[:,j] + M[:,x]) % z

    T, Q = SNF(M[1:m,1:n], z=z, S=S[1:m,1:n])
    M[1:m,1:n] = T
    S[1:m,1:n] = Q
    
    return M, S
