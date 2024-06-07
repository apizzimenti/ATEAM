
import numpy as np

def SNF(M):
    """
    Computes the Smith normal(/canonical) form of the matrix M.
    """
    nonzero = np.transpose((M == 1).nonzero())
    m,n = M.shape
    x = 0

    if m == 0 or n == 0: return M

    # Column and row swaps.
    for k, l in nonzero:
        if x <= k and x <= l:
            M[[x,k]] = M[[k,x]]
            M[:,[x,l]] = M[:,[l,x]]
            break

    # Adding rows and columns.
    if len(nonzero):
        for i in range(x+1, m):
            if M[i,x]:
                M[i] = M[i] + M[x]

        for j in range(x+1, n):
            if M[x,j]:
                M[:,j] = M[:,j] + M[:,x]
    
    M[1:m,1:n] = SNF(M[1:m,1:n])
    return M
