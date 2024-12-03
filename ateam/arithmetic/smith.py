
import numpy as np


def SNF(M, z=2):
    """
    Computes the Smith normal(/canonical) form of the matrix M.
    """
    m, n = M.shape

    S = np.zeros((2, max(m,n)), dtype=int)
    for i, L in enumerate([m, n]):
        for k in range(L):
            S[i][k] = k

    for x in range(min(m, n)):
        nonzero = np.transpose((M == 1).nonzero())

        # Column and row swaps. Keep track of these swaps, so we can identify
        # which cubes/faces end up in the bases of the groups.
        for k, l in nonzero:
            if x <= k and x <= l:
                M[[x,k]] = M[[k,x]]
                M[:,[x,l]] = M[:,[l,x]]

                sl = S[0][x]
                sr = S[0][k]
                S[0][x] = sr
                S[0][k] = sl

                sl = S[1][x]
                sr = S[1][l]
                S[1][x] = sr
                S[1][l] = sl
                break

        # Adding rows and columns.
        if len(nonzero):
            for i in range(x+1, m):
                if M[i,x]:
                    M[i] = (M[i] + M[x]) % z

            for j in range(x+1, n):
                if M[x,j]:
                    M[:,j] = (M[:,j] + M[:,x]) % z
    
    return M, S
