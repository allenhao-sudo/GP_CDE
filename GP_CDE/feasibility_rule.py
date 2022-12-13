import numpy as np


def feasibility_rule(P, dim, Q=None):
    if Q is None:
        best = P[0, :]
        for id in range(1, P.shape[0]):
            if P[id, dim + 1] < best[dim + 1]:
                best = P[id, :]
            elif P[id, dim + 1] == best[dim + 1] and P[id, dim] < best[dim]:
                best = P[id, :]
            else:
                pass
        return best

    else:
        pop = np.copy(Q)
        for id in range(P.shape[0]):
            if P[id, dim + 1] < Q[id, dim + 1]:
                pop[id, :] = P[id, :]
            elif P[id, dim + 1] == Q[id, dim + 1] and P[id, dim] < Q[id, dim]:
                pop[id, :] = P[id, :]
            else:
                pass
        return pop
