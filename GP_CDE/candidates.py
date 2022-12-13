import numpy as np
from feasibility_rule import feasibility_rule
from copy import deepcopy
from utils import  to_unit_cube
import gpytorch
import torch


def current_to_rand(Pt, problem, F, CR, Lambda, index=0, seed=0):
    np.random.seed(seed)
    pop_size = Pt.shape[0]
    pop2 = np.zeros([Lambda, problem.dim + 2])
    for i in range(Lambda):
        v1, v2, v3 = index, index, index
        while v1 == index:
            v1 = np.random.randint(pop_size)
        while (v2 == index) or (v2 == v1):
            v2 = np.random.randint(pop_size)
        while (v3 == index) or (v3 == v2) or (v3 == v1):
            v3 = np.random.randint(pop_size)

        # --- Mutation ---
        v = Pt[index, :problem.dim] + F * (Pt[v2, :problem.dim] - Pt[v3, :problem.dim]) \
            + F * (Pt[v1, :problem.dim] - Pt[index, :problem.dim])
        # --- Cross over ---
        co = np.random.rand(problem.dim)
        for j in range(problem.dim):
            if co[j] <= CR:
                pop2[i, j] = v[j]
            else:
                pop2[i, j] = Pt[index, j]

        # --- Forced crossing ---
        j = np.random.randint(problem.dim)
        pop2[i, j] = v[j]
    pop2[:, :problem.dim] = np.maximum(pop2[:, :problem.dim], problem.lb)
    pop2[:, :problem.dim] = np.minimum(pop2[:, :problem.dim], problem.ub)
    return pop2


def de_best_1(Pt, problem, F, CR, Lambda, index=0, seed=0, best=None):
    np.random.seed(seed)
    pop_size = Pt.shape[0]
    pop2 = np.zeros([Lambda, problem.dim + 2])
    for i in range(Lambda):
        v1, v2 = index, index
        while v1 == index:
            v1 = np.random.randint(pop_size)
        while (v2 == index) or (v2 == v1):
            v2 = np.random.randint(pop_size)
        # while (v3 == index) or (v3 == v2) or (v3 == v1):
        #     v3 = np.random.randint(npop)

        # --- Mutation ---
        v = best[:problem.dim] + F * (Pt[v1, :problem.dim] - Pt[v2, :problem.dim])

        # --- Cross over ---
        co = np.random.rand(problem.dim)
        for j in range(problem.dim):
            if co[j] <= CR:
                pop2[i, j] = v[j]
            else:
                pop2[i, j] = Pt[index, j]

        # --- Forced crossing ---
        j = np.random.randint(problem.dim)
        pop2[i, j] = v[j]
    pop2[:, :problem.dim] = np.maximum(pop2[:, :problem.dim], problem.lb)
    pop2[:, :problem.dim] = np.minimum(pop2[:, :problem.dim], problem.ub)
    return pop2


def candidates(problem, model, mu, sigma, Pt, F, CR, Lambda, seed=0, type="Phase I", best=None):
    np.random.seed(seed)
    npop = Pt.shape[0]
    Qt = np.zeros([npop, problem.dim + 2])

    if type == "Phase I":
        for i in range(npop):
            # generate trial vectors for each vector by DE/current_to_rand/1
            trial_vectors = current_to_rand(Pt, problem, F, CR, Lambda, index=i, seed=seed)
            trial_vectors[:, problem.dim] = problem.obj(trial_vectors[:, :problem.dim])
            trial_vectors[:, problem.dim + 1] = model_predict(deepcopy(trial_vectors[:, :problem.dim])
                                                              , model, mu, sigma, problem)
            trial_vectors[:, problem.dim + 1] = np.maximum(trial_vectors[:, problem.dim + 1], 0)
            Qt[i, :] = feasibility_rule(trial_vectors, problem.dim)

        return Qt

    elif type == "Phase II":
        for i in range(npop):
            # generate trial vectors for each vector by DE/best/1
            trial_vectors = de_best_1(Pt, problem, F, CR, Lambda, index=i, seed=seed, best=best)
            trial_vectors[:, problem.dim] = problem.obj(trial_vectors[:, :problem.dim])
            trial_vectors[:, problem.dim + 1] = model_predict(deepcopy(trial_vectors[:, :problem.dim])
                                                              , model, mu, sigma, problem)
            trial_vectors[:, problem.dim + 1] = np.maximum(trial_vectors[:, problem.dim + 1], 0)
            Qt[i, :] = feasibility_rule(trial_vectors, problem.dim)

        return Qt


def model_predict(X_cand, gp, mu, sigma, problem):
    X_cand = to_unit_cube(X_cand, problem.lb, problem.ub)

    # Figure out what device we are running on
    if len(X_cand) < 1024 * 5:
        device, dtype = torch.device("cpu"), torch.float64
    else:
        device, dtype = torch.device("cuda"), torch.float64

    # We may have to move the GP to a new device
    n_func = len(gp)
    n_cand = X_cand.shape[0]
    y_cand = np.zeros([n_cand, n_func])
    for i in range(n_func):
        gp[i] = gp[i].to(dtype=dtype, device=device)

    with torch.no_grad(), gpytorch.settings.max_cholesky_size(2000):
        X_cand_torch = torch.tensor(X_cand).to(device=device, dtype=dtype)
        for k in range(n_func):
            y_cand[:, k] = gp[k].likelihood(gp[k](X_cand_torch)).mean.cpu().detach().numpy()

    # De-standardize the sampled values
    y_cand = mu + sigma * y_cand
    y_cand = np.maximum(y_cand, 0)
    sum_g = np.sum(y_cand, axis=1)

    return sum_g
