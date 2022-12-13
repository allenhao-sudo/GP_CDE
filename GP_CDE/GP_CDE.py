import multiprocessing
import os
import time
import numpy as np
from candidates import candidates
from CEC2010 import CEC2010
from CEC2017 import CEC2017
import torch
from gp import train_gp
from utils import to_unit_cube
from copy import deepcopy
import gpytorch
import matlab.engine
from feasibility_rule import feasibility_rule


class SymmetricLatinHypercube(object):
    """Symmetric Latin Hypercube experimental design

    :param dim: Number of dimensions
    :type dim: int
    :param npts: Number of desired sampling points
    :type npts: int

    :ivar dim: Number of dimensions
    :ivar npts: Number of desired sampling points
    """

    def __init__(self, dim, npts):
        self.dim = dim
        self.npts = npts

    def _slhd(self):
        """Generate a matrix with the initial sample points,
        scaled to the unit hypercube

        :return: Symmetric Latin hypercube design in the unit cube of size npts x dim
        :rtype: numpy.array
        """

        # Generate a one-dimensional array based on sample number
        points = np.zeros([self.npts, self.dim])
        points[:, 0] = np.arange(1, self.npts + 1)

        # Get the last index of the row in the top half of the hypercube
        middleind = self.npts // 2

        # special manipulation if odd number of rows
        if self.npts % 2 == 1:
            points[middleind, :] = middleind + 1

        # Generate the top half of the hypercube matrix
        for j in range(1, self.dim):
            for i in range(middleind):
                if np.random.random() < 0.5:
                    points[i, j] = self.npts - i
                else:
                    points[i, j] = i + 1
            np.random.shuffle(points[:middleind, j])

        # Generate the bottom half of the hypercube matrix
        for i in range(middleind, self.npts):
            points[i, :] = self.npts + 1 - points[self.npts - 1 - i, :]

        return points / self.npts

    def generate_points(self):
        """Generate a matrix with the initial sample points,
        scaled to the unit hypercube

        :return: Symmetric Latin hypercube design in the unit cube of size npts x dim
            that is of full rank
        :rtype: numpy.array
        :raises ValueError: Unable to find an SLHD of rank at least dim + 1
        """

        rank_pmat = 0
        pmat = np.ones((self.npts, self.dim + 1))
        xsample = None
        max_tries = 100
        counter = 0
        while rank_pmat != self.dim + 1:
            xsample = self._slhd()
            pmat[:, 1:] = xsample
            rank_pmat = np.linalg.matrix_rank(pmat)
            counter += 1
            if counter == max_tries:
                raise ValueError("Unable to find a SLHD of rank at least dim + 1, is npts too smal?")
        return xsample


def GP_CDE(run_index, pop_size, Lambda, func_id, D, maxeval):
    seed = np.random.randint(1e6)
    np.random.seed(seed)
    dirname = os.path.dirname(__file__)
    eng = matlab.engine.start_matlab()  # start the matlab engine since the benchmarks are implemented in matlab
    eng.cd(dirname+'/CEC2017/', nargout=0)
    problem = CEC2017.CEC2017(int(func_id), D, eng)
    dim, lb, ub, ng = problem.dim, problem.lb, problem.ub, problem.ng  # get basic info of problem
    pop_size = pop_size  # num for initial LHS design
    num_eval = 0
    max_eval = maxeval  # max num of real function eval
    fbounds = (0.25, 0.75)  # the difference amplification factor. Values of 0.5-0.8 are good in most cases.
    cbounds = (0.25, 1)  # The cross-over probability. Use 0.9 to test for fast convergence, and smaller
    # values (~0.1) for a more elaborate search.
    database = np.empty([0, dim + 2])
    his_g = np.empty([0, ng])
    t1 = time.time()
    train_time = 0  # Recording the time cost for GP training
    pop_g = np.zeros([pop_size, ng])

    # initial design with LHS
    Pt = np.zeros([pop_size, dim + 2])
    exp_des = SymmetricLatinHypercube(dim, pop_size)
    Pt[:, :dim] = lb + exp_des.generate_points() * (ub - lb)

    for i in range(pop_size):
        Pt[i, dim] = problem.obj(Pt[i, :dim])  # Objective are treat as cheap for the CEC benchmarks
        # All constraints are treated as expensive for the CEC benchmarks
        pop_g[i, :] = problem.expensive_con(Pt[i, :dim])

        # Add the constraint violation of pop into global dataset
    his_g = np.vstack((his_g, pop_g))

    # compute the total constraint violation
    pop_g = np.maximum(pop_g, 0)
    Pt[:, dim + 1] = np.sum(pop_g, axis=1, keepdims=False)
    num_eval += pop_size
    database = np.vstack((database, Pt))

    # main loop
    while num_eval < max_eval:
        # Warp inputs
        X = to_unit_cube(deepcopy(database[:, :dim]), lb, ub)

        # Standardize values
        fX = deepcopy(his_g)

        # build global gp for expensive constraints
        tic = time.time()
        GP, mu, sigma = model_fit(X, fX, n_training_steps=25)
        train_time = train_time + (time.time() - tic)

        F = np.random.uniform(*fbounds)  # The difference amplification
        CR = np.random.uniform(*cbounds)  # The crossover rate
        fea_index = np.where(database[:, dim + 1] <= 1e-6)[0]

        if np.sum(fea_index) == 0:
            # Phase I using the DE/rand/1 mutation
            type = "Phase I"
            print('**process--{}**: evaluation:{} , no feasible solution'.format(run_index, num_eval))
            Qt = candidates(problem, GP, mu, sigma, Pt, F, CR, Lambda, seed, type)
        else:
            # Phase II using the DE/best/1 mutation
            type = "Phase II"
            best = feasibility_rule(Pt, dim)
            print('**process--{}**: evaluation:{} , best solution:{}'.format(run_index, num_eval, best[dim]))
            Qt = candidates(problem, GP, mu, sigma, Pt, F, CR, Lambda, seed, type, best)

        # Evaluating the pop Qt in expensive functions
        for i in range(pop_size):
            pop_g[i, :] = problem.expensive_con(Qt[i, :dim])

        # add the constraint violation of pop into global dataset
        his_g = np.vstack((his_g, pop_g))

        # compute the total constraint violation
        pop_g = np.maximum(pop_g, 0)
        Qt[:, dim + 1] = np.sum(pop_g, axis=1, keepdims=False)

        # update the training data and number of evaluations
        database = np.vstack((database, Qt))
        num_eval += Qt.shape[0]

        # update the pop based on feasibility rule
        Pt = feasibility_rule(Qt, dim, Pt)

    # save the result
    time_cost = time.time() - t1  # the time cost of the entire process
    fea_index = np.where(database[:, dim + 1] <= 1e-6)[0]
    fea_x = np.empty([0, dim])
    fea_y = np.empty([0, 1])
    best_x = np.empty([0, dim])
    best_y = np.empty([0, 1])
    if np.sum(fea_index) == 0:
        print('**process--{}**: evaluation:{} , no feasible solution'.format(run_index, num_eval))
        print("time_cost:", time_cost)
    else:
        minrows = np.argmin(database[fea_index, dim])
        minrows = fea_index[minrows]
        fea_x = database[fea_index, :dim]
        fea_y = database[fea_index, dim]
        best_x = database[minrows, :dim]
        best_y = database[minrows, dim]
        print('**process--{}**: evaluation:{} , best solution:{}'.format(run_index, num_eval,
                                                                         database[minrows, :].flatten()))
        print("time_cost:", time_cost)
    result = {'database': database.flatten(), 'time': time_cost, 'train_time': train_time,
              'feasible_x': fea_x.flatten(), 'feasible_y': fea_y.flatten(),
              "best_x": best_x.flatten(), 'best_y': best_y.flatten()}

    filename = dirname + '/results/CEC2010/Function_' + str(func_id + 1) + '/dim_' + str(D) + '/Maxeval_' + str(
        maxeval) + '/'
    if not os.path.isdir(filename):
        os.makedirs(filename)
    # np.save(filename + str(run_index) + ".npy", result)


def model_fit(X, fX, n_training_steps):
    # Standardize function values.
    mu, sigma = np.median(fX, axis=0), np.std(fX, axis=0)

    sigma = np.where(sigma < 1e-6, 1, sigma)

    fX = (deepcopy(fX) - mu) / sigma

    # Figure out what device we are running on
    if len(X) < 1024 * 5:
        device, dtype = torch.device("cpu"), torch.float64
    else:
        device, dtype = torch.device("cuda"), torch.float64

    # We use CG + Lanczos for training if we have enough data
    n_func = fX.shape[1]
    gp = [0] * n_func
    with gpytorch.settings.max_cholesky_size(2000):
        X_torch = torch.tensor(X).to(device=device, dtype=dtype)
        for i in range(n_func):
            y_torch = torch.tensor(fX[:, i]).to(device=device, dtype=dtype)
            gp[i] = train_gp(
                train_x=X_torch, train_y=y_torch, use_ard=True, num_steps=n_training_steps, hypers={}
            )
    del X_torch, y_torch
    return gp, mu, sigma


def run_parallel(func_set):
    Lambda = 3000
    pop_size = 25
    dim = 10
    maxeval = 500
    num_run = len(func_set) * 30
    para = np.zeros([num_run, 2], dtype='int32')
    for i in range(len(func_set)):
        for j in range(30):
            para[i * 30 + j, :] = [func_set[i], j]

    for m in range(int(num_run / 5)):
        Processes = []
        for n in range(5):
            p = multiprocessing.Process(target=GP_CDE,
                                        args=(para[m * 5 + n, 1], pop_size, Lambda, para[m * 5 + n, 0], dim, maxeval))
            Processes.append(p)
            p.start()
        for process in Processes:
            process.join()


if __name__ == "__main__":
    CEC2010_set = [0, 6, 7, 12]
    CEC2017_set = [0, 1, 3, 4, 12, 19]
    GP_CDE(1, 25, 3000, CEC2017_set[2], 10, 500)
    # run_parallel(CEC2010_set)
