import numpy as np
import matlab


class CEC2010:
    def __init__(self, fun_id, D, engine):
        Xmin1 = 0 * np.ones([1, D])
        Xmax1 = +10 * np.ones([1, D])
        Xmin2 = -5.12 * np.ones([1, D])
        Xmax2 = +5.12 * np.ones([1, D])
        Xmin3 = -1000 * np.ones([1, D])
        Xmax3 = +1000 * np.ones([1, D])
        Xmin4 = -50 * np.ones([1, D])
        Xmax4 = +50 * np.ones([1, D])
        Xmin5 = -600 * np.ones([1, D])
        Xmax5 = +600 * np.ones([1, D])
        Xmin6 = -600 * np.ones([1, D])
        Xmax6 = +600 * np.ones([1, D])
        Xmin7 = -140 * np.ones([1, D])
        Xmax7 = +140 * np.ones([1, D])
        Xmin8 = -140 * np.ones([1, D])
        Xmax8 = +140 * np.ones([1, D])
        Xmin9 = -500 * np.ones([1, D])
        Xmax9 = +500 * np.ones([1, D])
        Xmin10 = -500 * np.ones([1, D])
        Xmax10 = +500 * np.ones([1, D])
        Xmin11 = -100 * np.ones([1, D])
        Xmax11 = +100 * np.ones([1, D])
        Xmin12 = -1000 * np.ones([1, D])
        Xmax12 = +1000 * np.ones([1, D])
        Xmin13 = -500 * np.ones([1, D])
        Xmax13 = +500 * np.ones([1, D])
        Xmin14 = -1000 * np.ones([1, D])
        Xmax14 = +1000 * np.ones([1, D])
        Xmin15 = -1000 * np.ones([1, D])
        Xmax15 = +1000 * np.ones([1, D])
        Xmin16 = -10 * np.ones([1, D])
        Xmax16 = +10 * np.ones([1, D])
        Xmin17 = -10 * np.ones([1, D])
        Xmax17 = +10 * np.ones([1, D])
        Xmin18 = -50 * np.ones([1, D])
        Xmax18 = +50 * np.ones([1, D])
        lb = [Xmin1, Xmin2, Xmin3, Xmin4, Xmin5, Xmin6, Xmin7, Xmin8, Xmin9,
              Xmin10, Xmin11, Xmin12, Xmin13, Xmin14, Xmin15, Xmin16, Xmin17, Xmin18]
        ub = [Xmax1, Xmax2, Xmax3, Xmax4, Xmax5, Xmax6, Xmax7, Xmax8, Xmax9,
              Xmax10, Xmax11, Xmax12, Xmax13, Xmax14, Xmax15, Xmax16, Xmax17, Xmax18]
        hn = [0, 1, 1, 4, 2, 2, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1]
        gn = [2, 2, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 3, 3, 3, 2, 2, 1]
        self.dim = D
        self.lb = lb[fun_id].ravel()
        self.ub = ub[fun_id].ravel()
        self.bounds = np.hstack((lb[fun_id].T, ub[fun_id].T))
        self.ng = gn[fun_id]
        self.nh = hn[fun_id]
        self.func_id = fun_id+1
        self.engine = engine

    def obj(self, pop):
        pop = matlab.double(pop.tolist())
        [val, _, _] = self.engine.pyTEC(pop, self.func_id, nargout=3)
        val = np.array(val)
        return val.ravel()

    def expensive_con(self, pop):
        pop = matlab.double(pop.tolist())
        [_, g, _] = self.engine.pyTEC(pop, self.func_id, nargout=3)
        g = np.array(g)
        return g
