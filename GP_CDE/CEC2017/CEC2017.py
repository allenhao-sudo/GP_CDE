import numpy as np
import matlab


class CEC2017:
    def __init__(self, fun_id, D, engine):
        Xmin1 = -100 * np.ones([1, D])
        Xmax1 = +100 * np.ones([1, D])
        Xmin2 = -100 * np.ones([1, D])
        Xmax2 = +100 * np.ones([1, D])
        Xmin3 = -100 * np.ones([1, D])
        Xmax3 = +100 * np.ones([1, D])
        Xmin4 = -10 * np.ones([1, D])
        Xmax4 = +10 * np.ones([1, D])
        Xmin5 = -10 * np.ones([1, D])
        Xmax5 = +10 * np.ones([1, D])
        Xmin6 = -20 * np.ones([1, D])
        Xmax6 = +20 * np.ones([1, D])
        Xmin7 = -50 * np.ones([1, D])
        Xmax7 = +50 * np.ones([1, D])
        Xmin8 = -100 * np.ones([1, D])
        Xmax8 = +100 * np.ones([1, D])
        Xmin9 = -10 * np.ones([1, D])
        Xmax9 = +10 * np.ones([1, D])
        Xmin10 = -100 * np.ones([1, D])
        Xmax10 = +100 * np.ones([1, D])
        Xmin11 = -100 * np.ones([1, D])
        Xmax11 = +100 * np.ones([1, D])
        Xmin12 = -100 * np.ones([1, D])
        Xmax12 = +100 * np.ones([1, D])
        Xmin13 = -100 * np.ones([1, D])
        Xmax13 = +100 * np.ones([1, D])
        Xmin14 = -100 * np.ones([1, D])
        Xmax14 = +100 * np.ones([1, D])
        Xmin15 = -100 * np.ones([1, D])
        Xmax15 = +100 * np.ones([1, D])
        Xmin16 = -100 * np.ones([1, D])
        Xmax16 = +100 * np.ones([1, D])
        Xmin17 = -100 * np.ones([1, D])
        Xmax17 = +100 * np.ones([1, D])
        Xmin18 = -100 * np.ones([1, D])
        Xmax18 = +100 * np.ones([1, D])
        Xmin19 = -50 * np.ones([1, D])
        Xmax19 = +50 * np.ones([1, D])
        Xmin20 = -100 * np.ones([1, D])
        Xmax20 = +100 * np.ones([1, D])
        Xmin21 = -100 * np.ones([1, D])
        Xmax21 = +100 * np.ones([1, D])
        Xmin22 = -100 * np.ones([1, D])
        Xmax22 = +100 * np.ones([1, D])
        Xmin23 = -100 * np.ones([1, D])
        Xmax23 = +100 * np.ones([1, D])
        Xmin24 = -100 * np.ones([1, D])
        Xmax24 = +100 * np.ones([1, D])
        Xmin25 = -100 * np.ones([1, D])
        Xmax25 = +100 * np.ones([1, D])
        Xmin26 = -100 * np.ones([1, D])
        Xmax26 = +100 * np.ones([1, D])
        Xmin27 = -100 * np.ones([1, D])
        Xmax27 = +100 * np.ones([1, D])
        Xmin28 = -50 * np.ones([1, D])
        Xmax28 = +50 * np.ones([1, D])
        lb = [Xmin1, Xmin2, Xmin3, Xmin4, Xmin5, Xmin6, Xmin7, Xmin8, Xmin9,
              Xmin10, Xmin11, Xmin12, Xmin13, Xmin14, Xmin15, Xmin16, Xmin17, Xmin18,
              Xmin19, Xmin20, Xmin21, Xmin22, Xmin23, Xmin24, Xmin25, Xmin26, Xmin27, Xmin28]
        ub = [Xmax1, Xmax2, Xmax3, Xmax4, Xmax5, Xmax6, Xmax7, Xmax8, Xmax9,
              Xmax10, Xmax11, Xmax12, Xmax13, Xmax14, Xmax15, Xmax16, Xmax17, Xmax18,
              Xmax19, Xmax20, Xmax21, Xmax22, Xmax23, Xmax24, Xmax25, Xmax26, Xmax27, Xmax28]
        hn = [0, 0, 1, 0, 0, 6, 2, 2, 2, 2, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0]
        gn = [1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 1, 2, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 1, 1, 1, 1, 2, 2]
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
        val = self.engine.pyCEC2017(pop, self.func_id)
        val = np.array(val).ravel()
        return val

    def expensive_con(self, pop):
        pop = matlab.double(pop.tolist())
        [_, g, _] = self.engine.pyCEC2017(pop, self.func_id, nargout=3)
        g = np.array(g).reshape(-1, self.ng)
        g = np.maximum(g, 0)
        sum_g = np.sum(g, axis=1, keepdims=False)
        return sum_g