from coptpy import *
from scipy import sparse
import numpy as np
import pandas as pd
import os
import cosmopy as cosmo



#### first call julia is slow
c1_P = sparse.csc_matrix([[4., 1], [1, 2]])
c1_q = np.array([1., 1])
c1_A = sparse.csc_matrix([[1., 1], [1, 0], [0, 1], [-1., -1], [-1, 0], [0, -1]])
c1_b = np.array([1, 0.7, 0.7, -1, 0, 0])
c1_cone = {"l" : 6 }

# create a solver model
c1_model = cosmo.Model()

# setup the model with the problem data and some optional solver settings
c1_model.setup(P = c1_P, q = c1_q, A = c1_A, b = c1_b, cone = c1_cone, verbose = False, eps_abs = 1e-5, max_iter = 4000)

# optional: warm starting of x
c1_x = np.array([1., 0])
c1_model.warm_start(x = c1_x)

# solve the problem
for i in range(100):
    c1_model.optimize()
########################################################







cbf_files = []   

copt_time_list = []
copt_iter_list = []
copt_obj_list = []
copt_status_list = []

cosmo_time_list = []
cosmo_iter_list = []
cosmo_obj_list = []
cosmo_status_list = []


path = 'path to your cblib dataset' 
dirs = os.listdir(path)   
for i in dirs:                            
    if os.path.splitext(i)[1] == ".gz":   
        cbf_files.append(i)


inf = COPT.INFINITY
eps = 1e-4
time_limit = 100.0
iter_limit = int(1e10)

for file_name in cbf_files:

    env = Envr()
    model = env.createModel("test_cbf")

    model.readCbf(path + file_name)

    vars = model.getVars().getAll()

    sense = model.getAttr(COPT.Attr.ObjSense)


    blc = model.getInfo(COPT.Info.LB, model.getConstrs())
    buc = model.getInfo(COPT.Info.UB, model.getConstrs())
    blx = [var.LB for var in vars]
    bux = [var.UB for var in vars]


    ind_blc = []
    ind_buc = []

    for i in range(len(blc)):

        if blc[i] != -inf:
            ind_blc.append(i)
        if buc[i] != inf:
            ind_buc.append(i)


    ind_blx = []
    ind_bux = []

    for i in range(len(blx)):

        if blx[i] != -inf:
            ind_blx.append(i)
        if bux[i] != inf:
            ind_bux.append(i)


    blc = [blc[i] for i in ind_blc]
    buc = [buc[i] for i in ind_buc]
    blx = [blx[i] for i in ind_blx]
    bux = [bux[i] for i in ind_bux]


    len_blc = len(blc)
    len_buc = len(buc)
    len_blx = len(blx)
    len_bux = len(bux)


    A = model.getA()
    c = model.getObjective()
    m,n = A.shape

    obj = [0]*n

    for i in range(c.getSize()):
        obj[c.getVar(i).getIdx()] = c.getCoeff(i)
    c = obj

    cone_builder_array = model.getConeBuilders()
    cone_num = cone_builder_array.getSize()


    soc = []
    cones_sub = []

    for i in range(cone_num):

        cone_builder = cone_builder_array.getBuilder(i)
        cone_size = cone_builder.getSize()
        cone_type = cone_builder.getType()
        cone_vars = cone_builder.getVars().getAll()
        cone_idx = [var.getIdx() for var in cone_vars]
        cones_sub += cone_idx
        soc.append(cone_size)


    row1 = sparse.hstack([sparse.vstack([A[ind_blc,:], sparse.coo_matrix(([1]*len_blx, (list(range(len_blx)), ind_blx)), shape=(len_blx,n))], format='csc'), -sparse.eye(len_blc+len_blx, format='csc'), sparse.csc_matrix((len_blc+len_blx, len_buc+len_bux))])
    row2 = sparse.hstack([sparse.vstack([A[ind_buc,:], sparse.coo_matrix(([1]*len_bux, (list(range(len_bux)), ind_bux)), shape=(len_bux,n))], format='csc'), sparse.csc_matrix((len_buc+len_bux, len_blc+len_blx)), sparse.eye(len_buc+len_bux, format='csc')])

    A_new = sparse.vstack([row1, row2], format='csc')
    b_new = np.array(blc+blx+buc+bux)
    c_new = np.array(c+[0]*(len_blc+len_blx+len_buc+len_bux))

    free_cone = list(set(list(range(n))) - set(cones_sub))

    cones = {"f":len(free_cone), "l": len_blc+len_blx+len_buc+len_bux, "q":soc}

    A_new = sparse.hstack([A_new[:,free_cone], A_new[:,n:], A_new[:,cones_sub]], format='csc').transpose()
    c_new = np.hstack((c_new[free_cone], c_new[n:], c_new[cones_sub]))


    if sense == COPT.MAXIMIZE:
        c_new *= -1

    
    model.setParam(COPT.Param.TimeLimit, time_limit)
    model.setParam(COPT.Param.FeasTol, eps)
    model.setParam(COPT.Param.DualTol, eps)
    model.setParam(COPT.Param.RelGap, eps)
    model.setParam(COPT.Param.AbsGap, eps)
    model.setParam(COPT.Param.BarIterLimit, 50000)


    model.solve()
    copt_time= model.solvingtime
    copt_val = model.objval
    copt_iter= model.BarrierIter
    copt_status =  model.status

    copt_time_list.append(copt_time)
    copt_obj_list.append(copt_val)
    copt_iter_list.append(copt_iter)
    copt_status_list.append(copt_status)



    cosmo_model = cosmo.Model()
    n = len(b_new)
    cosmo_model.setup(P = sparse.csc_matrix((n,n)), q = -b_new, A = A_new.tocsc(), b = c_new, cone = cones, verbose = True, eps_abs = eps, eps_rel = eps, max_iter=iter_limit, time_limit=time_limit)
    cosmo_model.optimize()

    # query solution info
    cosmo_obj_val = -cosmo_model.get_objective_value() # optimal objective vale
    cosmo_status = cosmo_model.get_status() # solution status
    cosmo_iter = cosmo_model.get_iter() # number of iterations until convergence
    cosmo_time = cosmo_model.get_times()["solver_time"] # a dictionary of timings 

    cosmo_time_list.append(cosmo_time)
    cosmo_iter_list.append(cosmo_iter)
    cosmo_obj_list.append(cosmo_obj_val)
    cosmo_status_list.append(cosmo_status)



df = pd.DataFrame({'data':cbf_files,'copt time':copt_time_list, 
                    'copt iter': copt_iter_list, 'copt_obj':copt_obj_list, 
                    'copt status':copt_status_list, 'cosmo time':cosmo_time_list, 
                    'cosmo iter': cosmo_iter_list, 'cosmo_obj':cosmo_obj_list, 
                    'cosmo status':cosmo_status_list})


csv_name = 'cblib_cosmo_copt.csv'
df.to_csv(csv_name)