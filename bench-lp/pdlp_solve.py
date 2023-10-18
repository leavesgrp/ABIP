##
# run pdlp and crossover script by W. Gao.
##
import cvxpy as cp
import numpy as np
import argparse
import gurobipy as gp
import json
from coptpy import *

import ortools
from ortools.pdlp.solvers_pb2 import PrimalDualHybridGradientParams, TerminationCriteria
from ortools.pdlp import solve_log_pb2
from ortools.pdlp import solvers_pb2
from ortools.pdlp.python import pywrap_pdlp
from ortools.init import pywrapinit


def copt_solve(fname, timelimit=3600):
    env = Envr()
    model = env.createModel()
    model.read(fname)  # type: Model

    model.solveLP()
    model.setParam(COPT.Param.LpMethod, 0)
    model.setParam(COPT.Param.TimeLimit, 3)

    return


def copt_crossover(fname, A, b, x, y, s, tm=1000):
    env = Envr()
    model = env.createModel()
    model.read(fname)  # type: Model

    model.setLpSolution(x, b - A.dot(x), y, s)
    model.setParam(COPT.Param.LpMethod, 3)
    model.setParam(COPT.Param.TimeLimit, tm)

    model.solveLP()

    return model


def get_lp_data(fname):
    model = gp.read(fname)
    A = model.getA()
    b = np.asarray(model.rhs)
    c = np.asarray(model.obj)
    lb = np.asarray(model.lb)
    mm = pywrap_pdlp.read_quadratic_program_or_die(fname)

    return [A, b, c, lb], mm


def solve_pdhg_lp(A, b, c, lb, mm, tol=1e-08, tm=3600):
    """
    Solve equality constrained LP

    min  c' * x
    s.t.  A * x = b
         l <= x <= u

    Return (x, y, s, w)
    """

    # x = cp.Variable(len(c))
    # constr = [A @ x == b, x >= lb]
    # prob = cp.Problem(cp.Minimize(c.T @ x), constr)
    prob = mm

    criteria = TerminationCriteria(
        eps_optimal_absolute=tol,
        eps_optimal_relative=tol,
        eps_primal_infeasible=tol,
        eps_dual_infeasible=tol,
        time_sec_limit=tm,
    )

    p = solvers_pb2.PrimalDualHybridGradientParams(termination_criteria=criteria)
    p.verbosity_level = 2
    p.presolve_options.use_glop = False

    print(p)
    # prob.solve(solver='PDLP', verbose=True, parameters_proto=p)
    result = pywrap_pdlp.primal_dual_hybrid_gradient(prob, p.SerializeToString())
    solver_log = solve_log_pb2.SolveLog.FromString(result.solve_log_str)
    print("Primal solution:", result.primal_solution)
    print("Dual solution:", result.dual_solution)
    print("Reduced costs:", result.reduced_costs)

    return (
        prob,
        p,
        result,
        solver_log,
        [result.primal_solution, result.dual_solution, result.reduced_costs],
    )


parser = argparse.ArgumentParser()
parser.add_argument("--tol", type=float, default=1e-08)
parser.add_argument("--solve_tm", type=float, default=3600)
parser.add_argument("--cross_tm", type=float, default=-1)
parser.add_argument("--file", type=str, default="")
parser.add_argument("--output", type=str, default=".")


if __name__ == "__main__":
    args = parser.parse_args()
    fname = args.file
    tol = args.tol

    info, mm = get_lp_data(fname)
    prob, p, result, solver_log, sol = solve_pdhg_lp(
        # info[0], info[1], info[2], info[3], tol, tm=args.solve_tm
        info[0],
        info[1],
        info[2],
        info[3],
        mm,
        tol,
        tm=args.solve_tm,
    )
    print("Iterations:", solver_log.iteration_count)
    print("Solve time (sec):", solver_log.solve_time_sec)

    name = fname.split("/")[-1].replace(".mps", "").replace(".gz", "")
    with open(f"{args.output}/{name}.json", "w") as f:
        from google.protobuf.json_format import MessageToJson

        f.write(MessageToJson(solver_log))

    if args.cross_tm > 0:
        model = copt_crossover(
            fname, info[0], info[1], sol[0], sol[1], sol[2], tm=args.cross_tm
        )
        content = {}
        content["name"] = name
        content["Status"] = model.LpStatus
        content["ObjVal"] = model.ObjVal
        content["IterCount"] = model.SimplexIter
        content["SolvingTime"] = model.SolvingTime

        print(content)
        json.dump(content, open(f"{args.output}/{name}.crs.json", "w"))
