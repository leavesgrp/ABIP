import json
import pandas as pd


def google_pdhg_string_to_result(fpath):
  with open(fpath, 'r') as f:
    content = json.load(f)
    stats_solution = content['solutionStats']
    res_primal = stats_solution['convergenceInformation'][0]["lInfPrimalResidual"]
    res_dual = stats_solution['convergenceInformation'][0]["lInfDualResidual"]
    sol_time = content['solveTimeSec']
    sol_status = "optimal" if "optimal" in content['terminationReason'].lower() else "unfinished"
    val_primal = stats_solution['convergenceInformation'][0]["primalObjective"]
    val_dual = stats_solution['convergenceInformation'][0]["dualObjective"]

    return dict(iteration_num=content['iterationCount'],
                res_primal=res_primal.__round__(8),
                res_dual=res_dual.__round__(8),
                sol_time=sol_time.__round__(4),
                val_primal=val_primal.__round__(4),
                val_dual=val_dual.__round__(4),
                sol_status=sol_status,
                name=fpath.split("/")[-1].split(".")[0])


if __name__ == '__main__':
  import sys
  result = google_pdhg_string_to_result(
    sys.argv[1])
  print(result)
