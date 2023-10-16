import json
import pandas as pd


def google_pdhg_string_to_result(fpath):
  with open(fpath, 'r') as f:
    content = json.load(f)
    stats_solution = content['solution_stats']
    stats_conv = stats_solution['convergence_information'][0]
    res_primal = stats_conv['relative_l2_primal_residual']
    res_dual = stats_conv['relative_l2_dual_residual']
    sol_time = content['solve_time_sec']
    sol_status = content['termination_string']
    val_primal = stats_conv['primal_objective']
    val_dual = stats_conv['primal_objective']
    name = content['instance_name']
    return dict(iteration_num=stats_solution['iteration_number'],
                res_primal=res_primal.__round__(8),
                res_dual=res_dual.__round__(8),
                sol_time=sol_time.__round__(4),
                val_primal=val_primal.__round__(4),
                val_dual=val_dual.__round__(4),
                sol_status=sol_status,
                name=name)


if __name__ == '__main__':
  import sys
  result = google_pdhg_string_to_result(
    "../miplib2017-bench/pdhg_sol/pre_wachplan_summary.json")
  print(result)
