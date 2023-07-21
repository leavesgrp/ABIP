"""typical abip schema
{
}
"""
import json
import pandas as pd


def abip_string_to_result(fpath):
  with open(fpath, 'r') as f:
    content = json.load(f)
    res_primal = content['pres']
    res_dual = content['dres']
    sol_time = content['time']
    sol_status = content['status']
    val_primal = content['pobj']
    val_dual = content['dobj']
    name = fpath.split("/")[-1].split(".")[0]
    return dict(iteration_num=content['admm_iter'],
                ipm_num=content['ipm_iter'],
                res_primal=res_primal,
                res_dual=res_dual,
                sol_time=sol_time,
                val_primal=val_primal,
                val_dual=val_dual,
                sol_status=sol_status,
                name=name)
