import json
import pandas as pd


def copt_string_to_result(fpath):
    name = (
        fpath.split("/")[-1]
        .replace("s_pre", "pre")
        .replace(".mps", "")
        .replace(".gz", "")
        .replace(".json", "")
    )
    with open(fpath, "r") as f:
        content = json.load(f)
        res_primal = content["PrimalRes"]
        res_dual = content["DualRes"]
        sol_time = float(content["Runtime"])
        status = content["Status"]
        ipm_num = content.get("BarrierIter", 0)
        simplex_num = content.get("SimplexIter", 0)

        if status not in [2, 3]:
            val_primal = float(content["ObjVal"])
            val_dual = float(content["ObjVal"])
            sol_status = "OPTIMAL" if status == 1 else "UNFINISHED"
        elif status == 2:
            sol_status = "INFEASIBLE"
            val_primal = -1e6
            val_dual = -1e6
        elif status == 3:
            sol_status = "UNBOUNDED"
            val_primal = -1e6
            val_dual = -1e6
        else:
            sol_status = "UNFINISHED"
            val_primal = -1e6
            val_dual = -1e6

    return dict(
        res_primal=res_primal.__round__(8),
        res_dual=res_dual.__round__(8),
        sol_time=sol_time.__round__(4),
        val_primal=val_primal.__round__(4),
        val_dual=val_dual.__round__(4),
        sol_status=sol_status,
        ipm_num=ipm_num,
        iteration_num=simplex_num,
        name=name,
    )


if __name__ == "__main__":
    import sys

    result = copt_string_to_result(sys.argv[1])
    print(result)
