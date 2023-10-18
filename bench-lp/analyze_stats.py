import numpy as np
import pandas as pd
from scipy import stats
import sys

status = {
    "-": "-",
    "optimal": "optimal",
    "time_limit": "time_limit",
    "failure": "failed",
    "unfinished": "failed",
    "dual_infeasible": "dual_inf",
    "primal_infeasible": "prim_inf",
    "solved": "optimal",
    "solved/inaccurate": "optimal",
    "opt": "optimal",
}

if len(sys.argv) < 4 or sys.argv[3] not in {"sol_time", "matvec"}:
    print(
        """
    please indicate the metric to compute SGM;
    try one in {'sol_time', 'matvec'}
    """
    )
    exit(-1)

# SGM_TARGET = "sol_time"
# SGM_TARGET = "matvec"
SGM_TARGET = sys.argv[3]

pd.set_option("display.max_columns", None)
fname = sys.argv[1]
key = sys.argv[2]

fname_noaffix = fname.split(".")[0]


def bool_success(row):
    k = row.name[0]
    xl = row["sol_status"]
    x = xl.lower()
    if row[SGM_TARGET] == 15000 or np.isnan(row[SGM_TARGET]):
        return 0
    if "solved" in x or "optimal" in x:
        return 1
    if xl == df_std["sol_status"][k]:
        return 1
    else:
        return 0


df = pd.read_excel(fname).fillna({SGM_TARGET: 15000})
df["sol_status"] = df["sol_status"].apply(str.lower).apply(lambda x: status.get(x, x))

df.set_index(["name", "method"], inplace=True)
df_std = df.query(f"method == '{key}'")
df_std = df_std.droplevel(1)


def cal_dev_perc(row, k):
    try:
        return (row[k] - df_std[k].get(row.name[0])) / df_std[k].get(row.name[0])
    except:
        return np.nan


def cal_dev(row, k):
    try:
        return row[k] - df_std[k].get(row.name[0])
    except:
        return np.nan


def cal_sol_time(row):
    k = row.name[0]
    if row["success"]:
        return row[SGM_TARGET]
    return 15000


def scaled_gmean(arr):
    return stats.gmean(arr.fillna(15000) + 10) - 10


df[["iteration_num", SGM_TARGET]] = df[["iteration_num", SGM_TARGET]].apply(
    pd.to_numeric, errors="coerce"
)
df["success"] = df.apply(bool_success, axis=1)
df[SGM_TARGET] = df.apply(cal_sol_time, axis=1)
# df['more_iter'] = pd.to_numeric(df.apply(lambda row: cal_dev(row, 'iteration_num'), axis=1), errors='coerce')
df["more_time"] = pd.to_numeric(
    df.apply(lambda row: cal_dev(row, SGM_TARGET), axis=1), errors="coerce"
)
df["std_time"] = df.apply(lambda row: df_std[SGM_TARGET].get(row.name[0]), axis=1)
# df['std_iter'] = df.apply(lambda row: df_std['iteration_num'].get(row.name[0]), axis=1)
# df['more_iter_perc'] = pd.to_numeric(df['more_iter'].divide(df['std_iter']), errors='coerce')
# df['more_time_perc'] = pd.to_numeric(df['more_time'].divide(df['std_time']), errors='coerce')
dfa = df.groupby(level=1).aggregate(
    {
        # "more_iter": 'mean',
        # "more_iter_perc": 'mean',
        # "iteration_num": 'mean',
        # "more_time": 'mean',
        # "more_time_perc": 'mean',
        SGM_TARGET: scaled_gmean,
        "success": "sum",
    }
)
df_sol = dfa[SGM_TARGET]
df_sol_m = df_sol / df_sol[key]
dfc = pd.DataFrame.from_dict(
    {
        "mean": dfa[SGM_TARGET].to_dict(),
        "scale": df_sol_m.to_dict(),
        "solved": dfa["success"].to_dict(),
    },
    orient="index",
)

df_final = pd.concat([dfc, df[SGM_TARGET].unstack(level=1)], sort=False)
print(df_final.to_latex(longtable=True, multirow=True, multicolumn=True))
print(dfa.to_latex(longtable=True, multirow=True, multicolumn=True))
df_final.to_excel(f"{fname_noaffix}.analysis.xlsx")
dfa.to_excel(f"{fname_noaffix}.analysis-sum.xlsx")
