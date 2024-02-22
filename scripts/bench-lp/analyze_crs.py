# test script to summarize all results
# used for collect crossover results.
# @author Chuwen Zhang
#
import json
import logging
import os
import time

import pandas as pd

DEFAULT_CONF = "./conf.analyze.json"

FORMAT = "%(asctime)s %(levelname)s %(name)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

logger = logging.getLogger("@analyze")
logger.setLevel(logging.DEBUG)


def copt_crossover_to_dict(fpath):
    """
    {
      "name": "pre_25fv47",
      "Status": 1,
      "ObjVal": 5501.845888286746,
      "IterCount": 2802,
      "SolvingTime": 0.1771249771118164
    }
    """
    with open(fpath, "r") as f:
        content = json.load(f)
        iternum = content["IterCount"]
        soltime = content["SolvingTime"]
        status = content["Status"]
        objval = content["ObjVal"]
        name = content["name"]
        return dict(
            name=name, iteration=iternum, time=soltime, status=status, objval=objval,
        )


def analyze(fpath=DEFAULT_CONF):
    config = json.load(open(fpath, "r"))
    instance_path = config["directory"]
    dfs = []
    logger.info(f"registered methods {config['methods']}")
    for m in config["methods"]:
        name = m["name"]
        affix = m["affix"]
        solution_path = os.path.join(instance_path, m["solution_dir"])
        results = []
        logger.debug(f"try {m}")
        if not os.path.exists(solution_path):
            logger.debug(f"method {m} does not exist")
            continue
        logger.debug(f"analyze {m} @ {solution_path}")
        for _fp in os.listdir(solution_path):
            fp = os.path.join(solution_path, _fp)
            if not fp.endswith(affix):
                continue
            try:
                r = copt_crossover_to_dict(fp)
                results.append(r)
            except Exception as e:
                logging.exception(e)
                logging.error(f"failed to analyze {name} @ {fp}")
        if len(results) > 0:
            df = (
                pd.DataFrame.from_records(results)
                .assign(
                    method=name,
                    name=lambda df: df["name"].apply(
                        lambda x: x.replace(".mps", "").replace(".gz", "")
                    ),
                )
                .drop_duplicates(subset=["name", "method"])
            )
            dfs.append(df)

        logger.info(f"{name} solution finished")
    df_agg = pd.concat(dfs).fillna(0)

    instances = df_agg["name"].unique()
    method_names = [n["name"] for n in config["methods"]]
    index = pd.MultiIndex.from_tuples(
        list((i, m) for m in method_names for i in instances.tolist()),
        names=("name", "method"),
    )
    df_agg = (
        df_agg.set_index(["name", "method"])
        .reindex(index, fill_value="-")
        .reset_index(drop=False)
        .sort_values(["name", "method"])
        .assign(real_name=lambda df: df["name"].apply(lambda x: x.replace("pre_", "")))
    )

    return df_agg


def scaled_gmean(arr, scale=1):
    """compute scaled geometric mean,
    the scale=10 mimics the Hans's benchmark
    """
    return stats.gmean(arr + scale) - scale


if __name__ == "__main__":
    import argparse
    from scipy import stats

    pd.set_option("display.max_columns", None)
    parser = argparse.ArgumentParser()
    parser.add_argument("--conf", type=str)
    # only output the instances with prefix
    parser.add_argument(
        "--prefix", type=str, default="pre", help="only analyze results with prefix"
    )

    args = parser.parse_args()
    df_agg = analyze(args.conf)
    result_stamp = int(time.time())
    df_agg.to_excel(f"result_{result_stamp}.xlsx", index=False)
    idx_df_agg_with_prefix = df_agg["name"].apply(lambda x: x.startswith(args.prefix))
    df_agg_pre = df_agg[idx_df_agg_with_prefix]
    # is successful? todo, may have to add a ground truth
    bool_fail = lambda x: x != 1
    df_agg_pre = df_agg_pre.assign(bool_fail=lambda df: df.status.apply(bool_fail))
    df_agg_pre.to_excel(f"result_{result_stamp}.{args.prefix}.xlsx", index=False)
    # do some summary stats
    logger.info(f"analyze finished to with stamp {result_stamp}")

    df_agg_suc = df_agg_pre[~df_agg_pre.bool_fail]
    df_agg_suc = df_agg_suc[(df_agg_suc != "-")].dropna()
    # df_sum = df_agg_suc.groupby(
    #     'method'
    # )['name'].count().rename('solved').reset_index()
    df_sum = (
        df_agg_suc.replace(
            {"time": {"-": 1e3}, "iteration": {"-": 1e4}, "status": {"-": 0}}
        )
        .groupby("method")
        .agg({"status": sum, "iteration": scaled_gmean, "time": scaled_gmean})
    )
    df_sum.to_excel(f"result_{result_stamp}.summary.xlsx", index=True)
    print(
        r"""
\documentclass{article}
\usepackage{lscape,longtable,multirow}
\usepackage{booktabs,caption}
\begin{document}
    """
        + f"""
{df_sum.to_latex(multirow=True,longtable=True, multicolumn=True, caption='Number of solved instances')}
    """
        + r"""
\end{document}
    """
    )
