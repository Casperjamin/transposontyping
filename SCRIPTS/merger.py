import pandas as pd


def frame_merge(list_of_dfs):
    """take location of all dataframes, returns merged dataframe"""
    df = pd.concat([pd.read_csv(x, sep = "\t", index_col = 0) for x in list_of_dfs])
    return df
