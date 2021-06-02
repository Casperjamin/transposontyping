import pandas as pd

def frame_merge(list_of_dfs, stretchdel, output_csv):
    frames_to_merge = []
    for i in list_of_dfs[0].split(' '):
        if i.split('_')[-3] == stretchdel:
            frames_to_merge.append(i)
    df = pd.concat([pd.read_csv(x, sep = "\t", index_col = 0) for x in frames_to_merge], sort=False)
    df.to_csv(output_csv, sep = "\t")
    return output_csv
