import yaml
import pandas as pd
import glob
import io


failed_samples = []

def selection(kma_res):
    sample_id = kma_res.split('/')[1]
    df = pd.read_csv(kma_res, sep = '\t')
    if int(df['Template_Identity']) <= 50:
        failed_samples.append(sample_id)

def filter(input_yaml, output_yaml):
    with open(input_yaml, 'r') as open_input:
        loaded_yaml = yaml.safe_load(open_input)
        for i in failed_samples:
            print(f'ignoring sample: {i}')
            del loaded_yaml['SAMPLES'][i]
    with io.open(output_yaml, 'w', encoding='utf8') as outfile:
        yaml.dump(loaded_yaml, outfile, default_flow_style=False, allow_unicode=True)
