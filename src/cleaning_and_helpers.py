import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def clean_rbcl_data(data_filepath,
                    save_filepath):
    df = pd.read_csv(data_filepath, low_memory=False)
    df['GenSp'] = df['Genus'] + '.' + df['Species']
    rbcl_cols = [col for col in df if col.startswith('RBCL')]
    my_cols = ['UniqueID', 'Family', 'Genus', 'GenSp', 'Site', 'Lat', 'Long'] + rbcl_cols
    df = df.reindex(columns=my_cols)
    df = df.dropna(subset=[i for i in rbcl_cols], how='all')
    df = df.fillna(0)
    df.to_csv(save_filepath)
    return df