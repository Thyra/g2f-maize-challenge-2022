import os
from os.path import exists
from pathlib import Path
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID" 
import time

# for data wrangling
import pandas as pd
import numpy as np
import patsy as pt
from re import search

# for plots
import matplotlib.pyplot as plt

# for CV
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import random

# for tuning
import keras_tuner as kt

# for file flux
import pickle
import json

# for logging and debugging
import logging
import warnings

# Common functions ---------------------------------------------------------------------------------------------------------------------
def scale_data(to_transform, pd_cols, pd_index):
    scaler = MinMaxScaler((0,1))
    data_scaled = scaler.fit_transform(to_transform)
    data_scaled_df = pd.DataFrame(data_scaled, columns = pd_cols, index = pd_index)
    return data_scaled_df, scaler

def read_pkl(path):
    with open(path, "rb") as fp:   # Unpickling
        data = pickle.load(fp)
    return data

def write_pkl(data, path):
    with open(path, "wb") as fp:   # pickling
        pickle.dump(data, fp)
    return print("Done")

def read_json(path):
    with open(path, encoding = "utf8") as json_file:
        data = json.load(json_file)
    return data
def write_json(data, path):
    with open(path, "w") as fp:   
        json.dump(data, fp)
    return print("Done")

def print_function(func_name):
    lines = inspect.getsource(func_name)
    return print(lines)

def purge_excess_missing(data, thr = 0.2, plot = False, id_cols = []):
    data_missing = data[data.columns.difference(id_cols)].isna().sum()/data.shape[0]
    if plot:
        data_missing.reset_index(name="n").plot.barh(x='index', y='n',figsize=(20, 10)) # the columns with > 20 % values will be deleted. there are non numeric values in this data.
    to_include = data_missing[data_missing <= thr].index.values.tolist()
    to_include = id_cols + to_include
    data_subset = data.loc[:, to_include]
    return data_subset

def set_dirs(base_dir_path, verbose = True, run_id = None):
    if run_id is None:
        run_id = time.strftime("run_%Y_%m_%d")
        base_folder = base_dir_path + '/' + run_id
    else:
        base_folder = base_dir_path + '/' + f'{str(run_id)}'
    cb_at = base_folder + '/callback_data'
    tb_cb = cb_at + '/tb_cb'
    mc_cb = cb_at + '/mc_cb/'
    pred_at = base_folder + '/pred'
    model_at = base_folder + '/model'
    tmp_at = base_folder + '/tmp_data'
    if(not os.path.isdir(base_folder)):
        os.system(f'mkdir -p  {base_folder} {pred_at} {model_at} {cb_at} {tb_cb} {mc_cb} {tmp_at}')
    if (verbose):
        print(f'base folder at {base_folder}, \ncallbacks at {cb_at}, \npredictions at {pred_at}, \nmodel at {model_at}, \ntmp at {tmp_at}')
    # output
    out = {}
    out['base_folder'] = base_folder
    out['tb_cb'] = tb_cb
    out['mc_cb'] = mc_cb
    out['pred_at'] = pred_at
    out['model_at'] = model_at
    out['tmp_at'] = tmp_at
    return out 
