import os
from os.path import exists
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

# for DL
from tensorflow.keras.models import Sequential
from tensorflow.keras import Model
from tensorflow.keras.layers import Dense, Dropout, Conv1D, add, Flatten, LSTM, AveragePooling1D, concatenate, MaxPool1D
from tensorflow.keras.callbacks import EarlyStopping, Callback, ModelCheckpoint, ReduceLROnPlateau, TensorBoard
from tensorflow.keras.regularizers import l1, l2, L1L2
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K
from tensorflow.keras import Input
from tensorflow.keras.utils import plot_model 
from tensorboard.plugins.hparams import api as hp

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