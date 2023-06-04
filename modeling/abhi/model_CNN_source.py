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

# Model functions ---------------------------------------------------------------------------------------------------------------------
def model_tuner_two_streams(hp):
    # todo
    # optimize for activation function of dense layers

    # CNN model genetics--------------------------------------------------------------------------------------     
    markers_no = 45443
    CNN_input_g = Input(shape=(markers_no, 1), name = "CNN_in_g")
    
    ### First layer
    CNN_g = Conv1D(filters = hp.Int("CNN_g_f_fl", min_value=64, max_value=512, step = 64, default = 512), 
                 kernel_size = hp.Int("CNN_g_ks_fl", min_value=3, max_value=36, step = 3), 
                 padding='valid', activation='relu', name="CNN_g_fl")(CNN_input_g) # filter_size and  kernel size
    CNN_g = AveragePooling1D(pool_size =  hp.Int("CNN_g_ap_fl", min_value=2, max_value=32, default=16, step = 4), 
                           strides = 3, padding='same', name="CNN_g_ap_fl")(CNN_g) # pool size and strides
    ### Variable layers
    for i in range(hp.Int("CNN_g_num_vl", min_value = 2, max_value = 4)):
        CNN_g = Conv1D(filters = hp.Int(f'CNN_g_f_vl_{i}', min_value=64, max_value=512, step = 32, default = 256), 
                     kernel_size = hp.Int(f'CNN_g_ks_vl_{i}', min_value=3, max_value=36, step = 3), 
                     padding='valid', activation='relu', name=f'CNN_g_Conv_{i}')(CNN_g)
        CNN_g = AveragePooling1D(pool_size = hp.Int(f'CNN_g_ap_vl_{i}', min_value=2, max_value=32, default=16, step = 4), 
                               strides = 3, padding='same', name=f'CNN_g_ap_{i}')(CNN_g) # pool size and strides
    ### Flattening layer
    CNN_g_output = Flatten(name="CNN_g_flatten")(CNN_g)
    
    # CNN model ec--------------------------------------------------------------------------------------     
    ec_no = 765
    CNN_input_ec = Input(shape=(ec_no, 1), name = "CNN_in_ec")
    
    ### Fixed layers
    CNN_ec = Conv1D(filters = 64, kernel_size = 4, padding='valid', activation='relu', name=f'CNN_ce_Conv_0')(CNN_input_ec)
    CNN_ec = AveragePooling1D(pool_size = 4, strides = 1, padding='same', name=f'CNN_ec_ap_0')(CNN_ec)
    
    CNN_ec = Conv1D(filters = 32, kernel_size = 4, padding='valid', activation='relu', name=f'CNN_ce_Conv_1')(CNN_ec)
    CNN_ec = AveragePooling1D(pool_size = 4, strides = 1, padding='same', name=f'CNN_ec_ap_1')(CNN_ec)
    
    CNN_ec = Conv1D(filters = 24, kernel_size = 2, padding='valid', activation='relu', name=f'CNN_ce_Conv_2')(CNN_ec)
    CNN_ec = AveragePooling1D(pool_size = 2, strides = 1, padding='same', name=f'CNN_ec_ap_2')(CNN_ec)
    
    #for i in range(hp.Int("CNN_ec_num_vl", min_value = 1, max_value = 4)):
    #    CNN_ec = Conv1D(filters = hp.Int(f'CNN_ec_f_vl_{i}', min_value=8, max_value=32, step = 2, default = 16), 
    #                 kernel_size = hp.Int(f'CNN_ec_ks_vl_{i}', min_value=1, max_value=4, step = 1), 
    #                 padding='valid', activation='relu', name=f'CNN_ce_Conv_{i}')(CNN_ec)
    #    CNN_ec = AveragePooling1D(pool_size = hp.Int(f'CNN_ec_ap_vl_{i}', min_value=1, max_value=4, default=2, step = 1), 
    #                           strides = 2, padding='same', name=f'CNN_ec_ap_{i}')(CNN_ec) # pool size and strides
    
    ### Flattening layer
    CNN_ec_output = Flatten(name="CNN_ec_flatten")(CNN_ec)
    
    # concatenate models --------------------------------------------------------------------------------------
    model_concat = concatenate([CNN_g_output, CNN_ec_output], name="concat_in")
    
    ### Variable layers
    for i in range(hp.Int("concat_num_vl", min_value = 2, max_value = 8)):
        model_concat = Dense(units = hp.Int(f'concat_unit_vl_{i}', min_value=64, max_value=512, step = 32, default = 256), 
                             activation='relu', name=f'concat_vl_{i}')(model_concat)
        model_concat = Dropout(rate = hp.Float(f'concat_drop_rate_vl_{i}', min_value=0.1, max_value=0.5, step = 0.01), 
                               name=f'concat_drop_vl_{i}')(model_concat)
    ### Final layer
    model_concat_out = Dense(1, activation='tanh',name="concat_out")(model_concat) # 1 unit since its a regression model
    
    # compile model --------------------------------------------------------------------------------------------
    ## hyperparameters
    l_rate = hp.Float("compile_l_rate", min_value=1e-5, max_value=1e-2,
                             sampling="log")
    beta_val_1 = hp.Float("compile_beta_val_1", min_value=0, max_value=1)
    beta_val_2 = hp.Float("compile_beta_val_2", min_value=0, max_value=1)
    
    ## model arch
    compiled_model = Model(inputs=[CNN_input_g, CNN_input_ec], outputs = model_concat_out)
    compiled_model.compile(loss = 'mean_absolute_error',
                           optimizer = Adam(learning_rate = l_rate, #lr 
                                            beta_1 = beta_val_1, #beta1
                                            beta_2 = beta_val_2), # beta2
                           metrics = ['mean_squared_error'])
    
    # output model
    return(compiled_model)