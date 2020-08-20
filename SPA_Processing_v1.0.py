# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 14:08:53 2020

@author: -
"""


import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
import csv
from scipy.optimize import brentq
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import xlrd





def select_fldr():
    root = tk.Tk()
    root.withdraw()
    fld_selected = filedialog.askdirectory()
    
    if fld_selected == '':
        quit()
    
    return fld_selected



def select_files():
    root = tk.Tk()
    root.withdraw()
    files_selected = filedialog.askopenfilenames()
    
    if files_selected == '':
        quit()
    
    return files_selected


def select_file():
    root = tk.Tk()
    root.withdraw()
    files_selected = filedialog.askopenfilename()
    
    if files_selected == '':
        quit()
    
    return files_selected


def get_start_t_num(file_path):
    ''' Return the row index for data start (at scan header), time stamp, number of data points
    
    Returns list:
        int: [0] 
        
    '''
    
    
    test_idx = []
    data_here = False
    test_start = 0
    
    with open(file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj, delimiter='\t')
        for row in csv_reader:
                          
            if data_here and csv_reader.line_num == test_start + 1:
                time_stamp = float(row[0].split()[1]) # seconds
            elif data_here and csv_reader.line_num == test_start + 2:
                temp = float(row[0].split()[1])
            elif data_here and csv_reader.line_num == test_start + 3:
                suns = float(row[0].split()[1])
            elif data_here and csv_reader.line_num == test_start + 4:
                rh = float(row[0].split()[1])
            elif data_here and csv_reader.line_num > test_start + 7 and len(row) != 3:
                test_end = csv_reader.line_num
                n_pts = test_end - test_start - 7
                test_idx.append([test_start, n_pts, time_stamp, temp, suns, rh])
                data_here = False
            else:
                for row_str in row:
                    if 'START TEST HEADER' in row_str:
                        test_start = csv_reader.line_num
                        # This returns a base-1 line_num
                        data_here = True
                        
                        break
                    
            
        
        if data_here:   # File ends before switch condition encountered (last data set)
            test_end = csv_reader.line_num + 1
            n_pts = test_end - test_start - 7
            test_idx.append([test_start, n_pts, time_stamp, temp, suns, rh])
    
    return test_idx


def get_sample_prg(file_name):
    ''' Extract sample info from file name
    ex: .../1245/1245_pxB/M1245_pXB_IV_Dark16.txt
    
    Returns list:
        str: [0] SampleID
        str: [1] pixel
        str: [2] program ID
    '''
    sample = file_name.split('/')[-1].split("_")[0] # 
    pixel = file_name.split('/')[-1].split("_")[1][-1]
    prg = 'p' + file_name.split('/')[-1].split("_")[3][4:-4]
    
    return [sample, pixel, prg]


def get_ivt(file_name, row_idx, n_pts):
    ''' Read an individual test scan from data file
    
    Returns list:
        np.array:[0] Current
        np.array:[1] Voltage
        np.array:[2] time (seconds from test start)
    '''
    scan_df = pd.read_csv(
        file_name,
        sep = '\t',
        skiprows = row_idx+6,
        names = ['current', 'voltage', 'time'],
        nrows = n_pts
        )
    
    current = scan_df[['current']].values.flatten()
    voltage = scan_df[['voltage']].values.flatten()
    time = scan_df[['time']].values.flatten()
    
    return [current, voltage, time]


def get_ivt_all(file_path):
    ''' Compiles iv scans from SPA file into dataframe
    Function expects a full path
    
    Returns pd.DataFrame columns=[t, voltage, i@t0, i@t1, ...]
    '''
    
    scan_idx = get_start_t_num(file_path)
    for i, testinfo in enumerate(scan_idx):
        current, voltage, time_srs = get_ivt(file_path, testinfo[0], testinfo[1])
        
        if i==0:
            d = {'time':time_srs,'voltage':voltage}
            spa_df = pd.DataFrame(data=d)
            
        spa_df[testinfo[1]] = current
            
    return spa_df


def get_tir(file_name, row_idx):
    ''' Read an individual test scan header from data file
    
    Returns list:
        np.array:[0] temp
        np.array:[1] irradiance (suns)
        np.array:[2] RH (V)
    '''
    
    with open(file_name) as csvfile:
        readCSV = list(csv.reader(csvfile, delimiter='\t'))
        temp = float(readCSV[row_idx + 1][0].split(' ')[1])
        irr = float(readCSV[row_idx + 2][0].split(' ')[1])
        rh = float(readCSV[row_idx + 3][0].split(' ')[1])

    
    return [temp, irr, rh]



def find_pce(voltage: np.array, current, suns, pin: bool):
    """voltage [V] and current [ma/cm2] must be arrays"""
    power = current * voltage
    if pin:
        pce = abs(power.min())
    else:
        pce = power.max()
    
    return pce



def find_jsc(voltage, current):
    """voltage [V] and current [ma/cm2] must be arrays"""
    
    # Check that voltage array is ascending
    if not voltage[1]>voltage[0]:
        voltage = np.flip(voltage,0)
        current = np.flip(current,0)
    jsc = np.interp(0, voltage, current)
    
    return jsc



def find_voc(voltage, current):
    """voltage [V] and current [ma/cm2] must be arrays with voltage ascending"""
    

    
    # Check that voltage array is ascending
    if not voltage[1]>voltage[0]:
        voltage = np.flip(voltage,0)
        current = np.flip(current,0)
    
    # Create an interpolation function and use it to find root (i.e. Voc)        
    f_interp = interpolate.interp1d(voltage, current)
    
    if (f_interp(voltage.min())>=0) and (f_interp(voltage.max())>= 0): # Catch dead pixel
        voc = np.nan
        return voc
    
    voc = brentq(f_interp, voltage.min(), voltage.max())
        
    return voc



def find_ff(voltage, current, suns, pin):
    """voltage [V] and current [ma/cm2] must be arrays"""
    pce = find_pce(voltage, current, suns, pin)
    jsc = find_jsc(voltage, current)
    voc = find_voc(voltage, current)
    
    ff = pce / (jsc * abs(voc))
    
    return ff


def calc_fom(spa_data, fom_params):    
    '''
    

    Parameters
    ----------
    spa_data : pd.DataFrame columns=[t, voltage, i@t0, i@t1, ...]
    fom_params : list [suns, active area]

    Returns
    -------
    fom_df : [t, PCE, Voc, Jsc, FF]

    '''
    suns, active_area = fom_params
    timestamps = np.array(list(spa_data)[2:])
    elapsed_time = (timestamps-timestamps[0])/3660 #convert sec to hr
    
    fom_df = pd.DataFrame(columns=['t','PCE','Voc','Jsc','FF'])
    fom_df['t'] = elapsed_time
    v = spa_data['voltage']
    
    for i, thisscan in enumerate(timestamps):
        iv = spa_data[thisscan]*1000/active_area
        fom_df.at[i,'PCE'] = find_pce(v,iv,suns,True)
        fom_df.at[i,'Voc'] = find_voc(v,iv)
        fom_df.at[i,'Jsc'] = find_jsc(v,iv)
        fom_df.at[i,'FF'] = find_ff(v,iv,suns,True)
    
    return fom_df


def read_date(date):
    ''' Convert timestamp (seconds, 1904) to datetime
    '''
    return xlrd.xldate.xldate_as_datetime(date,1)




def make_tir_df(scanidx):
    '''
    Creates dataframe with time, temp, suns, RH from scan idx

    Parameters
    ----------
    scanidx : list [start row, n pts, timestamp (s, 1904), temp, suns, rh]

    Returns
    -------
    tir_df : [datetime, temp, suns, rh]

    '''
    tir_df = pd.DataFrame([c[2] for c in scanidx], columns=['t'])
    tir_df['t'] = tir_df['t']/86400 #convert from seconds to days
    tir_df['t']=pd.to_datetime(tir_df['t'].apply(read_date), errors='coerce')
    
    tir_df['Temp'] = [c[3] for c in scanidx]
    tir_df['Suns'] = [c[4] for c in scanidx]
    tir_df['RH'] = [c[5] for c in scanidx]

#     timestamps = [c[1] for c in scanidx]
#     dt = xlrd.xldate.xldate_as_datetime(timestamps/86400,0)
#     tir_df = pd.DataFrame(columns=['t','Temp','Irradiance','RH'])
#     for i, thisscan in enumerate(timestamps):
#         tir_df.at[i,'Temp'] = get_tir(scanfile,scanidx[i])
    
    return tir_df



def norm_fom_df(df):
    
    imax = df['PCE'].astype(float).argmax()
    ndf = df.copy()
    ndf['PCE'] = ndf['PCE'] / ndf['PCE'][imax]
    ndf['Voc'] = ndf['Voc'] / ndf['Voc'][imax]
    ndf['Jsc'] = ndf['Jsc'] / ndf['Jsc'][imax]
    ndf['FF'] = ndf['FF'] / ndf['FF'][imax]
    
    return ndf


def plot_fom(df):
    
    
    #TODO: Generate a new plot each time
    #TODO: Add Title
    
    plt.plot(df['t'], df['PCE'], 'o-', label='PCE')
    plt.plot(df['t'], df['Voc'], 'o-', label='Voc')
    plt.plot(df['t'], df['Jsc'], 'o-', label='Jsc')
    plt.plot(df['t'], df['FF'], 'o-', label='FF')
    
    plt.xlabel('Time (hr)')
    plt.ylabel('Normalized FOM')
    plt.legend()
    
    plt.show()
    
    pass
    
    
def get_fom_table():
    
    thisfile = select_file()
    this_df = pd.read_csv(thisfile)
    
    
    return this_df    
    

'''
Obejective: Output stability data from the SPA
    Source:   IVt Datasets
    
    Output:   FOM v. time plot
              FOM v. time table (raw)
              IV database - could be a pointer (sample, time) > (path, idx)
    Format:   FOM v. time plot
              FOM v. time table
              
Steps
1. Get sample info
2. Extract time series
3. Extract IV
4. Calc FOM at t_i
5. Generate raw FOM v. t table
6. Normalize FOM
7. Plot
8. Save
'''


def spa():
    active_area = 0.1
    suns = 1.5
    
    save_fld = select_fldr()
    
    # Get Files
    iv_files = select_files()
    
    # Iterate through files to pull parameters
    for i, thisfile in enumerate(iv_files):
        
        #sample_info = get_sample_prg(thisfile)
        thisname = thisfile.split('/')[-1].split('.')[0]
        this_data = get_ivt_all(thisfile)
        
        fom_df = calc_fom(this_data,[suns, active_area])
        norm_fom = norm_fom_df(fom_df)
        
        plot_fom(norm_fom)
        
        this_data.to_csv(save_fld + '/' + thisname + '.csv', index=False)
        fom_df.to_csv(save_fld + '/' + thisname + '_fom.csv', index=False)
        
        
        
    pass
    

