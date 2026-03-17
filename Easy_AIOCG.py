# -*- coding: utf-8 -*-
"""
@author: Nicolas Piette-Lauziere
Created on 2021/01/11
Updated on 2026/03/17

Script to plot AIOCG diagram with bar codes from Montreuil et al.(2013).
This example uses data from Park et al.
Can. J. Earth Sci. 51: 1–24 (2014) dx.doi.org/10.1139/cjes-2013-0120
IMPORTANT: Adjust path to the folder containing data and the AIOCG.png background at lines 245 and 258.
"""
import os
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%% Functions
def weigth2molar(df):
    """
    Calculate whole rock molar % from oxide wt%. 
    Imput must be a dataframe structured like: df.Al2O3, df.FeO, df.K2O etc.
    with a column named df.Specimen to identify the data
    
    Output is a new dataframe structured like: df.Al, df. Fe, df.K
    If both Fe2O3 and FeO are present, output is: df.FeII, df.FeIII
    Additional columns are ignored.
    """
    wSi = 28.085
    wTi = 47.867 
    wAl = 26.982
    wFe = 55.845
    wMn = 54.938
    wMg = 24.305
    wCa = 40.078
    wNa = 22.990
    wK = 39.098
    wP = 30.974
    wO = 15.999
    #Create a new dataframe called DF
    DF = df.filter(items=['Specimen'])

    DF['Si'] = (df.SiO2*wSi/(wSi+2*wO))/wSi
    DF['Ti'] = (df.TiO2*wTi/(wTi+2*wO))/wTi
    DF['Al'] = (df.Al2O3*2*wAl/(2*wAl+3*wO))/wAl
    if 'FeO' in df.columns:
        DF['Fe'] = (df.FeO*wFe/(wFe+wO))/wFe
    elif 'Fe2O3' in df.columns:
        DF['Fe'] = (df.Fe2O3*2*wFe/(2*wFe+3*wO))/wFe
    elif 'Fe2O3'in df.columns and 'FeO' in df.columns:
        DF['FeIII'] = (df.Fe2O3*2*wFe/(2*wFe+3*wO))/wFe
        DF['FeII'] = (df.FeO*wFe/(wFe+wO))/wFe
    DF['Mn'] = (df.MnO*wMn/(wMn+wO))/wMn
    DF['Mg'] = (df.MgO*wMg/(wMg+wO))/wMg
    DF['Ca'] = (df.CaO*wCa/(wCa+wO))/wCa
    DF['Na'] = (df.Na2O*2*wNa/(2*wNa+wO))/wNa
    DF['K'] = (df.K2O*2*wK/(2*wK+wO))/wK
    DF['P'] = (df.P2O5*2*wP/(2*wP+5*wO))/wP
    
    sum = DF.sum(axis=1)
    
    DF['Si'] = DF.Si/sum*100
    DF['Ti'] = DF.Ti/sum*100
    DF['Al'] = DF.Al/sum*100
    DF['Fe'] = DF.Fe/sum*100
    DF['Mn'] = DF.Mn/sum*100
    DF['Mg'] = DF.Mg/sum*100
    DF['Ca'] = DF.Ca/sum*100
    DF['Na'] = DF.Na/sum*100
    DF['K'] = DF.K/sum*100
    DF['P'] = DF.P/sum*100
    if 'Fe2O3' in df.columns and 'FeO' in df.columns:
        DF = DF.filter(items=['Si', 'Ti', 'Al', 'FeII','FeIII', 'Mn', 'Mg', 'Ca', 
                          'Na', 'K', 'P'])
        print('both')
    else:
        DF = DF.filter(items=['Si', 'Ti', 'Al', 'Fe','Mn', 'Mg', 'Ca', 
                          'Na', 'K', 'P'])
    return(DF)

def molar2weight(df, kwarg):
    """
    Calculate whole rock oxide weigth % from molar %. 
    Imput must be a dataframe structured like: df.Al, df. Fe, df.K
    If both Fe2O3 and FeO are present, input is: df.FeII, df.FeIII
    kwarg="FeIII" or "FeII" or "both" to choose how Fe will be expressed in wt%
    
    Output is a new dataframe structured like: df.Al2O3, df.FeO, df.K2O etc.
    If both Fe2O3 and FeO are present, output is: df.FeO, df.Fe2O3

    Additional columns are ignored.
    """
    wSi = 28.085
    wTi = 47.867 
    wAl = 26.982
    wFe = 55.845
    wMn = 54.938
    wMg = 24.305
    wCa = 40.078
    wNa = 22.990
    wK = 39.098
    wP = 30.974
    wO = 15.999

    DF = df.filter(items=['Specimen'])
    
    DF['SiO2'] = df.Si*wSi/(wSi/(wSi+2*wO))
    DF['TiO2'] = df.Ti*wTi/(wTi/(wTi+2*wO))
    DF['Al2O3'] = df.Al*wAl/(2*wAl/(2*wAl+3*wO))
    if 'FeII' in kwarg:
        DF['FeO'] = df.Fe*wFe/(wFe/(wFe+wO))
    elif 'FeIII' in kwarg:
        DF['Fe2O3'] = df.Fe*wFe/(2*wFe/(2*wFe+3*wO))
    elif 'both' in kwarg:
        DF['Fe2O3'] = df.FeIII*wFe/(2*wFe/(2*wFe+3*wO))
        DF['FeO'] = df.FeII*wFe/(wFe/(wFe+wO))
    DF['MnO'] = df.Mn*wMn/(wMn/(wMn+wO))
    DF['MgO'] = df.Mg*wMg/(wMg/(wMg+wO))
    DF['CaO'] = df.Ca*wCa/(wCa/(wCa+wO))
    DF['Na2O'] = df.Na*wNa/(2*wNa/(2*wNa+wO))
    DF['K2O'] = df.K*wK/(2*wK/(2*wK+wO))
    DF['P2O5'] = df.P*wP/(2*wP/(2*wP+5*wO))
    
    sum = DF.sum(axis=1)
    
    DF['SiO2'] = DF.SiO2/sum*100
    DF['TiO2'] = DF.TiO2/sum*100
    DF['Al2O3'] = DF.Al2O3/sum*100
    if 'FeII' in kwarg:
        DF['FeO'] = DF.FeO/sum*100
    elif 'FeIII' in kwarg:
        DF['Fe2O3'] = DF.Fe2O3/sum*100
    elif 'both' in kwarg:
        DF['Fe2O3'] = DF.Fe2O3/sum*100
        DF['FeO'] = DF.FeO/sum*100
    DF['MnO'] = DF.MnO/sum*100
    DF['MgO'] = DF.MgO/sum*100
    DF['CaO'] = DF.CaO/sum*100
    DF['Na2O'] = DF.Na2O/sum*100
    DF['K2O'] = DF.K2O/sum*100
    DF['P2O5'] = DF.P2O5/sum*100
    
    return (DF)

def AIOCG(df):
    """

    Parameters
    ----------
    df : Whole rock chemitry molar%
        e.g.
        df.Fe df.Ca etc
    Returns
    X,Y coordinates for AIOCG diagram of Montreuil et al.(2013)
    -------
    """
    x = (2*df.Ca+5*df.Fe+2*df.Mn)/(2*df.Ca+5*df.Fe+2*df.Mn+df.Mg+df.Si)
    y = df.K/(df.K+df.Na+0.5*df.Ca)
    return(x,y)

def AIOCG_Labels(df, kwarg):
    """
    input:  df: dataframe with whole rock chemistry in molar % 
            kwarg: type of bar code
                krwarg = 'Mg' for a Na-Ca-Fe-K-Mg label
                kwarg = 'SiAl' for a Na-Ca-Fe-K-(Si+Al)/10 label
    ouput: DF with one column per pie chart component. Each specimen gets an
            array with coordinates per columns.
    """
    if 'Mg' in kwarg:
        df = df.filter(items=['Na', 'Ca', 'Fe', 'K', 'Mg'])
        sum = df.sum(axis=1)
        df = df.divide(sum, axis='index') #normalize on Na+Ca+Fe+F+Mg
        df = df.multiply(100, axis='index') #normalize on Na+Ca+Fe+F+Mg
        x = []
        y = []
        xy = []
        DF = df
        df = df.astype('object')
        for i in range(len(df)):
            y0 = 0
            for r1 in df.columns:
                y_val = y0 + df.at[df.index[i], r1]
                df.at[df.index[i], r1] = [[0, y0],[0, y_val]]   # ✅ direct label-based assignment
                y0 = y_val
        
    elif 'SiAl' in kwarg:
        df['SiAl']= (df.Si+df.Al)/10
        df = df.filter(items=['Na', 'Ca', 'Fe', 'K', 'SiAl'])    
        sum = df.sum(axis=1)
        df = df.divide(sum, axis='index') #normalize on Na+Ca+Fe+F+Mg
        df = df.multiply(100, axis='index') #normalize on Na+Ca+Fe+F+Mg
        x = []
        y = []
        xy = []
        DF = df
        df = df.astype('object')
        for i in range(len(df)):
            y0 = 0
            for r1 in df.columns:
                y_val = y0 + df.at[df.index[i], r1]
                df.at[df.index[i], r1] = [[0, y0],[0, y_val]]   # ✅ direct label-based assignment
                y0 = y_val
    return(df, DF)

def plotAIOCG(df_molar, kwarg, ipath, ax):
    img = plt.imread(ipath)
    sy, sx, other = img.shape
    #adjust plotting space for margin
    #           left     right
    sx = sx - 159.412 - (sx - 1288.58)
    #           top      bottom
    sy = sy - 29.3133 - (sy - 781.424)
    
    x, y = AIOCG(df_molar)
    t, DF = AIOCG_Labels(df_molar, kwarg)
    
    # Change x,y coordinates to match the dimension of the image
    x = x*sx
    y = sy - (y*sy) #Origin at the upper left of the image
    
    # Adjust x,y for the image margins
    x = x+159.412
    y = y+29.3133
    
    size = 5000
    mw=10
    ax.imshow(img)

    for i in range(len(df_molar)):
        ax.scatter(x[i], y[i], marker=t.iloc[i].Na, 
                s=size, linewidth = mw, facecolor='hotpink')
        ax.scatter(x[i], y[i], marker=t.iloc[i].Ca, 
                s=size, linewidth = mw, facecolor='green')
        ax.scatter(x[i], y[i], marker=t.iloc[i].Fe, 
                s=size, linewidth = mw, facecolor='k')
        ax.scatter(x[i], y[i], linewidth = mw, marker=t.iloc[i].K, 
                s=size, facecolor='red')
        if 'Mg' in kwarg:
            ax.scatter(x[i], y[i], marker=t.iloc[i].Mg, 
                       s=size, linewidth = mw, facecolor='lime')
        elif 'SiAl' in kwarg:
            ax.scatter(x[i], y[i], marker=t.iloc[i].SiAl, 
                    s=size, linewidth = mw, facecolor='yellow')    
    return(ax)
#%%
ipath = "C:\\AdjustPath\\Park2014.csv"

df = pd.read_csv(ipath, index_col=0)
df = df.set_index('Specimen')
#%%
#Cleanup
df = df.filter(items=['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'MnO', 'MgO', 'CaO', 
                      'Na2O', 'K2O','P2O5', 'Specimen',  'Unit'
                      ,'MetZone'
                      ])
#%%
df_molar = weigth2molar(df)

ipath = "C:\\AdjustPath\\AIOCG.png"

fig, (ax, ax2) = plt.subplots(1,2, figsize=(20, 10)) 
plotAIOCG(df_molar, 'Mg', ipath, ax)
plotAIOCG(df_molar, 'SiAl', ipath, ax2)

ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax2.axes.xaxis.set_ticks([])
ax2.axes.yaxis.set_ticks([])
plt.show()

#ipath = "AdjustPath\\Park2014.png"
#plt.savefig(ipath)

