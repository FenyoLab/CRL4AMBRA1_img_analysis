#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 15:04:50 2017

@author: keriabermudez and sarahkeegan
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']

mpl.rcParams['xtick.labelsize'] = 7 
mpl.rcParams['ytick.labelsize'] = 7 
mpl.rcParams['axes.linewidth'] = 1 

mpl.rcParams['font.size'] = 7

sns.set()
sns.set_context("paper")
sns.set_style("white")
sns.set_style("ticks",{"xtick.major.size": 4, "ytick.major.size": 4,"ytick.minor.direction": 'in'})

import numpy as np

#%%
def plot_grid(table, column):
    kws = dict(s=1)
    g = sns.FacetGrid(table, col="Sample_Number",hue ='Label',col_wrap=3, height =2.5)
    g.map(plt.scatter, "offset_z", column, alpha=0.5,**kws)
    return g

#%%    
def plot_combined(table,column):
    g = sns.lmplot(x="offset_z", y=column, hue="Sample_Number_Label", data=table,fit_reg = False, height=3, legend = True,scatter_kws={"s": 1});
    return g

#%%
parent_path = "./PCNA_CYD1_Videos/"
folder_cd1 = "/results_cd1_Blad1/"
folder_pcna = "/results_pcna_Blad1/"
        
#%%
#### CD1 WT ####
    
min_consec_hour = 10
min_hour = -1 #-5 #-1
max_hour = 12 #24

group = "WT"

path =parent_path + group+folder_cd1

table_wt = pd.read_csv(path+"tracked_cells.csv")

table_wt['Sample_Number_Label'] =  table_wt['Sample_Number'].astype('str') +'_'+table_wt['Label'].astype('str')

table_wt['group'] = group
table_wt['hours'] = table_wt['offset_z']*10/60

path_table = parent_path+group+ "/"

xlsx = pd.ExcelFile(path_table + "WT_table.xlsx")
wt_tracking = pd.read_excel(xlsx, 'Blad1')

wt_tracking['Sample_Number_Label'] =  wt_tracking['Sample_Number'].astype('str') +'_'+wt_tracking['Label'].astype('str')
wt_tracking_well = wt_tracking[wt_tracking.ix[:,'Tracked Well'] == True]

wt_to_analyze =  wt_tracking_well.Sample_Number_Label.unique()

table_wt =table_wt[table_wt.Sample_Number_Label.isin(wt_to_analyze)]

wt_hours = table_wt.groupby(by='Sample_Number_Label').max()

# filter table for all the cells tracked from 0 to more than min_consec_hour

wt_to_analyze =  wt_hours[wt_hours.hours > min_consec_hour].index

table_wt =table_wt[table_wt.Sample_Number_Label.isin(wt_to_analyze)]

nuclei_ids = table_wt['Sample_Number_Label'].unique()
table_wt['intensity_change']=0
for nuclei_id in nuclei_ids:
    intensity_at_t0 = table_wt[(table_wt['offset_z']==0) & (table_wt['Sample_Number_Label']==nuclei_id)].iloc[0]['cyd1_mean_intensity']
    table_wt['intensity_change'] = np.where(table_wt['Sample_Number_Label']==nuclei_id,(table_wt['cyd1_mean_intensity']-intensity_at_t0)/intensity_at_t0,table_wt['intensity_change'])
    table_wt['perc_intensity_change'] = table_wt['intensity_change'] * 100

#%%
#### CD1 KO ####

group = "KO"
path = parent_path+group+folder_cd1

table_ko = pd.read_csv(path+"tracked_cells.csv")
table_ko['Sample_Number_Label'] =  table_ko['Sample_Number'].astype('str') +'_'+table_ko['Label'].astype('str')

table_ko['group'] = group
table_ko['hours'] = table_ko['offset_z']*10/60

path_table = parent_path+group+ "/"

xlsx = pd.ExcelFile(path_table + "KO_table.xlsx")
ko_tracking = pd.read_excel(xlsx, 'Blad1')

ko_tracking['Sample_Number_Label'] =  ko_tracking['Sample_Number'].astype('str') +'_'+ko_tracking['Label'].astype('str')
ko_tracking_well = ko_tracking[ko_tracking.ix[:,'Tracked Well'] == True]

ko_to_analyze =  ko_tracking_well.Sample_Number_Label.unique()

table_ko = table_ko[table_ko.Sample_Number_Label.isin(ko_to_analyze)]

ko_hours = table_ko.groupby(by='Sample_Number_Label').max()

# filter table for all the cells tracked from 0 to more than min_consec_hour

ko_to_analyze =  ko_hours[ko_hours.hours > min_consec_hour].index

table_ko =table_ko[table_ko.Sample_Number_Label.isin(ko_to_analyze)]

nuclei_ids = table_ko['Sample_Number_Label'].unique()
table_ko['intensity_change']=0
for nuclei_id in nuclei_ids:
    intensity_at_t0 = table_ko[(table_ko['offset_z']==0) & (table_ko['Sample_Number_Label']==nuclei_id)].iloc[0]['cyd1_mean_intensity']
    table_ko['intensity_change'] = np.where(table_ko['Sample_Number_Label']==nuclei_id,(table_ko['cyd1_mean_intensity']-intensity_at_t0)/intensity_at_t0,table_ko['intensity_change'])
    table_ko['perc_intensity_change'] = table_ko['intensity_change'] * 100

#%%

fig ,ax1 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300)

both_tables = [table_ko, table_wt]
both_tables = pd.concat(both_tables)
both_tables = both_tables[(both_tables.hours > min_hour) & (both_tables.hours <=max_hour)] # SK: changed from -1, 10 to min_hour,max_hour

sns.tsplot(data=both_tables, err_style="unit_traces",interpolate = True, time="hours", unit="Sample_Number_Label", value="cyd1_mean_intensity" , condition = "group", ax = ax1 )

fig.savefig(parent_path+'CYD1_unit_traces.pdf')

#%%
# Filtering wt table where the hours are greater than -1 and less than 24 hours

table_wt = table_wt[(table_wt.hours > min_hour) & (table_wt.hours <max_hour)]

wt_median = table_wt.groupby(by = 'hours')['cyd1_mean_intensity'].median()
wt_mean =  table_wt.groupby(by = 'hours')['cyd1_mean_intensity'].mean()
wt_std = table_wt.groupby(by = 'hours')['cyd1_mean_intensity'].std()

wt_quantile_25 = table_wt.groupby(by = 'hours')['cyd1_mean_intensity'].quantile(q = 0.25)
wt_quantile_75 = table_wt.groupby(by = 'hours')['cyd1_mean_intensity'].quantile(q = 0.75)

table_ko = table_ko[(table_ko.hours > min_hour) & (table_ko.hours < max_hour)]

ko_median = table_ko.groupby(by = 'hours')['cyd1_mean_intensity'].median()
ko_mean =  table_ko.groupby(by = 'hours')['cyd1_mean_intensity'].mean()
ko_std = table_ko.groupby(by = 'hours')['cyd1_mean_intensity'].std()

ko_quantile_25 = table_ko.groupby(by = 'hours')['cyd1_mean_intensity'].quantile(q = 0.25)
ko_quantile_75 = table_ko.groupby(by = 'hours')['cyd1_mean_intensity'].quantile(q = 0.75)


# TODO - try plotting percent change from time==0 for these?
# TODO - shift KO down so mean(KO) == mean(WT) at time == 0, and show 95% C.I.

## Median Plot
fig ,ax1 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))

ax1.plot(wt_median.index, wt_median, 'b', label = 'WT ' +'n = '+str(len(wt_to_analyze)))
ax1.plot(ko_median.index, ko_median, 'g', label = 'KO ' +'n = '+str(len(ko_to_analyze)))
ax1.fill_between(wt_median.index, wt_quantile_25,wt_quantile_75 , color='b', alpha=0.2)
ax1.fill_between(ko_median.index, ko_quantile_25,ko_quantile_75 , color='g', alpha=0.2)
#ax1.set_xticks([-120,0,120,300,450,600,750])
ax1.set_xlabel('time after foci formation (hours)')
ax1.set_title('Median Nuclear CYD1')
ax1.set_ylabel('nuclear CYD1 (a.u.)')
ax1.legend(loc='best')
fig.tight_layout()
fig.savefig(parent_path+'Median_CYD1.pdf')

## Mean Plot
fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))
ax2.plot(wt_mean.index, wt_mean, 'b', label = 'WT ' +'n = '+str(len(wt_to_analyze)))
ax2.plot(ko_mean.index, ko_mean, 'g', label = 'KO '+'n = '+str(len(ko_to_analyze)))
ax2.fill_between(wt_std.index,wt_mean- wt_std,wt_mean+wt_std , color='b', alpha=0.2)
ax2.fill_between(ko_std.index, ko_mean- ko_std,ko_mean+ko_std , color='g', alpha=0.2)
#ax2.set_xticks([-150,0,150,300,450,600,750])
ax2.set_xlabel('time after foci formation (hours)')
ax2.set_ylabel('nuclear CYD1 (a.u.)')
ax2.set_title('Mean Nuclear CYD1')
ax2.legend(loc='best')
fig2.tight_layout()
fig2.savefig(parent_path+'Mean_CYD1.pdf')
fig2.clf()

# Mean Plot - KO shifted down
wt_mean_at_t0 = wt_mean.loc[0]
ko_mean_at_t0 = ko_mean.loc[0]

ko_norm_mean = ko_mean - (ko_mean_at_t0-wt_mean_at_t0)

fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))
ax2.plot(wt_mean.index, wt_mean, 'b', label = 'WT ' +'n = '+str(len(wt_to_analyze)))
ax2.plot(ko_norm_mean.index, ko_norm_mean, 'g', label = 'KO '+'n = '+str(len(ko_to_analyze)))
ax2.fill_between(wt_std.index,wt_mean- wt_std,wt_mean+wt_std , color='b', alpha=0.2)
ax2.fill_between(ko_std.index, ko_norm_mean-ko_std, ko_norm_mean+ko_std , color='g', alpha=0.2)
#ax2.set_xticks([-150,0,150,300,450,600,750])
ax2.set_xlabel('time after foci formation (hours)')
ax2.set_ylabel('nuclear CYD1 (a.u.)')
ax2.set_title('Mean Nuclear CYD1')
ax2.legend(loc='best')
fig2.tight_layout()
fig2.savefig(parent_path+'Mean_CYD1_KO_shifted.pdf')
fig2.clf()

#Mean plot - CI
#use seaborn lineplot
#need the raw intensities

fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))

both_tables['norm_cyd1_mean_intensity']=both_tables['cyd1_mean_intensity']
both_tables['norm_cyd1_mean_intensity'] = np.where(both_tables['group']=='KO', both_tables['cyd1_mean_intensity']-(ko_mean_at_t0-wt_mean_at_t0), both_tables['cyd1_mean_intensity'])

sns.lineplot(x='hours',y='norm_cyd1_mean_intensity', data=both_tables, hue='group', ax=ax2)

ax2.set_xlabel('time after foci formation (hours)')
ax2.set_ylabel('nuclear CYD1 (a.u.)')
ax2.set_title('Mean Nuclear CYD1')
ax2.legend(loc='best')
fig2.tight_layout()
fig2.savefig(parent_path+'Mean_CYD1_KO_shifted_CI.pdf')
fig2.clf()

#mean plot + sd
fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))

both_tables['norm_cyd1_mean_intensity']=both_tables['cyd1_mean_intensity']
both_tables['norm_cyd1_mean_intensity'] = np.where(both_tables['group']=='KO', both_tables['cyd1_mean_intensity']-(ko_mean_at_t0-wt_mean_at_t0), both_tables['cyd1_mean_intensity'])

sns.lineplot(x='hours',y='norm_cyd1_mean_intensity', data=both_tables, hue='group', ci='sd', legend=False, ax=ax2)

ax2.set_xlabel('time after foci formation (hours)')
ax2.set_ylabel('nuclear CYD1 (a.u.)')
ax2.set_title('Mean Nuclear CYD1')
ax2.legend(loc='best')
fig2.tight_layout()
fig2.savefig(parent_path+'Mean_CYD1_KO_shifted_SD.pdf')
fig2.clf()

#change in intensity - CI
fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))


sns.lineplot(x='hours',y='perc_intensity_change', data=both_tables, hue='group', legend=False, ax=ax2)

children = ax2.get_children()
polys={}
polys['KO'] = children[0]
polys['WT'] = children[1]

ci_table={}
for group in polys.keys():
    ci_table[group] = {}

    paths = polys[group].get_paths()
    p1=paths[0]
    len_ = len(p1.vertices)

    first_len = int(len_/2)
    v1 = p1.vertices[1:first_len]
    v2 = p1.vertices[first_len+1:-1]
    v2=v2[::-1]

    for v_i,v in enumerate(v1):

        ci_table[group][v[0]]=[v[1],v2[v_i][1]]

ax2.set_xlabel('Time After Foci Formation (hours)')
ax2.set_ylabel('Change in Mean Nuclear CYD1 (%)')
#ax2.set_title('Change in Mean Nuclear CYD1')
#ax2.legend(loc='best')
ax2.set_xlim(0,12)
fig2.tight_layout()
fig2.savefig(parent_path+'Change_in_Mean_CYD1_CI.pdf')
fig2.clf()

#change in intensity - sd
fig2 ,ax2 = plt.subplots(nrows= 1,ncols= 1, figsize = (3,3), dpi=300) #(2.0,3))

sns.lineplot(x='hours',y='perc_intensity_change', data=both_tables, hue='group', ci='sd', ax=ax2)

ax2.set_xlabel('Time After Foci Formation (hours)')
ax2.set_ylabel('Change in Mean Nuclear CYD1 (%)')
#ax2.set_title('Change in Mean Nuclear CYD1')
#ax2.legend(loc='best')
ax2.set_xlim(0,12)
fig2.tight_layout()
fig2.savefig(parent_path+'Change_in_Mean_CYD1_SD.pdf')
fig2.clf()

##########
#save raw data for Daniele to include in the paper
both_tables.to_csv(parent_path+'All_Data_Change_in_Mean_CYD1.txt', columns=['offset_z', 'hours','group','Sample_Number_Label', 'cyd1_mean_intensity','perc_intensity_change'], sep='\t') #,index=False)
wt_to_save = pd.DataFrame(columns=['hours','mean','ci_high','ci_low'])
for time in table_wt['hours'].unique():
    if(time in ci_table['WT']):
        cur_table = table_wt[table_wt['hours'] == time]
        new_row = {}
        new_row['hours'] = time
        new_row['mean'] = cur_table['perc_intensity_change'].mean()
        new_row['group'] = 'WT'
        new_row['ci_high'] = ci_table['WT'][time][1]
        new_row['ci_low'] = ci_table['WT'][time][0]

        wt_to_save = wt_to_save.append(new_row, ignore_index=True)

ko_to_save = pd.DataFrame(columns=['hours','mean','ci_high','ci_low'])
for time in table_ko['hours'].unique():
    if(time in ci_table['KO']):
        cur_table = table_ko[table_ko['hours'] == time]
        new_row = {}
        new_row['hours'] = time
        new_row['mean'] = cur_table['perc_intensity_change'].mean()
        new_row['group']='KO'
        new_row['ci_high'] = ci_table['WT'][time][1]
        new_row['ci_low'] = ci_table['WT'][time][0]

        ko_to_save = ko_to_save.append(new_row, ignore_index=True)

both_to_save = pd.concat([wt_to_save,ko_to_save])
both_to_save.index = range(len(both_to_save))

both_to_save.to_csv(parent_path+'Plot_Data_Change_in_Mean_CYD1_CI.txt', sep='\t')

#%%
### Individual Plots for WT
group = "WT"
g = plot_grid(table_wt, "cyd1_mean_intensity")
g.savefig(parent_path+group +"_cyd1_mean_intensity.pdf")

g = plot_combined(table_wt,"cyd1_mean_intensity")
g.savefig(parent_path+group +"_combined_cyd1_mean_intensity.pdf")


#%%
### Individual Plots for KO

group = "KO"
g = plot_grid(table_ko, "cyd1_mean_intensity")
g.savefig(parent_path+group +"_cyd1_mean_intensity.pdf")

g = plot_combined(table_ko,"cyd1_mean_intensity")
g.savefig(parent_path+group +"_combined_cyd1_mean_intensity.pdf")



#%%
#### PCNA WT ####
### Individual Plots for WT

group = "WT"
path = parent_path+ "/"+ group+folder_pcna
table = pd.read_csv(path+"tracked_cells.csv")
table['Sample_Number_Label'] =  table['Sample_Number'].astype('str') +'_'+table['Label'].astype('str')

table = table[table.Sample_Number_Label.isin(wt_to_analyze)]

g = plot_grid(table, "pcna_mean_intensity")
g.savefig(parent_path+group +"_pcna_mean_intensity.pdf")

g = plot_combined(table,"pcna_mean_intensity")
g.savefig(parent_path+group +"_combined_pcna_mean_intensity.pdf")


#%%
#### PCNA KO ####
### Individual Plots for KO

group = "KO"
path = parent_path+ "/"+ group+folder_pcna
table = pd.read_csv(path+"tracked_cells.csv")
table['Sample_Number_Label'] =  table['Sample_Number'].astype('str') +'_'+table['Label'].astype('str')

table[table.Sample_Number_Label.isin(ko_to_analyze)]

g = plot_grid(table, "pcna_mean_intensity")
g.savefig(parent_path+group+"_pcna_mean_intensity.pdf")

g = plot_combined(table,"pcna_mean_intensity")
g.savefig(parent_path+group+"_combined_pcna_mean_intensity.pdf")


