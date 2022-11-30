## Do the imports

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


## define directories

outfold = '/Users/bn/Desktop/bhack_td-main/results'


## load data into 3D array: models *cross-vals * channels

raw = scipy.io.loadmat(outfold + '/bhack_task_04_output_GLMs.mat')
data = np.squeeze(np.append(raw['Stats'][0,0], raw['Stats'][1,0], axis=1))
print(data.shape)


## convert 3D array to DataFrame (+ change order of dim)

df = pd.DataFrame(columns=['model','channel','cross_vals'])

for i in range( data.shape[0] ):
    for j in range( data.shape[2] ):
        df_tmp = pd.DataFrame(data=data[i,:,j], columns=['cross_vals'])
        df_tmp['model'] = i
        df_tmp['channel'] = j
        df = pd.concat((df, df_tmp), ignore_index=True)


## plot split violin plot

fig, ax = plt.subplots(sharey=True, facecolor='w', figsize=(5, 3)) # dpi=300, 

sns.violinplot(data=df, x='channel', y='cross_vals', hue='model', split=True, 
               inner= 'box', bw=0.5, gridsize=55, linewidth=1, linecolor='black')
sns.swarmplot (data=df, x='channel', y='cross_vals', hue='model', 
               dodge=True, size=3, edgecolor='black', linewidth=0.5, alpha=1);

ax.axhline(0, color='black', ls='--', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([-.1, .25])
ax.set_yticks([-.1, 0, .1, .2, .3])

ax.set_xlabel('iEEG channels')
ax.set_ylabel('R2_CV');
ax.tick_params(axis='both', labelsize=10, which='both')
ax.legend().set_visible(False)
ax.set_title('The decoding ultimate analysis'); # fontweight='bold'


## 1. non-split: only 1 model - entire pipeline

# load the data
raw = scipy.io.loadmat(outfold + '/bhack_task_04_output_GLMs.mat')
# print(raw.keys()) # or: for key in raw.keys(): print(key)
data = [raw[key] for key in raw.keys()]
data = data[-1].copy() # generic
#data[0,0][0,0,:,:].shape

# extract data of interest
y = data[0,0][0,0,:,:]
y2 = data[1,0][0,0,:,:]
x = np.nonzero(y)[-1]
#out1.shape
y = y.flatten()
y2 = y2.flatten()

# figure
fig, ax = plt.subplots(sharey=True, facecolor='w', figsize=(5, 3)) # dpi=300, 

sns.violinplot(x=x, y=y2, inner= 'box', bw=.2, gridsize=55, linewidth=2, linecolor='black', axis=ax) #split=True
sns.swarmplot (x=x, y=y2, color="black", size=3, edgecolor='black', linewidth=0.5, alpha=0.5)

ax.axhline(0, color='black', ls='--', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim([-.05, .25]);
ax.set_xlabel('iEEG channels')
ax.set_ylabel('R2_CV');
ax.tick_params(axis='both', labelsize=10, which='both')
ax.set_title('Only one model'); # fontweight='bold'

fig.savefig(outfold + '/RTviolin_split.png', format='png', dpi=1000)

##### Examples

# non split violin plot (Jeremy Giroud)
import pandas as pd
import numpy as np
import os
from natsort import natsorted 
import string
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

data = pd.read_csv("C:\\Users\\Jeremy\\Documents\\enfin_channel_cap\\exp_2&ctrl.csv")

exp_sent = data.loc[data['experiment']=='exp_2']['scores'].values
exp_ctrl = data.loc[data['experiment']=='exp_ctrl']['scores'].values
exp_ctrl = exp_ctrl[~np.isnan(exp_ctrl)]

fig,ax = plt.subplots(1,1, facecolor='w', dpi=300)

sns.violinplot(x="experiment", y="scores",inner= 'box',bw = 0.8, color= "white",data=data,gridsize = 55, linewidth = 5,linecolor = 'black',axis =ax)
#sns.violinplot(x="experiment", y="scores",bw = 0.6,inner = "points", data=raw_df,gridsize = 25, linewidth = 2, saturation = 0.75,axis =ax)
ax.tick_params(axis='both', labelsize=15, which='both')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0,1);

sns.swarmplot(x="experiment", y="scores", size=12, linewidth= 1,data=data, ax=ax);
ax.collections[4].set_alpha(0.75)
ax.collections[5].set_alpha(0.75)
ax.collections[2].set_color("navajowhite")

ax.collections[2].set_edgecolor("#ff7f0e")
ax.collections[0].set_color("lightblue")
ax.collections[0].set_edgecolor("#1f77b4")
ax.collections[0].set_alpha(0.75)
ax.collections[2].set_alpha(0.75)
ax.set_ylabel("Performance",fontsize = 20);
ax.set_xlabel("");
ax.set_xticklabels({"Experiment 2\n (gate n°4)","Experiment 3\n (gate n°1)"});
ax.annotate('n.s.', (.5, 1), va='center', ha='center', alpha=1,fontsize = 20);
ax.annotate('______', (.485, .98), va='center', ha='center', alpha=1,fontsize = 20);


# split violin plots (J Duprez)

# create a copy of the data to put IKI in hz
RTmat3 = RTmat2.copy()
RTmat3 = RTmat3[RTmat3.Precision != 'Correct trials']
RTmat4 = RTmat3.copy()
ax.ylim = ([100, 1500])

colors = ['#1f77b4', '#9467bd', '#ff7f0e', '#d62728']

plt.subplot(1, 2, 1)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

ax = [None] * (1 + 1)
fig = plt.figure(figsize=(12, 5))
gs = mpl.gridspec.GridSpec(nrows=1,
                           ncols=2,
                           figure=fig,
                           wspace=0.15, hspace=0.05
                           )

ax[0] = fig.add_subplot(gs[0, 0])
sns.violinplot(x="Condition", y="RT", data=RTmat2[RTmat2.Precision == 'Correct trials'][RTmat2.RT != 0])
sns.swarmplot(x='Condition', y='RT', data=RTmat2[RTmat2.Precision == 'Correct trials'][RTmat2.RT != 0], color="black",
              size=3, edgecolor='black', linewidth=0.5, alpha=0.5)

ax[0].collections[0].set_facecolor('tab:blue')
ax[0].collections[2].set_facecolor('tab:purple')
ax[0].collections[4].set_facecolor('tab:orange')
ax[0].collections[6].set_facecolor('tab:red')
ax[0].set_ylim([0, 2000])
ax[0].set_yticks([0, 500, 1000, 1500, 2000])
ax[0].set_ylabel('RT (ms)')
plt.xticks(fontsize=9)
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[0].set_xlabel('')
ax[0].set_title('Correct trials', fontweight='bold')


# split violin plots
ax[1] = fig.add_subplot(gs[0, 1])
sns.violinplot(x="Condition", y="RT", hue="Precision", data=RTmat3[RTmat3.RT != 0], split=True)
sns.swarmplot(x='Condition', y='RT', hue="Precision", data=RTmat3[RTmat3.RT != 0], dodge=True, color="white", size=3,
              edgecolor='black', linewidth=0.5, alpha=0.5)

ax[1].spines['top'].set_visible(False)
ax[1].spines['left'].set_visible(False)
ax[1].axes.get_yaxis().set_visible(False)
ax[1].legend().set_visible(False)
plt.xticks(fontsize=9)
ax[1].collections[0].set_facecolor('tab:blue')
ax[1].collections[1].set_facecolor('lightsteelblue')
ax[1].collections[3].set_facecolor('tab:purple')
ax[1].collections[4].set_facecolor('plum')
ax[1].collections[6].set_facecolor('tab:orange')
ax[1].collections[7].set_facecolor('wheat')
ax[1].collections[9].set_facecolor('tab:red')
ax[1].collections[10].set_facecolor('salmon')

ax[1].set_xlabel('')
ax[1].set_title('Corrected errors (darker colors)\nand other errors (lighter colors)', fontweight='bold')

# fig.savefig(outfold + '/RTviolin_split.png', format='png', dpi=1000)

# -*- coding: utf-8 -*-

