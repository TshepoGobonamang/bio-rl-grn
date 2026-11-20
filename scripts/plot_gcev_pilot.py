"""Visualize observed controlled synthetic pilot only; never infer real biology."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot(source:Path,output:Path):
    frame=pd.read_csv(source)
    chosen=['conditional','fixed_sources','global_source','no_graph_state','gcev']
    table=frame.loc[frame.method.isin(chosen)].pivot(index='seed',columns='method',values='ap')
    table=table[chosen]
    output.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(9.4,4.7),constrained_layout=True)
    pos=np.arange(len(table))
    width=.16
    for i,method in enumerate(chosen):
        ax.bar(pos+(i-2)*width,table[method],width,label=method.replace('_',' '))
    ax.set_xticks(pos,[f'Network seed {s}' for s in table.index])
    ax.set_ylabel('Average precision (noninterpolated)')
    ax.set_ylim(0,1)
    ax.set_title('Synthetic pilot: GCEV against fixed and graph-independent source weighting')
    ax.grid(axis='y',alpha=.2)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),ncol=3,fontsize=8,frameon=False)
    fig.text(.5,.005,'ENGINEERING PILOT ONLY: simulated graph and evidence generated using ground truth; not independent biology',ha='center',fontsize=8)
    for suffix in ('png','svg','pdf'):
        fig.savefig(output/f'GCEV_AP_PILOT_NOT_BIOLOGICAL.{suffix}',dpi=240)
    plt.close(fig)

    loss=frame.loc[frame.method.isin(['fixed_sources','gcev','no_graph_state'])].pivot(index='seed',columns='method',values='test_nmse')
    fig,ax=plt.subplots(figsize=(8.3,4.4),constrained_layout=True)
    pos=np.arange(len(loss));width=.23
    for i,method in enumerate(loss.columns):
        ax.bar(pos+(i-1)*width,loss[method],width,label=method.replace('_',' '))
    ax.set_xticks(pos,[f'Network seed {s}' for s in loss.index])
    ax.set_ylabel('Untouched expression-test NMSE (lower is better)')
    ax.set_title('Synthetic pilot: independent test-expression generalization')
    ax.legend(fontsize=9,frameon=False)
    ax.grid(axis='y',alpha=.2)
    fig.text(.5,.005,'ENGINEERING PILOT ONLY — three generated networks; no biological assays',ha='center',fontsize=8)
    for suffix in ('png','svg','pdf'):
        fig.savefig(output/f'GCEV_TEST_LOSS_PILOT_NOT_BIOLOGICAL.{suffix}',dpi=240)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--metrics',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args();plot(a.metrics,a.out)
