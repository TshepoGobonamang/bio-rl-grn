"""Generate descriptive figures from predeclared network-level bootstrap outputs.

Figures show synthetic-only estimation; not independent biological confirmation.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCENARIOS=['clean','dropout50','corrupt25','combined']
LABELS=['Clean evidence','50% missing','25% inverted scores','Missing + inverted']

def main(root:Path,out:Path):
    data=pd.read_csv(root/'network_bootstrap_95.csv')
    paired=pd.read_csv(root/'paired_network_effects.csv')
    out.mkdir(parents=True,exist_ok=True)
    for metric,stem,title in [('delta_ap_vs_no_graph_state','F01_GRAPH_STATE_INCREMENTAL_AP','GCEV vs graph-independent adaptive weighting'),
                              ('delta_ap_vs_fixed_sources','F02_GCEV_VS_FIXED_AP','GCEV vs fixed source weighting')]:
        rows=data[data.comparison==metric].set_index('scenario').loc[SCENARIOS]
        y=rows['mean'].to_numpy(float)
        lo=y-rows['ci95_low'].to_numpy(float)
        hi=rows['ci95_high'].to_numpy(float)-y
        fig,ax=plt.subplots(figsize=(8.5,4.8))
        xx=np.arange(4)
        ax.errorbar(xx,y,yerr=[lo,hi],fmt='o',capsize=5,linewidth=1.4)
        ax.axhline(0,linewidth=.9,linestyle='--')
        ax.set_xticks(xx,LABELS,rotation=13,ha='right')
        ax.set_ylabel('Mean paired AP difference')
        ax.set_title(title+'\n48 independently generated DAGs; synthetic only')
        ax.grid(axis='y',alpha=.2)
        fig.tight_layout()
        for ext in ['png','svg','pdf']:
            fig.savefig(out/f'{stem}.{ext}',dpi=220)
        plt.close(fig)
    clean=paired[paired.scenario=='clean']
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    rng=np.random.default_rng(20261125)
    fams=['linear','nonlinear']
    for i,fam in enumerate(fams):
        vals=clean[clean.family==fam].delta_ap_vs_no_graph_state.to_numpy(float)
        ax.scatter(rng.normal(loc=i,scale=.06,size=len(vals)),vals,alpha=.63,marker='o')
        ax.scatter([i],[vals.mean()],marker='_',s=350,linewidths=2)
    ax.set_xticks([0,1],['Linear synthetic DAGs\n24 networks','Nonlinear synthetic DAGs\n24 networks'])
    ax.set_ylabel('GCEV − no_graph_state AP per network')
    ax.set_title('Graph-state incremental value varies by graph generator\nUnplanned subgroup interpretation; synthetic only')
    ax.axhline(0,linewidth=.9,linestyle='--')
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout()
    for ext in ['png','svg','pdf']:
        fig.savefig(out/f'F03_GENERATOR_HETEROGENEITY.{ext}',dpi=220)
    plt.close(fig)

if __name__=='__main__':
    a=argparse.ArgumentParser()
    a.add_argument('--results',type=Path,required=True)
    a.add_argument('--output',type=Path,required=True)
    b=a.parse_args();main(b.results,b.output)
