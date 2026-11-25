"""Publication-style, accurately labeled exploratory-in-silico charts from frozen mechanistic results."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scripts.gcev_mechanism_factorial import bootstrap

ORDER=['linear','saturation','interaction','combined_mechanism']
LABEL=['Linear','Saturation','Pairwise interaction','Combined']

def export(fig,out,name):
    fig.savefig(out/f'{name}.png',dpi=300,bbox_inches='tight')
    fig.savefig(out/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

def plot(root):
    root=Path(root);out=root/'figures';out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(root/'per_dag_mechanism_comparisons.csv')
    p=pd.read_csv(root/'per_dag_factorial_effects.csv')
    diag=pd.read_csv(root/'per_condition_source_decisions.csv')
    clean=d[d.scenario=='clean'];x=np.arange(4)
    fig,ax=plt.subplots(figsize=(9,5.2));rng=np.random.default_rng(142)
    for i,name in enumerate(ORDER):
        q=clean[clean.mechanism==name];v=q.dAP_graph.to_numpy();strata=q.n_genes.astype(str)+'_'+q.evidence_regime
        b=bootstrap(v,strata.to_numpy())
        ax.scatter(i+rng.uniform(-.18,.18,len(v)),v,s=14,alpha=.35)
        ax.errorbar(i,b['mean'],yerr=[[b['mean']-b['low95']],[b['high95']-b['mean']]],
                    fmt='o',capsize=5,markersize=7,linewidth=1.8,color='black')
    ax.axhline(0,color='gray',linestyle='--',linewidth=.9)
    ax.set_xticks(x,LABEL);ax.set_ylabel('GCEV − graph-independent AP')
    ax.set_title('Graph-state contribution across matched mechanisms')
    ax.set_ylim(min(-.24,clean.dAP_graph.min()-.03),max(.22,clean.dAP_graph.max()+.03))
    ax.text(.01,.02,'32 paired DAGs per mechanism; dots = DAGs; black = mean ± 95% bootstrap CI\nSynthetic oracle-derived evidence — not biological validation',transform=ax.transAxes,fontsize=8)
    fig.tight_layout();export(fig,out,'F01_GCEV_GRAPH_STATE_AP_BY_MECHANISM')

    fig,ax=plt.subplots(figsize=(9,5.2))
    cols=['ap_fixed_sources','ap_no_graph_state','ap_gcev']
    captions=['Fixed evidence','Adaptive, no graph state','GCEV']
    for j,col in enumerate(cols):
        m=[clean.loc[clean.mechanism==name,col].mean() for name in ORDER]
        ax.bar(x+(j-1)*.24,m,width=.23,label=captions[j])
    ax.set_xticks(x,LABEL);ax.set_ylabel('Mean directed-edge average precision (AP)')
    ax.set_title('Fixed weighting remains competitive under oracle-derived evidence')
    ax.set_ylim(0,.75);ax.legend(loc='upper right')
    fig.tight_layout();export(fig,out,'F02_ABSOLUTE_AP_BASELINES')

    fig,ax=plt.subplots(figsize=(9,5.2));labels=[]
    for i,(stratum,g) in enumerate(p.groupby('stratum',sort=True)):
        ax.scatter(np.full(len(g),i)+rng.uniform(-.16,.16,len(g)),g.theta_interaction,alpha=.55,s=30)
        ax.plot([i-.23,i+.23],[g.theta_interaction.mean()]*2,color='black',linewidth=2)
        labels.append(stratum.replace('_','\n'))
    b=bootstrap(p.theta_interaction,p.stratum.to_numpy());ax.axhline(b['mean'],linestyle='-',label=f"Pooled paired mean {b['mean']:+.4f}")
    ax.axhline(0,linestyle='--',color='gray')
    ax.set_xticks(range(len(labels)),labels);ax.set_ylabel('Paired interaction effect on incremental AP')
    ax.set_title('The preregistered interaction hypothesis is not supported')
    ax.legend(loc='upper right');fig.tight_layout();export(fig,out,'F03_PAIRED_INTERACTION_EFFECT_PER_DAG')

    fig,ax=plt.subplots(figsize=(9,5.2))
    for j,scenario in enumerate(['clean','source_stress']):
        vals=[diag.loc[(diag.mechanism==name)&(diag.scenario==scenario),'changed_edit_fraction'].mean() for name in ORDER]
        ax.bar(x+(j-.5)*.35,np.array(vals)*100,width=.34,label=scenario.replace('_',' '))
    ax.set_xticks(x,LABEL);ax.set_ylabel('Paired source-ablation decisions changing an edge (%)')
    ax.set_ylim(0,30);ax.set_title('Source evidence changes relatively few one-edit choices')
    ax.legend();fig.tight_layout();export(fig,out,'F04_SOURCE_INDUCED_DECISION_CHANGES')

    fig,ax=plt.subplots(figsize=(9,5.2))
    for i,name in enumerate(ORDER):
        q=clean[clean.mechanism==name]
        strata=(q.n_genes.astype(str)+'_'+q.evidence_regime).to_numpy()
        b=bootstrap(q.dNMSE_graph.to_numpy(),strata)
        ax.errorbar(i,b['mean']*100,yerr=[[(b['mean']-b['low95'])*100],[(b['high95']-b['mean'])*100]],
                    fmt='o',capsize=5,markersize=7,color='black')
    ax.axhline(0,color='gray',linestyle='--');ax.set_xticks(x,LABEL)
    ax.set_ylabel('Relative untouched-test NMSE: GCEV − no graph state (%)')
    ax.set_title('Generalization trade-off across paired mechanisms')
    fig.tight_layout();export(fig,out,'F05_UNTOUCHED_EXPRESSION_NMSE')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path)
    plot(p.parse_args().root)
