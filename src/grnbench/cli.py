"""Command-line interface for clean benchmark and explicitly separate gold evaluation."""
import argparse
from pathlib import Path
from .core import run_fit, evaluate, generate_smoke


def main(argv=None):
    p=argparse.ArgumentParser(prog='grnbench')
    cmd=p.add_subparsers(dest='command',required=True)
    g=cmd.add_parser('generate-smoke')
    g.add_argument('--output',type=Path,required=True)
    g.add_argument('--seed',type=int,default=42)
    g.add_argument('--cells',type=int,default=180)
    g.add_argument('--genes',type=int,default=8)
    f=cmd.add_parser('fit')
    f.add_argument('--expression',type=Path,required=True)
    f.add_argument('--output',type=Path,required=True)
    f.add_argument('--tf-list',type=Path)
    f.add_argument('--allow-all-genes',action='store_true')
    f.add_argument('--prior',type=Path)
    f.add_argument('--prior-secondary',type=Path)
    f.add_argument('--gcev',action='store_true',help='Fit nested paired evidence-value gates (experimental)')
    f.add_argument('--seed',type=int,default=42)
    f.add_argument('--max-genes',type=int,default=12)
    f.add_argument('--max-parents',type=int,default=2)
    f.add_argument('--top-k',type=int,default=10)
    f.add_argument('--alpha',type=float,default=1.)
    f.add_argument('--query-budget',type=int,default=5000)
    f.add_argument('--missing',type=float,default=0.)
    f.add_argument('--corrupt',type=float,default=0.)
    e=cmd.add_parser('evaluate')
    e.add_argument('--fitted',type=Path,required=True)
    e.add_argument('--reference',type=Path,required=True)
    e.add_argument('--output',type=Path,required=True)
    e.add_argument('--reference-complete',action='store_true')
    e.add_argument('--top-k',type=int,default=10)
    s=cmd.add_parser('smoke')
    s.add_argument('--output',type=Path,required=True)
    s.add_argument('--seed',type=int,default=42)
    s.add_argument('--gcev',action='store_true',help='Test GCEV using synthetic-only data')
    a=p.parse_args(argv)
    if a.command=='generate-smoke':
        print(generate_smoke(a.output,a.seed,a.cells,a.genes))
    elif a.command=='fit':
        v=run_fit(a.expression,a.output,tf_list=a.tf_list,prior=a.prior,seed=a.seed,
                  max_genes=a.max_genes,max_parents=a.max_parents,k=a.top_k,
                  alpha=a.alpha,graph_query_budget=a.query_budget,missing=a.missing,
                  corrupt=a.corrupt,allow_all_genes=a.allow_all_genes,
                  gcev=a.gcev,prior_secondary=a.prior_secondary)
        print(f"Sealed {v['n_candidate_edges']} candidate edges; reference not accessed")
    elif a.command=='evaluate':
        print(evaluate(a.fitted,a.reference,a.output,
                       reference_complete=a.reference_complete,k=a.top_k).to_string(index=False))
    elif a.command=='smoke':
        data=a.output/'data'
        generate_smoke(data,seed=a.seed)
        run_fit(data/'expression.csv',a.output/'fit',tf_list=data/'tfs.csv',
                prior=data/'TOY_ONLY_random_prior.csv',seed=a.seed,max_genes=8,k=7,
                graph_query_budget=1500,gcev=a.gcev)
        metrics=evaluate(a.output/'fit',data/'reference_complete.csv',
                         a.output/'evaluation',reference_complete=True,k=7)
        print(metrics.to_string(index=False))
        print('SMOKE TEST ONLY — not a biological or confirmatory result')


if __name__=='__main__':
    main()
