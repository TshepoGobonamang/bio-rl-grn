"""Publication-style descriptive BIO-V8 source-action plot, not accuracy."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.summary)
    frame = frame[(frame['transform'] == 'log1p_original') & (frame.variant == 'primary_mixed')].set_index('representation')
    keys = ['accessibility_only', 'motif', 'global_mean_accessibility_motif',
            'matched_accessibility_motif', 'permuted_accessibility_motif']
    labels = ['Accessibility alone', 'Motif', 'Globally attenuated motif',
              'Matched accessibility × motif', 'Permuted accessibility × motif']
    colors = ['#84949e', '#315f83', '#c08c36', '#087f83', '#8293b4']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.spines.left': False})
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9), gridspec_kw={'width_ratios': [1, 1]})
    for ax, field, title in zip(axes,
            ['changed_vs_expression', 'structurally_varying_targets'],
            ['A  Changed TF choices', 'B  Targets with varying source values']):
        values = frame.loc[keys, field].to_numpy(int)
        ax.barh(range(len(keys)), values, color=colors, height=.62)
        ax.set_yticks(range(len(keys)), labels if ax is axes[0] else [''] * len(keys))
        ax.invert_yaxis()
        ax.set_xlim(0, 72)
        ax.set_xticks([0, 16, 32, 48, 64])
        ax.grid(axis='x', alpha=.17)
        ax.set_axisbelow(True)
        ax.tick_params(axis='y', length=0)
        ax.set_title(title, loc='left', fontweight='bold', pad=16, fontsize=11)
        ax.set_xlabel('Targets (out of 64)')
        for j, value in enumerate(values):
            ax.text(value + 1.4, j, str(value), va='center', fontweight='bold', color='#203748')
    fig.suptitle('Which evidence can affect the graph-edit decision?', x=.025, ha='left',
                 fontsize=16, fontweight='bold', color='#183447')
    fig.text(.025, .90, 'Frozen original-scale fitting scores • primary accessibility processing • boost = 0.3', color='#566878', fontsize=10)
    fig.text(.025, .035, 'Decision influence only: no new predictive or biological accuracy was measured.\nMatched and permuted interactions change 8 choices each, but differ from each other on 8 targets.',
             color='#566878', fontsize=9)
    fig.subplots_adjust(left=.29, right=.98, top=.79, bottom=.23, wspace=.16)
    a.output.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'svg', 'pdf']:
        fig.savefig(a.output / ('F03_BIOV8_Source_Action_Qualification.' + ext), dpi=220,
                    facecolor='white', metadata={'Creator': 'BIO-V8 source-action report'} if ext == 'pdf' else None)
    plt.close(fig)


if __name__ == '__main__':
    main()
