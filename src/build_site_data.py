"""Collect the numbers shown on the project page from the saved results into docs/site/assets/data.js.
Run from the project root: python src/build_site_data.py"""
import json
import re
import shutil

import pandas as pd
from scipy.stats import wilcoxon

R, D, OUT = 'results/', 'data/', 'docs/site/assets/'
C = ['NEG', 'BCR/ABL', 'ALL1/AF4', 'E2A/PBX1']
M = ['accuracy', 'balanced_accuracy']


def by_workflow(df):
    df = df.assign(workflow=df.selection + '_' + df.model)
    return df, df.groupby('workflow')[M + C].mean()


labels = pd.read_csv(D + 'labels.csv', dtype=str)
n_genes = len(open(D + 'expression.csv').readline().split(',')) - 1
nb = json.load(open('notebooks/02_baseline.ipynb'))
text = ''.join(''.join(o.get('text', [])) for c in nb['cells'] if c['cell_type'] == 'code' for o in c['outputs'])
n_anova = int(re.search(r'ANOVA genes on all \d+ patients: (\d+)', text).group(1))

base, b = by_workflow(pd.read_csv(R + 'baseline_scores.csv'))
leak, l = by_workflow(pd.read_csv(R + 'leak_scores.csv'))
sh_in, s_in = by_workflow(pd.read_csv(R + 'shuffle_inside_scores.csv'))
sh_out, s_out = by_workflow(pd.read_csv(R + 'shuffle_outside_scores.csv'))
sd = base.groupby('workflow')['accuracy'].std()

acc_p = lambda d: d.pivot(index='round', columns='workflow', values='accuracy')
bp, lp = acc_p(base), acc_p(leak)
gap = l['accuracy'] - b['accuracy']
best_acc, best_bal = b['accuracy'].idxmax(), b['balanced_accuracy'].idxmax()
top10 = b['accuracy'].sort_values(ascending=False).index[:10]
p_top = [wilcoxon(bp[best_acc], bp[w]).pvalue for w in top10[1:]]
p_leak = [wilcoxon(lp[w], bp[w]).pvalue for w in b.index]
rank = b['accuracy'].rank(ascending=False)
is_rf = b.index.str.contains('_rf_')
fam = b['accuracy'].groupby(b.index.str.split('_').str[1].str.replace(r'\d+', '', regex=True)).mean()

conf = pd.read_csv(R + 'baseline_confusion.csv')
sel, model = best_acc.split('_', 1)
cm = conf[(conf.selection == sel) & (conf.model == model)].set_index('true')[C].loc[C]

within = pd.read_csv(R + 'b1_within_stage_scores.csv').groupby('selection')[M + ['ALL1/AF4', 'NEG']].mean()
stage = pd.read_csv(R + 'stage_prediction_scores.csv', index_col=[0, 1])
stage = {k: stage.loc[k].groupby('selection')[M].mean().loc['s2'] for k in ('five_stages', 'b1_vs_other')}
counts = labels['class'].value_counts()
b1 = labels[(labels.BT == 'B1') & (labels['class'] != 'BCR/ABL')]['class']

data = {
    'data': {'patients': len(labels), 'genes': n_genes, 'first_filter': len(pd.read_csv(D + 'genes_kept.csv')),
             'anova': n_anova, 'classes': counts.to_dict(), 'always_neg': counts.max() / len(labels),
             'stage_by_class': pd.crosstab(labels['BT'], labels['class']).to_dict('index')},
    # figures quoted from the book's chapter (docs/01_chapter7_walkthrough.md), not computed here
    'book': {'first_filter': 3943, 'anova': 746, 'top_genes': 30, 'workflows': 54, 'rounds': 100,
             'best_accuracy': 0.813, 'ensemble_members': 50, 'clusters': 30},
    'workflows': [{'name': w, 'selection': w[:2], 'model': w[3:], 'accuracy': b.loc[w, 'accuracy'],
                   'sd': sd[w], 'balanced': b.loc[w, 'balanced_accuracy'], 'recall': b.loc[w, C].to_dict(),
                   'leak_accuracy': l.loc[w, 'accuracy'], 'leak_balanced': l.loc[w, 'balanced_accuracy']}
                  for w in b.index],
    'summary': {
        'best_acc': {'name': best_acc, **b.loc[best_acc, M + C].to_dict()},
        'best_bal': {'name': best_bal, **b.loc[best_bal, M + C].to_dict()},
        'rf_best': {'name': b[is_rf]['accuracy'].idxmax(), 'accuracy': b[is_rf]['accuracy'].max(),
                    'rank': rank[b[is_rf]['accuracy'].idxmax()]},
        'rf_zero_e2a': {'name': 's1_rf_mtry5', **b.loc['s1_rf_mtry5', M + C].to_dict()},
        'svm_failing': [{'name': w, **b.loc[w, M + C].to_dict()} for w in b.index[(b['accuracy'] < 0.6)]],
        'family_accuracy': fam.to_dict(),
        'top10_p_below_005': int(sum(p < 0.05 for p in p_top)), 'top10_others': len(p_top),
        'leak': {'mean_gap': gap.mean(), 'mean_gap_balanced': (l['balanced_accuracy'] - b['balanced_accuracy']).mean(),
                 'gap_s1': gap[gap.index.str[:2] == 's1'].mean(), 'gap_s2': gap[gap.index.str[:2] == 's2'].mean(),
                 'max_gap': gap.max(), 'max_gap_name': gap.idxmax(),
                 'best_outside': l['accuracy'].max(), 'best_outside_name': l['accuracy'].idxmax(),
                 'significant': int(sum(p < 0.01 for p in p_leak)), 'n': len(p_leak),
                 'higher_outside': int((gap > 0).sum())},
        'shuffle': {'inside_mean': sh_in[M].mean().to_dict(), 'outside_mean': sh_out[M].mean().to_dict(),
                    'inside_best': s_in['accuracy'].max(), 'outside_best': s_out['accuracy'].max(),
                    'inside_best_balanced': s_in['balanced_accuracy'].max(),
                    'outside_best_balanced': s_out['balanced_accuracy'].max(),
                    'rounds': int(sh_in['round'].nunique())},
        'b1_within': {'patients': len(b1), 's1': within.loc['s1'].to_dict(), 's2': within.loc['s2'].to_dict(),
                      'reference': b1.value_counts().max() / len(b1)},
        'stage': {k: v.to_dict() for k, v in stage.items()},
        'rounds': int(base['round'].nunique()), 'holdout': 0.3},
    'confusion': {'name': best_acc, 'classes': C, 'matrix': cm.values.tolist()},
}

with open(OUT + 'data.js', 'w') as f:
    f.write('window.RESULTS = ' + json.dumps(data, default=float, separators=(',', ':')) + ';\n')
shutil.copy(R + 'data_check.png', OUT + 'data_check.png')
print('written', OUT)
