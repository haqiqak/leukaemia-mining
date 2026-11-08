import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, recall_score
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import RobustScaler, StandardScaler
from sklearn.svm import SVC

CLASSES = ['NEG', 'BCR/ABL', 'ALL1/AF4', 'E2A/PBX1']


def load_data(data_dir='../data'):
    """Patients x genes after the label-free first filter, and the class labels."""
    X = pd.read_csv(f'{data_dir}/expression.csv', index_col='patient', dtype={'patient': str})
    y = pd.read_csv(f'{data_dir}/labels.csv', index_col='patient', dtype={'patient': str})['class']
    kept = pd.read_csv(f'{data_dir}/genes_kept.csv')['probe']
    return X[kept], y


def anova_genes(X, y, p=0.01):
    """Genes whose mean differs between classes (one-way ANOVA, p < 0.01)."""
    pvals = f_classif(X, y)[1]
    return X.columns[pvals < p]


def rf_top_genes(X, y, n=30, seed=0):
    """The n genes with the highest random-forest importance."""
    rf = RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=-1).fit(X, y)
    return X.columns[np.argsort(rf.feature_importances_)[::-1][:n]]


def select(X, y, seed=0):
    """Gene sets for both strategies: s1 = all ANOVA genes, s2 = random-forest top 30 of those."""
    s1 = anova_genes(X, y)
    s2 = rf_top_genes(X[s1], y, seed=seed)
    return {'s1': s1, 's2': s2}


def make_models(seed=0):
    models = {}
    for m in (5, 15):
        models[f'rf_mtry{m}'] = RandomForestClassifier(n_estimators=500, max_features=m,
                                                       random_state=seed, n_jobs=-1)
    for C in (1, 100):
        for gamma in (0.01, 0.001):
            models[f'svm_C{C}_g{gamma}'] = make_pipeline(StandardScaler(), SVC(C=C, gamma=gamma))
    for k in (3, 5, 7):
        models[f'knn{k}_z'] = make_pipeline(StandardScaler(), KNeighborsClassifier(k))
        models[f'knn{k}_iqr'] = make_pipeline(RobustScaler(), KNeighborsClassifier(k))
    return models


def evaluate(X, y, n_rounds=100, test_size=0.3, seed=0, leak=False, models=None, classes=CLASSES):
    """Repeated stratified resampling. Selection is fitted on the training patients of each
    round (leak=False) or once on all patients before resampling (leak=True)."""
    splits = StratifiedShuffleSplit(n_splits=n_rounds, test_size=test_size, random_state=seed)
    genes_all = select(X, y, seed) if leak else None
    models = models or make_models(seed)
    rows, confusion = [], {}
    for r, (tr, te) in enumerate(splits.split(X, y)):
        genes = genes_all if leak else select(X.iloc[tr], y.iloc[tr], seed=r)
        for sel, cols in genes.items():
            for name, model in models.items():
                model.fit(X.iloc[tr][cols], y.iloc[tr])
                pred = model.predict(X.iloc[te][cols])
                true = y.iloc[te]
                row = {'round': r, 'selection': sel, 'model': name,
                       'accuracy': accuracy_score(true, pred),
                       'balanced_accuracy': balanced_accuracy_score(true, pred)}
                row.update(dict(zip(classes, recall_score(true, pred, labels=classes, average=None))))
                rows.append(row)
                cm = confusion_matrix(true, pred, labels=classes)
                confusion[(sel, name)] = confusion.get((sel, name), 0) + cm
    return pd.DataFrame(rows), confusion
