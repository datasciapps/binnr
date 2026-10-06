"""Hartman et al. (2023) BINN on the shipped septic-AKI proteomics matrix.

Trains the Python `binn` package (v0.1.1) with its documented defaults
(n_layers=4, tanh, dropout 0.2, Adam lr 1e-4, batch 32, 100 epochs) under
repeated stratified cross-validation, records test AUROC/accuracy, and
exports (a) the layered connectivity matrices and (b) SHAP node importances
so that binnr can fit the same topology and compare explanations.
"""
import json, sys, time, os, warnings
import numpy as np, pandas as pd, torch
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn import preprocessing
warnings.filterwarnings("ignore")

from binn import BINN, BINNTrainer, BINNExplainer

OUT = os.environ.get("BINNR_REPLICATION_OUT", "replication-out"); os.makedirs(OUT, exist_ok=True)
data_matrix = pd.read_csv("binn/data/  (clone of InfectionMedicineProteomics/BINN) sample_datamatrix.csv")
design = pd.read_csv("binn/data/  (clone of InfectionMedicineProteomics/BINN) sample_design_matrix.tsv", sep="\t")
N_REP, K, EPOCHS, LR, BATCH = int(sys.argv[1]) if len(sys.argv) > 1 else 5, 5, 100, 1e-4, 32

torch.manual_seed(0)
binn = BINN(data_matrix=data_matrix, network_source="reactome", n_layers=4, dropout=0.2)

# ---- export topology: one edge list per layer, inputs in binn.inputs order
edges = []
for li, mat in enumerate(binn.connectivity_matrices):
    m = mat.values.astype(bool)
    rows, cols = np.where(m)
    for r, c in zip(rows, cols):
        edges.append((li, mat.index[r], mat.columns[c]))
pd.DataFrame(edges, columns=["layer", "from", "to"]).to_csv(f"{OUT}/binn_edges.csv", index=False)
pd.Series(binn.inputs).to_csv(f"{OUT}/binn_inputs.csv", index=False, header=["protein"])
print("layers:", [m.shape for m in binn.connectivity_matrices], "inputs:", len(binn.inputs))

# ---- data in the model's input order (binn's dataloader does the same, then standardises)
X_all = data_matrix.set_index("Protein").loc[binn.inputs].T  # samples x proteins
# binn's own dataloader fills missing intensities with 0 before standardising (37% of cells)
X_all = X_all.loc[design["sample"].values].fillna(0).values.astype(np.float32)
y_all = (design["group"].values == 2).astype(int)  # group 2 = "AKI" per design file (123 vs 74)
print("n =", len(y_all), "positives:", y_all.sum())

def make_loader(X, y, shuffle):
    ds = torch.utils.data.TensorDataset(torch.tensor(X, dtype=torch.float32), torch.tensor(y, dtype=torch.long))
    return torch.utils.data.DataLoader(ds, batch_size=BATCH, shuffle=shuffle)

results, importances, pred_rows = [], [], []
for rep in range(N_REP):
    skf = StratifiedKFold(n_splits=K, shuffle=True, random_state=1000 + rep)
    for fold, (tr, te) in enumerate(skf.split(X_all, y_all)):
        scaler = preprocessing.StandardScaler().fit(X_all[tr])
        Xtr, Xte = scaler.transform(X_all[tr]), scaler.transform(X_all[te])
        torch.manual_seed(rep * 100 + fold)
        model = BINN(data_matrix=data_matrix, network_source="reactome", n_layers=4, dropout=0.2)
        loaders = {"train": make_loader(Xtr, y_all[tr], True), "val": make_loader(Xte, y_all[te], False)}
        t0 = time.time()
        trainer = BINNTrainer(model)
        import contextlib, io
        with contextlib.redirect_stdout(io.StringIO()):
            trainer.fit(dataloaders=loaders, num_epochs=EPOCHS, learning_rate=LR)
        model.eval()
        with torch.no_grad():
            logits = model(torch.tensor(Xte, dtype=torch.float32))
            prob = torch.softmax(logits, 1)[:, 1].numpy()
        auc = roc_auc_score(y_all[te], prob); acc = accuracy_score(y_all[te], prob > 0.5)
        results.append(dict(rep=rep + 1, fold=fold + 1, auroc=auc, accuracy=acc, seconds=time.time() - t0))
        for i, p in zip(te, prob):
            pred_rows.append(dict(rep=rep + 1, fold=fold + 1, sample=design["sample"].values[i], prob=p, y=int(y_all[i])))
        print(f"rep {rep+1} fold {fold+1}: AUROC {auc:.3f} acc {acc:.3f} ({time.time()-t0:.0f}s)", flush=True)
        # SHAP importances on the training data (background = training data), their defaults
        if rep == 0:
            with contextlib.redirect_stdout(io.StringIO()):
                expl = BINNExplainer(model).explain_single({"train": loaders["train"]}, split="train",
                                                           normalization_method="subgraph")
            expl["rep"], expl["fold"] = rep + 1, fold + 1
            importances.append(expl)

pd.DataFrame(results).to_csv(f"{OUT}/binn_cv.csv", index=False)
pd.DataFrame(pred_rows).to_csv(f"{OUT}/binn_pred.csv", index=False)
pd.concat(importances).to_csv(f"{OUT}/binn_shap.csv", index=False)
r = pd.DataFrame(results); print(r[["auroc", "accuracy"]].agg(["mean", "std"]))
