"""DOC-2-080F frozen analysis (PROTOCOL.md lock-1). Run once. Needs proteins.tsv, res_emb.npy."""
import json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import average_precision_score as ap, roc_auc_score as auc
rng = np.random.default_rng(12345)
P = pd.read_csv("proteins.tsv", sep="\t", keep_default_na=False); E = np.load("res_emb.npy")
AA = "ACDEFGHIKLMNPQRSTVWY"; ix = {a: i for i, a in enumerate(AA)}
pid, pos, y, unit, excl = [], [], [], [], []
off = 0; offs = []
for r in P.itertuples():
    ep = {int(v) for v in r.exp_pos.split(",") if v != ""}; nx = {int(v) for v in r.nonexp_pos.split(",") if v != ""}
    for i in range(r.length):
        if i in nx: continue
        pid.append(r.pid); pos.append(i); y.append(int(i in ep)); unit.append(r.unit); offs.append(off + i)
    off += r.length
pid = np.array(pid); pos = np.array(pos); y = np.array(y); unit = np.array(unit); offs = np.array(offs); n = len(y)
E = E[offs].astype(np.float32); seqs = P.set_index("pid").sequence.to_dict()
def window(w):
    X = np.zeros((n, (2 * w + 1) * 20), dtype=np.float32)
    for k, (p, q) in enumerate(zip(pid, pos)):
        s = seqs[p]
        for o in range(-w, w + 1):
            j = q + o
            if 0 <= j < len(s): X[k, (o + w) * 20 + ix[s[j]]] = 1
    return X
B0 = window(0); B2 = window(3)
units = sorted(set(unit)); fold_of = {u: i % 5 for i, u in enumerate(units)}; fid = np.array([fold_of[u] for u in unit])
info = dict(n_units=len(units), n_proteins=len(P), n_residues=int(n), n_positives=int(y.sum()), prevalence=float(y.mean()))
if len(units) < 100 or y.sum() < 400:
    info["LABEL"] = "INSUFFICIENT-DATA"; print("RESULT_JSON", json.dumps(info)); open("results.json", "w").write(json.dumps(info, indent=1)); raise SystemExit
def oof(X, yy):
    s = np.zeros(n)
    for k in range(5):
        a = np.where(fid != k)[0]; neg = a[yy[a] == 0]; keep = np.r_[a[yy[a] == 1], rng.choice(neg, min(len(neg), 15 * len(set(pid[a]))), replace=False)]
        m = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=3000)).fit(X[keep], yy[keep]); b = fid == k; s[b] = m.decision_function(X[b])
    return s
S = {"E": oof(E, y), "B0": oof(B0, y), "B2": oof(B2, y)}
yperm = y.copy()
for p in set(pid): m = pid == p; yperm[m] = rng.permutation(y[m])
sperm = oof(E, yperm); prev = float(y.mean()); res = dict(info, AUPRC={k: float(ap(y, v)) for k, v in S.items()}, AUROC={k: float(auc(y, v)) for k, v in S.items()})
res["G1"] = dict(perm_AUPRC_E=float(ap(yperm, sperm)), prevalence=prev, pass_=bool(ap(yperm, sperm) < prev + 0.01))
members = {u: np.where(unit == u)[0] for u in units}; bs = {"g2": [], "g2b": [], "e": []}
for _ in range(2000):
    idx = np.concatenate([members[u] for u in rng.choice(units, len(units))]); yy = y[idx]
    if yy.sum() == 0: continue
    a, b, c = ap(yy, S["E"][idx]), ap(yy, S["B2"][idx]), ap(yy, S["B0"][idx]); bs["g2"].append(a - b); bs["g2b"].append(a - c); bs["e"].append(a)
ci = lambda k: [float(x) for x in np.percentile(bs[k], [2.5, 97.5])]
d2 = res["AUPRC"]["E"] - res["AUPRC"]["B2"]; d2b = res["AUPRC"]["E"] - res["AUPRC"]["B0"]
res["G2"] = dict(delta_AUPRC_E_minus_B2=d2, ci=ci("g2"), pass_=bool(d2 >= 0.05 and ci("g2")[0] > 0.01))
res["G2b"] = dict(delta_AUPRC_E_minus_B0=d2b, ci=ci("g2b"), pass_=bool(d2b >= 0.05 and ci("g2b")[0] > 0.01))
top = []
for p in sorted(set(pid)):
    m = np.where(pid == p)[0]
    if y[m].sum() == 0: continue
    top.append(int(y[m][np.argsort(-S["E"][m])[:5]].sum() > 0))
res["reported_only"] = dict(E_AUPRC_ci=ci("e"), E_top5_per_protein_hit_rate=float(np.mean(top)), AUPRC_x_prevalence_E=res["AUPRC"]["E"] / prev)
res["LABEL"] = "INVALID" if not res["G1"]["pass_"] else ("EMBEDDING-BEATS-LOCAL-CONTEXT-AT-NATURAL-PREVALENCE" if res["G2"]["pass_"] and res["G2b"]["pass_"] else ("BEATS-IDENTITY-ONLY-NOT-WINDOW" if res["G2b"]["pass_"] else "HONEST NEGATIVE"))
print("RESULT_JSON", json.dumps(res, default=float)); open("results.json", "w").write(json.dumps(res, default=float, indent=1))
