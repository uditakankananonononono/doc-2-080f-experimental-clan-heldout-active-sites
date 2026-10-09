"""DOC-2-080F: build proteins.tsv (frozen rules in PROTOCOL.md). Needs act_raw.tsv, clans_raw.tsv."""
import re, numpy as np, pandas as pd
d = pd.read_csv("act_raw.tsv", sep="\t"); d.columns = ["accession", "pfam", "length", "sequence", "act"]
d["pf"] = d.pfam.fillna("").str.strip(";").str.split(";"); d = d[d.pf.str.len() == 1].copy(); d["family"] = d.pf.str[0]
d = d[~d.sequence.str.contains("[XBZUO]")].copy()
def feats(s):
    ex, nx = [], []
    for f in re.split(r"(?=ACT_SITE \d+;)", str(s)):
        m = re.match(r"ACT_SITE (\d+);", f)
        if m: (ex if "ECO:0000269" in f else nx).append(int(m.group(1)) - 1)
    return sorted(set(ex)), sorted(set(nx) - set(ex))
x = d.act.apply(feats); d["exp_pos"] = x.str[0]; d["nonexp_pos"] = x.str[1]
d = d[(d.exp_pos.apply(len) >= 1) & ((d.exp_pos.apply(len) + d.nonexp_pos.apply(len)) <= 5)].copy()
d = d[d.apply(lambda r: all(0 <= p < len(r.sequence) for p in r.exp_pos + r.nonexp_pos), axis=1)].sort_values("accession")
c = pd.read_csv("clans_raw.tsv", sep="\t", header=None, names=["family", "clan", "clan_name", "fam_name", "desc"]).drop_duplicates("family").set_index("family").clan
d["clan"] = d.family.map(c); d["unit"] = d.clan.fillna(d.family); d["has_clan"] = d.clan.notna()
rng = np.random.default_rng(81); out = []
for u in sorted(d.unit.unique()):
    g = d[d.unit == u]; out.append(g if len(g) <= 10 else g.iloc[np.sort(rng.choice(len(g), 10, replace=False))])
P = pd.concat(out).reset_index(drop=True); P["pid"] = P.index
P["exp_pos"] = P.exp_pos.apply(lambda v: ",".join(map(str, v))); P["nonexp_pos"] = P.nonexp_pos.apply(lambda v: ",".join(map(str, v)))
P[["pid", "accession", "family", "clan", "unit", "has_clan", "length", "sequence", "exp_pos", "nonexp_pos"]].to_csv("proteins.tsv", sep="\t", index=False)
npos = int(P.exp_pos.apply(lambda s: len(s.split(","))).sum())
print("BUILD_DONE units", P.unit.nunique(), "families", P.family.nunique(), "proteins", len(P), "positives", npos, "residues", int(P.length.sum()), "in_clan_proteins", int(P.has_clan.sum()))
