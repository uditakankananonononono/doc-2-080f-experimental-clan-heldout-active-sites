# DOC-2-080F RESULTS - independent gate verdict: SCOPED PASS (claim narrowed as below)

Label (set mechanically by analysis_080f.py): **EMBEDDING-BEATS-LOCAL-CONTEXT-AT-NATURAL-PREVALENCE** (as emitted; the label NAME is superseded by the scoped claim and note below). G1 pass, G2 pass, G2b pass.

## Scoped claim (gate verdict, verbatim)
"Across 136 Pfam clan/family units of reviewed UniProt enzymes with experimentally evidenced (ECO:0000269) active-site annotations, a linear probe on frozen ESM-2 35M residue embeddings ranks those residues above a +-3 AA window (AUPRC 0.187 vs 0.054) and an identity-only baseline (0.037) under 5-fold clan-held-out evaluation with every residue scored; absolute AUPRC is low (about 1 in 5 top-ranked is a positive in the PR sense), and the lift is relative to local baselines."

Label note (gate, verbatim): "The label NAME 'AT-NATURAL-PREVALENCE' overstates: 0.0057 is the prevalence inside the study population (annotated enzymes, at least 1 experimental site, at most 5 sites, up to 10 per unit, length 100-400). It is not proteome-wide prevalence, so no deployment-precision claim. Positives are the experimentally annotated sites; unannotated catalytic residues are scored as negatives."

Nothing beyond this claim is made. Where the table or text below says "natural prevalence", read it as prevalence inside that study population, not proteome-wide.

## Candidates for a NEW preregistered unit (not amendments to 080F)
Items the gate listed as uncontrolled here: (1) sequence-identity-clustered holdout, (2) a conservation baseline, (3) a pretraining-overlap check, (4) multi-seed runs. None was run in 080F. Prior art: residue-level functional signal in protein language models is known; this is a replication-style measurement with a 35M model, not a novelty claim.

| item | value |
|---|---|
| data | UniProt 2026_03, 65,746 reviewed entries with active-site annotation; Pfam-A.clans.tsv (Last-Modified 22 Jan 2026); md5s in DATA_HASHES.tsv and input_md5.txt |
| realised set | 136 holdout units (Pfam clans, or the family itself when it has no clan), 238 families, 396 proteins (333 in a clan, 63 not), 585 experimental (ECO:0000269) positives, 102,533 scored residues (all residues of kept proteins minus 84 excluded non-experimental sites), natural prevalence 0.00571 |
| expectation stated at lock-1 | about 583 proteins / 290 families / 892 positives before the per-unit cap of 10. The cap, sequence-letter filter and per-protein rules cut this to the realised numbers; the insufficiency floor (100 units, 400 positives) was passed |
| G1 control (within-protein permuted labels, E OOF AUPRC) | 0.00636 vs prevalence 0.00571; threshold < 0.01571; pass |
| AUPRC (pooled OOF, clan-held-out, all residues) | E (ESM-2 35M) **0.187**, B2 (+-3 window) 0.0538, B0 (identity only) 0.0371 |
| AUROC | E 0.940, B2 0.903, B0 0.895 |
| G2: AUPRC(E) minus AUPRC(B2) | **+0.133**, 95% unit-bootstrap CI [0.094, 0.179]; needed >= 0.05 and CI lb > 0.01: pass |
| G2b: AUPRC(E) minus AUPRC(B0) | **+0.150**, CI [0.110, 0.196]: pass |
| reported only | E AUPRC CI [0.147, 0.234]; E AUPRC is 32.7x prevalence; top-5 hit rate per protein 0.609 |

## What this shows and does not show
- Under clan-held-out folds, with experimental-only (ECO:0000269) positives and all residues scored at natural prevalence (about 0.6%), a linear probe on frozen ESM-2 35M residue embeddings ranks annotated active-site residues well above both a +-3 amino-acid window and an amino-acid-identity-only baseline. The absolute E AUPRC (0.19) is low in practical terms: most top-ranked residues are still not annotated active sites.
- This addresses the four elements the DOC-2-080 gate required (experimental-only labels, clan holdout, natural prevalence, identity-only baseline). It is a separate project and does not change the DOC-2-080 scoped wording.
- Clan holdout does not remove remote homology across clans, and 63 proteins (families without a clan) are held out only at family level. No sequence-identity clustering was done.
- Experimental labels are few and biased toward well-studied enzymes; 585 positives in 136 units. Evidence codes are curator-assigned. A feature with both ECO:0000269 and propagated evidence counts as experimental.
- Baselines are local amino-acid features only: no alignment, conservation or learned-sequence baseline. No claim against alignment-based tools.
- One model size, last layer, linear probe, single seed, CPU only; training uses 15 negatives per training protein, test uses all residues.
- Narrow target: annotated catalytic residues in enzymes. No claim about binding sites or other motifs. The result may include residues that are catalytic but unannotated and so scored as negatives.

## Disclosures
- Run once. All six locked files match the lock-1 tag (tag_tree_check.txt, tag commit d4255c87, tree 630a860c). No amendments.
- The first acquire_act.py attempt was killed by a tool timeout before completing and was rerun from scratch; the data used is the complete rerun (act_raw.tsv md5 ca3c4cc9). acquire_act.py as locked overwrites DATA_HASHES.tsv, so the clans_raw.tsv md5 line (saved in DATA_HASHES_clans.tsv) was appended by hand.
- run_log.txt records UTC times, versions and the command. The exit code of analysis_080f.py was not captured separately; results.json and the RESULT_JSON line in analysis.log are complete.
- Counts differ from the lock-1 expectation (see table); the protocol's rules were applied as written.
- The analysis script was smoke-tested only on synthetic random data in a scratch directory before the run; no record is kept and no result relies on it.
- Large files (act_raw.tsv, proteins.tsv, res_emb.npy as float16) are in the Drive folder; res_emb.npy may be split into parts if over the Drive per-file cap.
