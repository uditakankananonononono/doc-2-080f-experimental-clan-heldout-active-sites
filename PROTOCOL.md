# DOC-2-080F: experimental-evidence-only, clan-held-out, natural-prevalence follow-up to DOC-2-080 (frozen protocol, lock-1)

Written 2026-10-09 IST and committed BEFORE any data was downloaded, built, embedded or scored in this repo. "080F" names it as the follow-up to DOC-2-080 (ledger numbers 081+ are other topics); it is a NEW preregistered project, not an amendment to DOC-2-080.

## Why this exists
The independent gate on DOC-2-080 passed only a scoped claim (family-held-out, sampled prevalence ~10%, mostly similarity/rule-propagated labels) and required, before any wider claim, four elements: (1) experimental-evidence-only (ECO:0000269) labels, (2) clan-level holdout, (3) natural-prevalence evaluation, (4) an amino-acid-identity-only baseline. This project is exactly those four, nothing else. DOC-2-080 code is reused with the changes below. Prior art: residue-level functional signal in protein language models is known; no novelty claim.

## Data (frozen)
- UniProtKB reviewed, length 100-400, single Pfam id, annotated active site (acquire_act.py; same query as DOC-2-080; release and md5 in DATA_HASHES.tsv). Pfam clan map: Pfam-A.clans.tsv from the Pfam FTP current release at fetch time (acquire_clans.py; md5 recorded).
- Evidence rule: an ACT_SITE feature is experimental iff its own /evidence string contains ECO:0000269. A protein qualifies if it has >= 1 experimental site and at most 5 sites in total, no X/B/Z/U/O. Positives = experimental sites only. Non-experimental sites in a qualifying protein are EXCLUDED from scoring (neither positive nor negative). A feature carrying both ECO:0000269 and ECO:0000250/0000255 counts as experimental (limit: some of these may still be partly propagated).
- Holdout unit = Pfam clan; a family with no clan is its own unit. Up to 10 proteins per unit (seed 81); all units kept (build_080f.py).
- Scored residues = ALL residues of every kept protein (natural prevalence, not subsampled), except excluded non-experimental sites.
- Expectation (a cheap count of the DOC-2-080 download, taken before lock-1, no labels scored): about 580 qualifying proteins in about 290 families before the per-unit cap, about 890 experimental positives, prevalence on the order of 0.5%. Realised counts are reported against this.

## Features
- E: ESM-2 35M (esm2_t12_35M_UR50D) last-layer residue embedding (480), full-length forward pass, CPU, stored float16 (embed_all.py).
- B0 (new, identity-only): amino-acid identity of the residue alone (20 dims). B2: amino-acid identity in a +-3 window (140 dims). No relative-position baseline.

## Method
- Standardize then L2 logistic regression (C = 0.1, max_iter 3000). Training uses all positives plus 15 random negatives per training protein (seed 12345); TEST scoring uses every residue of held-out proteins. Units sorted, assigned round-robin to 5 folds (no unit spans folds); pooled out-of-fold scores.
- Metrics: AUPRC and AUROC pooled at natural prevalence. CI: unit-cluster bootstrap, 2,000 resamples, seed 12345, percentile 95%.

## Gates (analysis_080f.py, mechanical)
- Insufficient data: fewer than 100 units or fewer than 400 positives -> label INSUFFICIENT-DATA, no scoring.
- G1 (control): E trained on within-protein permuted labels has OOF AUPRC < prevalence + 0.01, else INVALID.
- G2: AUPRC(E) - AUPRC(B2) >= +0.05 and CI lower bound > +0.01.
- G2b: AUPRC(E) - AUPRC(B0) >= +0.05 and CI lower bound > +0.01.
- Labels: INVALID; EMBEDDING-BEATS-LOCAL-CONTEXT-AT-NATURAL-PREVALENCE (G2 and G2b pass); BEATS-IDENTITY-ONLY-NOT-WINDOW (G2b pass, G2 fail); HONEST NEGATIVE (otherwise).
- Reported only: E AUPRC CI, E AUPRC / prevalence, share of proteins whose top-5 scored residues contain a positive.

## Limits stated up front
- Clan holdout is stronger than family holdout but does not remove remote homology across clans or between unclassified families; no sequence-identity clustering is done (MMseqs not used here).
- Experimental labels are few and biased toward well-studied enzymes; results may not transfer. UniProt evidence codes are curator-assigned, not independently verified.
- Thresholds (+0.05 AUPRC) are absolute and set before seeing data; natural-prevalence AUPRC is small in scale, so a pass requires a clear absolute lift.
- One model size, last layer, linear probe, single seed, CPU only, training negatives subsampled. B2 is a local-window baseline, not an alignment or conservation baseline, so no claim against alignment tools.
- Narrow target: annotated catalytic residues in enzymes. No claim about binding sites or other motifs.
- Single run; crash fixes are dated AMENDMENT-N.md files committed before outcomes exist. No simulated data in results.
