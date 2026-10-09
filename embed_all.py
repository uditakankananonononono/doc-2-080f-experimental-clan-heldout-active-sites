"""DOC-2-080F: ESM-2 35M last-layer embeddings for ALL residues of every protein in proteins.tsv -> res_emb.npy (float16, row = cumulative offset + position)."""
import numpy as np, pandas as pd, torch, time
from transformers import AutoTokenizer, EsmModel
torch.set_num_threads(2); P = pd.read_csv("proteins.tsv", sep="\t")
tok = AutoTokenizer.from_pretrained("facebook/esm2_t12_35M_UR50D"); m = EsmModel.from_pretrained("facebook/esm2_t12_35M_UR50D").eval()
emb = np.zeros((int(P.length.sum()), m.config.hidden_size), dtype=np.float16); off = 0; t0 = time.time()
with torch.no_grad():
    for k, s in enumerate(P.sequence):
        h = m(**tok(s, return_tensors="pt")).last_hidden_state[0, 1:len(s) + 1]; emb[off:off + len(s)] = h.numpy().astype(np.float16); off += len(s)
        if k % 50 == 0: print(k, len(P), round(time.time() - t0), flush=True)
np.save("res_emb.npy", emb); print("EMBED_DONE", emb.shape)
