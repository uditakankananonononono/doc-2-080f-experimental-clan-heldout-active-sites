"""DOC-2-080F: fetch Pfam-A.clans.tsv (current release at fetch time). Writes clans_raw.tsv and appends md5 to DATA_HASHES.tsv."""
import gzip, hashlib, requests
r = requests.get("https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz", timeout=300); r.raise_for_status()
raw = gzip.decompress(r.content).decode(); open("clans_raw.tsv", "w").write(raw)
open("DATA_HASHES.tsv", "a").write(f"clans_raw.tsv\t{hashlib.md5(raw.encode()).hexdigest()}\t{r.headers.get('Last-Modified','NA')}\t{len(raw.splitlines())}\n")
print("CLANS_DONE", len(raw.splitlines()))
