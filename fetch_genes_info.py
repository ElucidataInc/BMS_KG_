"""
fetch_genes_info.py
───────────────────
Fetches all Cricetulus griseus (CHO) genes from NCBI and writes genes_info.json.

Usage:
    pip install requests
    python fetch_genes_info.py

Output:
    genes_info.json  ← place this in the same folder as your HTML on GitHub Pages

Runtime: ~60–90 minutes for ~17k genes (NCBI rate limit = 3 req/sec without API key,
         10 req/sec with key). Add your key below to go faster.

NCBI API key (free): https://www.ncbi.nlm.nih.gov/account/
"""

import requests, json, time, sys, os
from xml.etree import ElementTree as ET

# ── CONFIG ────────────────────────────────────────────────────────────────────
NCBI_API_KEY  = ""          # optional but recommended — speeds up 3x
OUTPUT_FILE   = "genes_info.json"
ORGANISM_TAXID = "10029"    # Cricetulus griseus (CHO)
BATCH_SIZE    = 200         # genes per efetch call
RATE_LIMIT    = 0.34        # seconds between requests (3/sec without key, 0.11 with key)
# ─────────────────────────────────────────────────────────────────────────────

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def build_params(extra=None):
    p = {"retmode": "json"}
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    if extra:
        p.update(extra)
    return p

def get(url, params, retries=4):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 429:
                wait = 5 * (attempt + 1)
                print(f"  Rate limited — waiting {wait}s…")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            print(f"  Request error (attempt {attempt+1}): {e}")
            time.sleep(3)
    raise RuntimeError(f"Failed after {retries} attempts: {url}")

# ── STEP 1: Get all gene IDs ───────────────────────────────────────────────
def fetch_all_ids():
    print("Step 1: Fetching all gene IDs for Cricetulus griseus…")
    params = build_params({
        "db": "gene",
        "term": f"txid{ORGANISM_TAXID}[Organism] AND alive[prop]",
        "retmax": 0,
    })
    r = get(f"{BASE}/esearch.fcgi", params)
    total = int(r.json()["esearchresult"]["count"])
    print(f"  Total genes found: {total}")

    all_ids = []
    for start in range(0, total, 500):
        params2 = build_params({
            "db": "gene",
            "term": f"txid{ORGANISM_TAXID}[Organism] AND alive[prop]",
            "retmax": 500,
            "retstart": start,
        })
        r2 = get(f"{BASE}/esearch.fcgi", params2)
        ids = r2.json()["esearchresult"]["idlist"]
        all_ids.extend(ids)
        print(f"  Fetched IDs {start+1}–{min(start+500, total)} / {total}", end="\r")
        time.sleep(RATE_LIMIT)

    print(f"\n  Total IDs collected: {len(all_ids)}")
    return all_ids

# ── STEP 2: Batch fetch gene details via esummary ─────────────────────────
def fetch_batch_summary(id_batch):
    params = build_params({
        "db": "gene",
        "id": ",".join(id_batch),
        "retmode": "json",
    })
    r = get(f"{BASE}/esummary.fcgi", params)
    return r.json().get("result", {})

# ── STEP 3: Parse esummary result into our schema ─────────────────────────
def parse_gene(uid, data):
    name      = data.get("name", "")
    full_name = data.get("description", "")
    org       = data.get("organism", {}).get("scientificname", "Cricetulus griseus")
    gene_type = data.get("type", "")
    chrom     = data.get("chromosome", "")
    map_loc   = data.get("maplocation", "")
    summary   = data.get("summary", "")
    aliases_raw = data.get("otheraliases", "")
    aliases   = [a.strip() for a in aliases_raw.split(",") if a.strip()] if aliases_raw else []

    # Human ortholog from linkstoortholog (not always in esummary — use symbol heuristic)
    # Will be enriched via elink in optional step
    human_orth = None

    # GO terms from genomicinfo / gene2refseq — not in esummary
    # We populate these from Gene Ontology annotation if available
    mol_fn = ""
    bio_proc = ""
    pathways = ""

    # Chromosome formatting
    chr_display = f"Chr {chrom}" if chrom and not chrom.lower().startswith("chr") else chrom

    return {
        "ncbi_id":            uid,
        "symbol":             name,
        "full_name":          full_name,
        "organism":           org,
        "gene_type":          gene_type,
        "chromosome":         chr_display,
        "map_location":       map_loc or None,
        "aliases":            aliases,
        "summary":            summary or "No summary available.",
        "molecular_function": mol_fn or "—",
        "biological_process": bio_proc or "—",
        "pathways":           pathways or "—",
        "human_ortholog":     human_orth,
        "ncbi_url":           f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
    }

# ── STEP 4 (optional): Enrich GO terms via gene2go ────────────────────────
def enrich_go_terms(genes_dict):
    """
    Downloads gene2go for taxid 10029 and adds GO term annotations.
    File size ~200MB uncompressed; skipped if you set SKIP_GO=True below.
    """
    SKIP_GO = False  # set True to skip GO enrichment
    if SKIP_GO:
        print("Skipping GO term enrichment (SKIP_GO=True)")
        return

    go_url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    gz_path = "gene2go.gz"
    txt_path = "gene2go.txt"

    if not os.path.exists(txt_path):
        print("Step 4: Downloading gene2go annotations (~200MB)…")
        import gzip, shutil
        r = requests.get(go_url, stream=True, timeout=120)
        with open(gz_path, "wb") as f:
            for chunk in r.iter_content(65536):
                f.write(chunk)
        with gzip.open(gz_path, "rb") as f_in, open(txt_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)
        print("  Downloaded and extracted gene2go.txt")

    print("  Parsing GO annotations for taxid 10029…")
    mf, bp, cc = {}, {}, {}
    with open(txt_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            tax, gid, go_id, ev, go_name, ns, _ = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
            if tax != ORGANISM_TAXID:
                continue
            if ns == "Function":
                mf.setdefault(gid, []).append(go_name)
            elif ns == "Process":
                bp.setdefault(gid, []).append(go_name)

    enriched = 0
    for uid, gene in genes_dict.items():
        if uid in mf:
            gene["molecular_function"] = "; ".join(mf[uid][:4])
            enriched += 1
        if uid in bp:
            gene["biological_process"] = "; ".join(bp[uid][:4])

    print(f"  GO terms added for {enriched} genes")

# ── MAIN ──────────────────────────────────────────────────────────────────
def main():
    # Resume from partial output if it exists
    genes = {}
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, encoding="utf-8") as f:
            existing = json.load(f)
        if isinstance(existing, list):
            genes = {g["ncbi_id"]: g for g in existing}
        else:
            genes = existing
        print(f"Resuming — {len(genes)} genes already in {OUTPUT_FILE}")

    all_ids = fetch_all_ids()
    remaining = [uid for uid in all_ids if uid not in genes]
    print(f"\nStep 2: Fetching details for {len(remaining)} genes in batches of {BATCH_SIZE}…")

    total = len(remaining)
    done = 0
    for i in range(0, total, BATCH_SIZE):
        batch = remaining[i : i + BATCH_SIZE]
        try:
            result = fetch_batch_summary(batch)
            uids_in_result = result.get("uids", batch)
            for uid in uids_in_result:
                if uid == "uids":
                    continue
                entry = result.get(uid, {})
                if entry:
                    genes[uid] = parse_gene(uid, entry)
            done += len(batch)
            pct = done / total * 100
            print(f"  {done}/{total} ({pct:.1f}%) — {len(genes)} genes parsed", end="\r")
        except Exception as e:
            print(f"\n  Error on batch {i}–{i+BATCH_SIZE}: {e} — skipping")
        time.sleep(RATE_LIMIT if not NCBI_API_KEY else 0.11)

        # Save checkpoint every 2000 genes
        if done % 2000 < BATCH_SIZE:
            with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
                json.dump(list(genes.values()), f, ensure_ascii=False)
            print(f"\n  Checkpoint saved ({len(genes)} genes)")

    print(f"\n\nStep 3: Enriching with GO terms…")
    enrich_go_terms(genes)

    print(f"\nStep 5: Writing {OUTPUT_FILE}…")
    out = list(genes.values())
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"))

    size_mb = os.path.getsize(OUTPUT_FILE) / 1024 / 1024
    print(f"\n✅ Done! {len(out)} genes written to {OUTPUT_FILE} ({size_mb:.1f} MB)")
    print(f"\nNext step: Copy genes_info.json to your GitHub Pages repo root (same folder as index.html)")

if __name__ == "__main__":
    main()
