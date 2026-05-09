#!/usr/bin/env python3
"""
Step 042: SPARC Database Ingestion
Downloads rotation curve data from the official SPARC archive.
Source: Lelli, McGaugh & Schombert 2016, AJ, 152, 157
URL: http://astroweb.cwru.edu/SPARC/
"""

import sys, hashlib, json, time, urllib.request, urllib.error
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.utils.logger import TEPLogger, set_step_logger, print_status

SPARC_URL = "http://astroweb.cwru.edu/SPARC"
DATA_DIR = PROJECT_ROOT / "data" / "raw" / "sparc"

# 175 SPARC galaxies from Lelli et al. 2016, Table 1
GALAXIES = [
    "DDO064","DDO154","DDO161","DDO168","DDO170",
    "ESO116-G012","ESO563-G021",
    "F563-1","F563-V1","F563-V2","F565-V1","F565-V2","F567-2",
    "F568-1","F568-3","F568-V1","F571-8","F571-V1","F574-1","F574-2",
    "F579-V1","F583-1","F583-4",
    "IC0257","IC0420","IC0750","IC1029","IC2233","IC2574",
    "KK98-250","KK98-251",
    "NGC0024","NGC0055","NGC0247","NGC0253","NGC0289",
    "NGC0300","NGC0428","NGC0598","NGC0628","NGC0801",
    "NGC0891","NGC0925","NGC1003","NGC1090","NGC1560",
    "NGC2366","NGC2403","NGC2683","NGC2841","NGC2903",
    "NGC2915","NGC2955","NGC2976","NGC2998","NGC3031",
    "NGC3109","NGC3198","NGC3319","NGC3351","NGC3521",
    "NGC3621","NGC3726","NGC3741","NGC3769","NGC3877",
    "NGC3893","NGC3917","NGC3949","NGC3953","NGC3972",
    "NGC3992","NGC4010","NGC4013","NGC4051","NGC4062",
    "NGC4085","NGC4088","NGC4100","NGC4138","NGC4157",
    "NGC4183","NGC4190","NGC4214","NGC4217","NGC4220",
    "NGC4244","NGC4258","NGC4302","NGC4389","NGC4414",
    "NGC4449","NGC4455","NGC4490","NGC4559","NGC4565",
    "NGC4605","NGC4631","NGC4656","NGC4736","NGC5005",
    "NGC5023","NGC5033","NGC5055","NGC5204","NGC5371",
    "NGC5474","NGC5585","NGC5907","NGC6015","NGC6195",
    "NGC6503","NGC6674","NGC6689","NGC6946","NGC7331",
    "NGC7793","NGC7814",
    "PGC51017",
    "UGC00128","UGC00191","UGC00634","UGC00731","UGC00849",
    "UGC01230","UGC01281","UGC02023","UGC02259","UGC02455",
    "UGC02885","UGC02916","UGC02953","UGC03205","UGC03546",
    "UGC03580","UGC04278","UGC04305","UGC04325","UGC04499",
    "UGC05253","UGC05414","UGC05716","UGC05721","UGC05750",
    "UGC05764","UGC05829","UGC05918","UGC05986","UGC06399",
    "UGC06446","UGC06614","UGC06628","UGC06786","UGC06787",
    "UGC06917","UGC06923","UGC06930","UGC06973","UGC06983",
    "UGC07089","UGC07125","UGC07151","UGC07261","UGC07323",
    "UGC07399","UGC07524","UGC07559","UGC07577","UGC07603",
    "UGC07608","UGC07690","UGC07866","UGC07950","UGC08091",
    "UGC08286","UGC08320","UGC08490","UGC08550","UGC08621",
    "UGC09037","UGC09133","UGC09215","UGC09992","UGC10310",
    "UGC12082","UGC12506","UGC12632","UGC12709","UGC12732",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while c := f.read(65536): h.update(c)
    return h.hexdigest()


def download(url, dest, logger):
    if dest.exists():
        logger.info(f"  Exists: {dest.name}")
        return True
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "TEP-UCD/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = r.read()
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        logger.info(f"  Saved: {dest.name} ({len(data)} B)")
        return True
    except Exception as e:
        logger.warning(f"  Failed: {dest.name} — {e}")
        return False


def main():
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = TEPLogger("step_042", str(log_dir / "step_042_sparc_data_ingestion.log"))
    set_step_logger(logger)

    print_status("SPARC Data Ingestion — downloading from astroweb.cwru.edu/SPARC/", "TITLE")
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    manifest = {"source": "Lelli et al. 2016, AJ, 152, 157", "base_url": SPARC_URL,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "files": {}}

    # Master table
    table = "SPARC_Lelli2016c.mrt"
    if download(f"{SPARC_URL}/{table}", DATA_DIR / table, logger):
        p = DATA_DIR / table
        manifest["files"][table] = {"sha256": sha256(p), "size": p.stat().st_size}

    # Rotation curves
    ok, fail = 0, 0
    for name in GALAXIES:
        fn = f"{name}_rotmod.dat"
        if download(f"{SPARC_URL}/{fn}", DATA_DIR / fn, logger):
            ok += 1
            p = DATA_DIR / fn
            manifest["files"][fn] = {"sha256": sha256(p), "size": p.stat().st_size}
        else:
            fail += 1
        time.sleep(0.05)

    with open(DATA_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    results = {
        "step_id": "step_042",
        "source": "Lelli et al. 2016, AJ, 152, 157",
        "source_url": SPARC_URL,
        "status": "PASS" if ok > 100 else "PARTIAL",
        "downloaded": ok, "failed": fail, "expected": len(GALAXIES),
    }

    logger.info(f"Downloaded {ok}/{len(GALAXIES)} rotation curves ({fail} failed)")
    logger.save_step_results(results, PROJECT_ROOT, "step_042_sparc_data_ingestion")


if __name__ == "__main__":
    main()
