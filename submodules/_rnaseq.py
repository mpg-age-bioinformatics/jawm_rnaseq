import jawm
import os
import re
import pandas as pd
import random
import string
from io import StringIO
from pathlib import Path
import requests
import xml.etree.ElementTree as ET
from typing import Optional, Tuple, List, Tuple, Dict, Union, Any
import gzip
import io
import csv
import warnings


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

AGEPY_IMAGE="mpgagebioinformatics/agepy:9acb5de"



def _bucket_dir_for_gse(gse: str) -> str:
    """
    GEO series directory bucketing:
      GSE127917 -> GSE127nnn
      GSE999    -> GSEnnn
    """
    n = int(gse[3:])
    prefix = str(n)[:-3]  # drop last 3 digits
    return f"GSE{prefix}nnn" if prefix else "GSEnnn"


def _download_family_soft_text(
    gse: str,
    *,
    timeout: int = 60,
    user_agent: str = "geo-soft-parser/1.0 (contact: you@example.com)",
) -> str:
    gse = gse.strip().upper()
    bucket = _bucket_dir_for_gse(gse)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/soft/{gse}_family.soft.gz"

    sess = requests.Session()
    sess.headers.update({"User-Agent": user_agent})
    r = sess.get(url, timeout=timeout)
    r.raise_for_status()

    with gzip.GzipFile(fileobj=io.BytesIO(r.content)) as gz:
        return gz.read().decode("utf-8", errors="replace")

# _PRJ_RE = re.compile(r"\b(PRJ[A-Z]{2}\d+)\b")  # PRJNA..., PRJEB..., PRJDB...

# _PRJ_RE = re.compile(r"(PRJ[A-Z]{2}\d+)")

def parse_bioproject_and_gsms_from_family_soft(soft_text: str) -> Tuple[Optional[str], List[Tuple[str, str]]]:
    """
    Returns:
      (bioproject_or_None, [(GSM, Sample_title), ...])
    """
    bioproject: Optional[str] = None
    samples: List[Tuple[str, str]] = []

    current_gsm: Optional[str] = None
    current_title: str = ""

    for line in soft_text.splitlines():
        # Series-level relation lines are where BioProject often appears
        if line.startswith("!Series_relation"):
            if bioproject is None:
                if "bioproject/" in line:
                    m = line.split("bioproject/")[-1]
                    if m:
                        bioproject = m.strip()

        # Start of a sample block
        if line.startswith("^SAMPLE"):
            # flush previous GSM if any
            if current_gsm is not None:
                samples.append((current_gsm, current_title))

            m = re.search(r"\b(GSM\d+)\b", line)
            current_gsm = m.group(1) if m else None
            current_title = ""
            continue

        # Within a sample block, grab Sample_title (your "N2_1", etc.)
        if current_gsm is not None and line.startswith("!Sample_title"):
            # Format: !Sample_title = N2_1
            parts = line.split("=", 1)
            current_title = parts[1].strip() if len(parts) == 2 else ""

    # flush last sample
    if current_gsm is not None:
        samples.append((current_gsm, current_title))

    return bioproject, samples


def geo_gse_to_bioproject_and_samples(
    gse: str,
    *,
    timeout: int = 60,
    user_agent: str = "geo-soft-parser/1.0 (contact: you@example.com)",
) -> Tuple[Optional[str], List[Tuple[str, str]]]:
    """
    Convenience wrapper: downloads + parses.
    """
    soft_text = _download_family_soft_text(gse, timeout=timeout, user_agent=user_agent)
    parsed=parse_bioproject_and_gsms_from_family_soft(soft_text)
    return parsed

_GSM_RE = re.compile(r"\bGSM\d+\b")
_PRJ_RE = re.compile(r"^PRJ[A-Z]{2}\d+$", re.IGNORECASE)


def sra_runinfo_for_bioproject_find_gsm_columns(
    bioproject: str,
    *,
    timeout: int = 60,
    user_agent: str = "sra-runinfo-fetch/1.0 (contact: you@example.com)",
    retmax: int = 200000,
) -> Tuple[Any, str, List[str]]:
    """
    Fetch SRA RunInfo metadata for a BioProject and identify columns containing GSM accessions.

    Returns:
      (table, run_column_name, gsm_columns)

    - table:
        If pandas is available -> pandas.DataFrame
        Else -> List[Dict[str, str]] rows
    - run_column_name:
        Usually 'Run' in RunInfo CSV
    - gsm_columns:
        Column names where at least one value matches /GSM\\d+/

    Notes:
    - Uses NCBI E-utilities: esearch (db=sra, term=<bioproject>, usehistory=y) -> efetch (rettype=runinfo).
    - RunInfo schema varies by experiment; GSM may appear in columns like 'GEO_Accession', etc.
    """
    bioproject = bioproject.strip()
    if not _PRJ_RE.match(bioproject):
        raise ValueError(f"Expected a BioProject accession like PRJNA592073, got: {bioproject!r}")

    sess = requests.Session()
    sess.headers.update({"User-Agent": user_agent})

    # 1) ESearch SRA with history
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {
        "db": "sra",
        "term": bioproject,
        "usehistory": "y",
        "retmode": "json",
        "retmax": str(retmax),
    }
    r = sess.get(esearch_url, params=esearch_params, timeout=timeout)
    r.raise_for_status()
    j = r.json()

    esr = j.get("esearchresult", {})
    count = int(esr.get("count", "0") or 0)
    webenv = esr.get("webenv")
    query_key = esr.get("querykey")

    if not count or not webenv or not query_key:
        # Could be private/removed, typo, or no SRA runs indexed for that BioProject
        raise RuntimeError(
            f"No SRA records found for {bioproject} (count={count}). "
            "If this is a GEO-only project or private, RunInfo may not be available."
        )

    # 2) EFetch RunInfo CSV using the search history
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    efetch_params = {
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "runinfo",
        "retmode": "text",
    }
    r2 = sess.get(efetch_url, params=efetch_params, timeout=timeout)
    r2.raise_for_status()
    text = r2.text

    # RunInfo sometimes returns an error page or empty content; guard lightly
    if "Error:" in text[:200] or len(text.strip()) == 0:
        raise RuntimeError("EFetch did not return a valid RunInfo CSV payload.")

    # Parse CSV
    reader = csv.DictReader(io.StringIO(text))
    rows: List[Dict[str, str]] = []
    for row in reader:
        # normalize None -> ""
        rows.append({k: (v if v is not None else "") for k, v in row.items()})

    if not rows:
        raise RuntimeError("RunInfo CSV parsed but contained no rows.")

    # Identify RUN column (most commonly 'Run', sometimes 'RUN' though rare)
    fieldnames = list(rows[0].keys())
    run_col = None
    for cand in ("Run", "RUN", "run"):
        if cand in fieldnames:
            run_col = cand
            break
    if run_col is None:
        # fallback: any column name containing 'run' exactly once
        for f in fieldnames:
            if f.lower() == "run":
                run_col = f
                break
    if run_col is None:
        raise RuntimeError(f"Could not find a RUN column in RunInfo fields: {fieldnames[:30]} ...")

    # Find columns containing GSM accessions
    gsm_cols: List[str] = []
    for col in fieldnames:
        # scan until first match for speed
        for row in rows:
            val = row.get(col, "") or ""
            if _GSM_RE.search(val):
                gsm_cols.append(col)
                break

    # Return as pandas DataFrame if available (nice for downstream work)
    try:
        import pandas as pd  # type: ignore
        df = pd.DataFrame(rows)
        return df, run_col, gsm_cols
    except Exception:
        return rows, run_col, gsm_cols


def clean_group(x):
    x_=x[:-1]
    x_=x_.replace("-","_").replace("rep","").replace("Rep","").replace("REP","").strip("_ ").replace(" ","_")
    return x_


def get_gse_groups(gse):

    bioproject, samples = geo_gse_to_bioproject_and_samples(gse)
    table, run_col, gsm_cols = sra_runinfo_for_bioproject_find_gsm_columns(bioproject)
    samples=pd.DataFrame(samples, columns=["GSM","Sample_Name"])
    df=table[[run_col,gsm_cols[0] ]]
    df.columns=["SRR","GSM"]
    df=pd.merge(df,samples,on=["GSM"],how="left")
    df = (
        df
        .groupby("Sample_Name", as_index=False)
        .agg(lambda x: ",".join(x.astype(str)))
    )
    df["group"]=df["Sample_Name"].apply( lambda x: clean_group(x) )
    df=df[["SRR","group"]]
    groups=df.to_csv(sep=";", index=False, header=False ).rstrip("\n").replace("\n", r"\n")
    return groups


def _first_text(root: ET.Element, xpath: str) -> Optional[str]:
    el = root.find(xpath)
    if el is not None and el.text:
        t = el.text.strip()
        return t if t else None
    return None


def srr_to_bioproject_and_organism_eutils(
    srr: str,
    *,
    email: str = "your_email@example.org",
    tool: str = "srr-metadata-script",
    timeout: int = 30,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Given an SRA Run accession (SRR...), return (bioproject_id, organism)
    by querying NCBI SRA via E-utilities.
    """

    # 1) esearch: SRR -> internal SRA UID
    r = requests.get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "sra",
            "term": srr,
            "retmode": "xml",
            "tool": tool,
            "email": email,
        },
        timeout=timeout,
    )
    r.raise_for_status()
    esearch_root = ET.fromstring(r.text)
    uid = _first_text(esearch_root, ".//IdList/Id")
    if not uid:
        raise ValueError(f"No SRA record found for run {srr}")

    # 2) efetch docsum: UID -> metadata XML
    r = requests.get(
        f"{EUTILS}/efetch.fcgi",
        params={
            "db": "sra",
            "id": uid,
            "rettype": "docsum",
            "retmode": "xml",
            "tool": tool,
            "email": email,
        },
        timeout=timeout,
    )
    r.raise_for_status()
    root = ET.fromstring(r.text)

    # Common direct fields (these names can differ by record)
    bioproject = None
    for name in ("Bioproject", "BioProject", "bioproject"):
        el = root.find(f".//Item[@Name='{name}']")
        if el is not None and el.text and el.text.strip():
            bioproject = el.text.strip()
            break

    organism = None
    for name in ("Organism", "ScientificName"):
        el = root.find(f".//Item[@Name='{name}']")
        if el is not None and el.text and el.text.strip():
            organism = el.text.strip()
            break

    # Fallback: parse ExpXml / Runs XML blob if needed
    if not bioproject or not organism:
        expxml_el = root.find(".//Item[@Name='ExpXml']")
        if expxml_el is not None and expxml_el.text:
            blob = expxml_el.text

            if not bioproject:
                # Look for PRJNA/PRJEB/PRJDB patterns
                m = re.search(r"\bPRJ(?:NA|EB|DB)\d+\b", blob)
                if m:
                    bioproject = m.group(0)

            if not organism:
                # Often present like: <Organism ScientificName="Homo sapiens" .../>
                m = re.search(r'ScientificName="([^"]+)"', blob)
                if m:
                    organism = m.group(1)

    return bioproject, organism.lower().replace(" ", "_")


def ccg_samples( groups, source_folder ):
    """
    Processing sample tables from the Cologne Center for Genomics
    """
    df = pd.read_csv(StringIO(s), sep=";", header=None, names=["read1", "read2", "group"])
    for c in df.columns.tolist() :
        df[c]=df[c].apply(lambda x: x.strip(" ") )
    folder = Path(source_folder)
    raw_files=list(folder.glob("*.fastq.gz")) + list(folder.glob("*.fq.gz"))
    read_1=[ s for s in raw_files if "_R1_" in s ]
    read_2=[ s for s in raw_files if "_R2_" in s ]
    read_1.sort()
    read_2.sort()
    df["read2"]=df["read1"].apply(lambda x: ",".join( [ s for s in read_2 if f"_{x}_" in s ]    ))
    df["read1"]=df["read1"].apply(lambda x: ",".join( [ s for s in read_1 if f"_{x}_" in s ]    ))
    return df.to_csv(sep=";", index=False, header=False).rstrip("\n")
    
def load_table(file_path):
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    if ext in [".xls", ".xlsx"]:
        # Excel file
        df=pd.read_excel(file_path)
    
    elif ext in [".tsv", ".csv"]:
        # TSV file (tab-separated)
        print("Reading file as tab separated file")
        try:
            df=pd.read_csv( file_path, sep="\\t" )
        except:
            raise ValueError(f"Could not read tab separated file")

    else:
        raise ValueError(f"Unsupported file type: {ext}")

    if len( df.columns.tolist() ) < 2 :
        raise ValueError(f"Imported table had less than 2 columns.")

    return df


def process_groups( groups ) :

    if os.path.isfile(groups) :
        groups=load_table(groups )
    
    else:
        groups = groups.rstrip(';')
        groups = groups.replace("\\n", "\n")
        if groups.endswith("\n"):
            groups = groups[:-1]
        groups = "\n".join(line.strip() for line in groups.split("\n"))
        groups=[ s.split(";") for s in groups.split("\n") ]
        groups=pd.DataFrame(groups)

    value=groups.iloc[0, 0]

    groups=groups.to_csv(sep=";", index=False, header=False ).rstrip("\n").replace("\n", r"\n")

    if value.startswith("SRR") and not value.endswith((".fastq", ".fq", ".fastq.gz")):
        bioproject, organism = srr_to_bioproject_and_organism_eutils( value.split(",")[0] )
        return groups, bioproject, organism
    else:
        return groups, None, None

def get_commit():
    if os.path.isfile(".commit") :
        with open(".commit", "r") as f:
            ver=f.readlines()[0][:7]
    elif os.path.isfile(".git/refs/heads/main"):
        with open(".git/refs/heads/main", "r") as f:
            ver=f.readlines()[0][:7]
    else:
        ver="<version>"

    return ver

def make_padj_hash_file( deseq2_output ):
    results_files=os.listdir(deseq2_output)
    results_files=[ s for s in results_files if s.endswith(".results.tsv") ]
    padj=[]
    for f in results_files :
        tmp=pd.read_csv(os.path.join( deseq2_output, f ), sep="\t" )
        tmp=tmp["padj"].tolist()
        padj=padj+tmp
    padj.sort()
    padj=pd.DataFrame({"padj":padj})
    padj.to_csv( os.path.join( deseq2_output, "padj.hash.tsv" ), sep="\t" ) 


def make_readme(var, md=None):
    ver=get_commit()
    f=os.path.basename(var["report_file"])
    org=var["organism"].replace("_", " ").capitalize()
    readme=f"""## Readme

Under {f} you will find a list of all relevant files from the analysis and to which section \
they belong to.

*fastqc*: fastqc results

*multiqc*: multiqc results

*deseq2*: differential gene expression results

*topgo*: topGO enrichment analysis of differentialy expressed genes (ie. padj < 0.05)

*qc_plots*: quality control of differential gene expression

*samples.mapping.xlsx* file to group and replicate mapping

## Material and methods

RNAseq analysis was performed using mpgagebioinformatics/jawm_rnaseq@{ver} and can be reproduced with: \
`jawm jawm_rnaseq@{ver} --global.var.geo <GEO ACCESSION>`. Shortly, \
reads were mapped to {org} reference genome (Ensembl {var["release"]}) using Kallisto (v0.48.0). \
After normalization of read counts by making use of the standard median-ratio for estimation \
of size factors, pair-wise differential gene expression was performed using DESeq2 (v1.38.0). \
After removal of genes with less then 10 overall reads log2 fold changes were shrank using \
approximate posterior estimation for GLM coefficients."""

    if md:
        with open(md, "w") as md_out:
            md_out.write(readme)

    return readme


download_link=jawm.Process( 
    name="download_link",
    when=lambda p: not os.path.isfile( os.path.join( p.var["source_folder"],"fastq.md5sums.check.log" ) ) ,
    script="""#!/bin/bash
mkdir -p {{source_folder}}
cd {{source_folder}}
{{download_link}}
find . -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
    -not -path "./*" -exec mv {} . \;
if [ -f fastq.md5sums.txt ] ; then m=fastq.md5sums.1.txt ; else  m=fastq.md5sums.txt ; fi
grep ".fastq.gz" {{md5sums_file}} > ${m}
grep ".fq.gz" {{md5sums_file}} >> ${m}
md5sum -c ${m} > fastq.md5sums.check.log
""",
    desc={
        "download_link":"Full command for download eg.: wget -r --http-user=NGS_25R480_doej_DX535 --http-passwd=your_delivery_password http://delivery.ccg.uni-koeln.de/downloads/NGS_25R480_doej_DX535",
        "source_folder": "Folder where the data will be downloaded to.",
    },
    # container=AGEPY_IMAGE,
)

preprocess=jawm.Process( 
    name="preprocess",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"],"preprocess.touch" ) ) ,
    script="""#!/usr/bin/env python3
import os
import gzip
import AGEpy as age
import re
from pathlib import Path
import pandas as pd
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

def concat_gz_files( files, output, concat="{{concat}}" ):

    if concat == "text": 
        with gzip.open(output, "wt") as outfile :
            for fname in files:
                with gzip.open(fname, "rt") as infile :
                    for line in infile:
                        outfile.write(line)

    elif concat == "stream" :
        with open(output, "wb") as outfile :
            for fname in files:
                with open(fname, "rb") as infile :
                    shutil.copyfileobj(infile, outfile, length=1024 * 1024)

def add_rep_suffix(groups):
    counts = {}
    output = []
    for g in groups:
        counts[g] = counts.get(g, 0) + 1
        output.append(f"rep_{counts[g]}")
    return output

groups = "{{groups}}".replace("\\x01", "" )
groups = groups.strip().strip(";")
groups = re.sub(r'\s*([;])\s*', r'\\1', groups)
groups = [s.split(";") for s in groups.split("\\n")]
groups = pd.DataFrame(groups, columns=["file_1", "file_2", "group"])
groups["group"] = groups["group"].apply(lambda x: age.safe_filename(x))
groups["rep"] = add_rep_suffix(groups["group"])


source_folder = "{{source_folder}}"
dest_folder = "{{raw_data}}"

if not os.path.isdir(dest_folder):
    os.makedirs(dest_folder)

groups.to_excel(os.path.join(dest_folder, "samples.mapping.xlsx"), index=False )


# We'll collect all concat jobs and run them with max 4 workers
futures = []

if int({{parallel}}) == 0 :
    w=1
else:
    w=int({{parallel}})

with ProcessPoolExecutor(max_workers=w  ) as executor:
    for i in groups.index.tolist():
        source_files_1 = groups.loc[i, "file_1"]
        source_files_2 = groups.loc[i, "file_2"]

        group = groups.loc[i, "group"]
        rep = groups.loc[i, "rep"]

        dest_path_1 = os.path.join(dest_folder, f"{group}.{rep}.read_1.fastq.gz")
        dest_path_2 = os.path.join(dest_folder, f"{group}.{rep}.read_2.fastq.gz")

        if os.path.exists(dest_path_1):
            print(f"{dest_path_1} already exists. File will not be overwritten. Using existing files.")
            continue

        # ---- Single-run case: just symlink/copy, no need for pool ----
        if "," not in source_files_1:
            # deal with one run per sample first
            source_name_1 = source_files_1
            source_path_1 = os.path.join(source_folder, source_files_1)

            if not os.path.isfile(source_path_1):
                raise FileNotFoundError(f"{source_path_1} not found")

            if str("{{copy}}").lower() != "true":
                os.symlink(source_path_1, dest_path_1)
            else:
                shutil.copy2(source_path_1, dest_path_1)

            if source_files_2:
                source_name_2 = source_files_2
                source_path_2 = os.path.join(source_folder, source_files_2)

                if not os.path.isfile(source_path_2):
                    raise FileNotFoundError(f"{source_path_2} not found")

                if str("{{copy}}").lower() != "true":
                    os.symlink(source_path_2, dest_path_2)
                else:
                    # NOTE: this used to copy source_path_1 by mistake
                    shutil.copy2(source_path_2, dest_path_2)

        # ---- Multi-run case: schedule concatenations in the pool ----
        else:
            # now deal with concatenations
            source_name_1 = [os.path.join(source_folder, s) for s in source_files_1.split(",")]
            futures.append(executor.submit(concat_gz_files, source_name_1, dest_path_1))

            if source_files_2:
                source_name_2 = [os.path.join(source_folder, s) for s in source_files_2.split(",")]
                futures.append(executor.submit(concat_gz_files, source_name_2, dest_path_2))

    # Wait for all concatenations to finish and raise exceptions if any
    for fut in as_completed(futures):
        fut.result()

# Only touch the file once everything is done
Path(os.path.join(dest_folder, "preprocess.touch")).touch()
""",
    var={
        "parallel": 4,
        "concat":"stream"
    },
    desc={
        "concat": "'stream or text for how you want to concat fastq.gz files when multiple files per sample exist. Default='stream' ",
        "groups": "A dictionary like structure '<file_read_1>:<file_read_2>:<group_name>;<file_read_1>:<file_read_2>:<group_name>;..' ",
        "source_folder": "Folder where the data is currently located.",
        "raw_data":"Samples with multiple runs will be concatenated into one file an located to this folder. Sample with only one file will be symlinked to this folder.", 
        "copy": "If True it will it will copy the original raw data file into its new location instead of creating symbolic links."
    },
    container=AGEPY_IMAGE,
)