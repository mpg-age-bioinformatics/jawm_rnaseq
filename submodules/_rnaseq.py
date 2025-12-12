import jawm
import os
import re
import pandas as pd
import random
import string
from io import StringIO
from pathlib import Path

AGEPY_IMAGE="mpgagebioinformatics/agepy:9acb5de"

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
    
def load_table(file_path, n):
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

    if len( df.columns.tolist() ) < n :
        raise ValueError(f"Imported table had less than 2 columns.")

    if len( df.columns.tolist() ) > n :

        cols = list(df.columns)

        # special case sra groups
        if n == 2 :
            
            cols_=[ s for s in cols if s == "Experiment" ]+[ s for s in cols if s in ["group", "groups"] ]

            if len(cols_) == n :
                df=df[cols_]

            exp_idx = cols.index("Experiment")

            # Check that a column exists after Experiment
            if exp_idx + 1 >= len(cols):
                raise ValueError("Column 'Experiment' exists but there is no column after it.")

            # Select Experiment and the following column
            df=df[ cols[exp_idx : exp_idx + 2] ]
    
            df.columns=["sample","group"]

        else :
            df=df[ cols[:n] ]

    return df.to_csv(sep=";", index=False, header=False ).rstrip("\n").replace("\n", r"\n")


def process_groups( groups, sra_download ) :

    if sra_download :
        # sra groups are given in the form
        # <samples list><groups_list> 
        n=2
    else:
        # own data groups are given in the form
        # <list of read 1 files><list of read 2 files><groups_list> 
        n=3

    if os.path.isfile(groups) :
        groups=load_table(groups, n)
    
    else:

        groups=groups.split("\n")
        groups=[ s.replace(";","_") for s in groups if s ]
        
        if n != 2:
            # own data
            if ( len(groups) % 2 == 0) :
                # single reads
                n=2
            
            if ( len(groups) % 3 == 0) :
                # paired reads
                n=3

        size = len(groups) // n
        groups = [ groups[ i*size : (i+1)*size ] for i in range(n) ]

        if sra_download:
            samples=groups[0]
            groups=groups[-1]
            groups=[ f"{a};{b}" for a, b in zip(samples, groups) ]
        else:
            samples=groups[0]
            samples_=groups[0]
            groups=groups[-1]
            groups=[ f"{a};{b};{c}" for a, b, c in zip( samples, samples_, groups ) ]

        groups=r"\n".join(groups)

    return groups

download_link=jawm.Process( 
    name="download_link",
    when=lambda p: not os.path.isfile( os.path.join( p.var["source_folder"],"fastq.md5sums.check.log" ) ) ,
    script="""#!/bin/bash
mkdir -p {{source_folder}}
cd {{source_folder}}
{{download_link}}
find . -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
    -not -path "./*" -exec mv {} . \;
grep ".fastq.gz" {{md5sums_file}} > fastq.md5sums.txt 
grep ".fq.gz" {{md5sums_file}} >> fastq.md5sums.txt
md5sum -c fastq.md5sums.txt > fastq.md5sums.check.log
""",
    desc={
        "download_link":"Full command for download eg.: wget -r --http-user=NGS_25R480_doej_DX535 --http-passwd=your_delivery_password http://delivery.ccg.uni-koeln.de/downloads/NGS_25R480_doej_DX535",
        "source_folder": "Folder where the data will be downloaded to.",
    },
    container=AGEPY_IMAGE,
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

            if str("{{copy}}").lower() != "true":
                os.symlink(source_path_1, dest_path_1)
            else:
                shutil.copy2(source_path_1, dest_path_1)

            if source_files_2:
                source_name_2 = source_files_2
                source_path_2 = os.path.join(source_folder, source_files_2)

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