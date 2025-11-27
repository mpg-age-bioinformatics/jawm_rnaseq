import jawm
import os
import re
import pandas as pd

preprocess=jawm.Process( 
    name="preprocess",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"],"preprocess.touch" ) ) ,
    script="""#!/usr/bin/env python3
import os
import AGEpy as age
import re
from pathlib import Path
import pandas as pd
import shutil

def concat_gz_files(files, output):
    with gzip.open(output, "wt") as outfile:
        for fname in files:
            with gzip.open(fname, "rt") as infile:
                for line in infile:
                    outfile.write(line)

def add_rep_suffix(groups):
    counts = {}
    output = []
    for g in groups:
        counts[g] = counts.get(g, 0) + 1
        output.append(f"rep_{counts[g]}")
    return output

groups="{{groups}}".replace("\\x01", "" )
groups=groups.strip().strip(";")
groups = re.sub(r'\s*([;])\s*', r'\\1', groups )
groups=[ s.split(";") for s in groups.split("\\n") ]
groups=pd.DataFrame(groups, columns=["file_1", "file_2","group"] )
groups["group"] = groups["group"].apply(lambda x: age.safe_filename(x) )
groups["rep"]=add_rep_suffix( groups["group"] )

source_folder="{{source_folder}}"
dest_folder="{{raw_data}}"

if not os.path.isdir(dest_folder):
    os.makedirs(dest_folder)

for i in groups.index.tolist():
    source_files_1=groups.loc[i,"file_1"]
    source_files_2=groups.loc[i,"file_2"]

    group=groups.loc[i,"group"]
    rep=groups.loc[i,"rep"]

    dest_path_1=os.path.join(dest_folder,f"{group}.{rep}.read_1.fastq.gz")
    dest_path_2=os.path.join(dest_folder,f"{group}.{rep}.read_2.fastq.gz")

    if os.path.exists(dest_path_1) : 
        print(f"{dest_path_1} already exists. File will not be overwritten. Using existing files.")
        continue

    if "," not in source_files_1:

        # deal with one run per sample first
        source_name_1=source_files_1
        source_path_1=os.path.join(source_folder,source_files_1)

        if str("{{copy}}").lower() != "true" :
            os.symlink(source_path_1, dest_path_1)
        else:
            shutil.copy2(source_path_1,dest_path_1)

    
        if source_files_2 :

            source_name_2=source_files_2
            source_path_2=os.path.join(source_folder,source_files_2)

            if str("{{copy}}").lower() != "true" :  
                os.symlink(source_path_2, dest_path_2)
            else:
                shutil.copy2(source_path_1,dest_path_1)
        
    else:

        # now deal with concatenations
        source_name_1=[ os.path.join(source_folder,s) for s in source_files_1.split(",") ]
        concat_gz_files(source_name_1, dest_path_1)
        if source_files_2 :
            source_name_2=[ os.path.join(source_folder,s) for s in source_files_2.split(",") ]
            concat_gz_files(source_name_2, dest_path_2)

Path( os.path.join( dest_folder,"preprocess.touch" ) ).touch()
""",
    desc={
        "groups": "A dictionary like structure '<file_read_1>:<file_read_2>:<group_name>;<file_read_1>:<file_read_2>:<group_name>;..' ",
        "source_folder": "Folder where the data is currently located.",
        "raw_data":"Samples with multiple runs will be concatenated into one file an located to this folder. Sample with only one file will be symlinked to this folder.", 
        "copy": "If True it will it will copy the original raw data file into its new location instead of creating symbolic links."
    },
    container="mpgagebioinformatics/agepy:441ee35",
)

if __name__ == "__main__":
    import sys
    from jawm.utils import workflow, load_modules, get_image, id_files
    from pathlib import Path

    # load modules from online git repos
    load_modules([
        "jawm_sra",
        "jawm_fastqc",
        "jawm_genomes",
        "jawm_kallisto",
        "jawm_featurecounts",
        "jawm_multiqc",
        "jawm_deseq2",
        ])

    workflows, var, args, unknown_args= jawm.utils.parse_arguments(["main","rnaseq","test"])

    # usage: 

    if workflow( ["main","rnaseq","test"], workflows ) :

        ############################
        # download required images #
        ############################

        images=get_image( [ 
            preprocess.container,
            sra.read_geo.container, 
            sra.prefetch.container,
            kallisto.get_erccs.container,
            kallisto.writecdna.container,
            kallisto.indexer.container,
            kallisto.gene_model.container,
            multiqc.multiqc.container
            ] )

        ####################
        # collect raw data #
        ####################

        # check if sra data must be download or if the user is analysing local data
        # all geo analysis must have sra.read_geo.var["geoacc"] defined

        sra_download=False
        if "geoacc" in sra.read_geo.var :
            if sra.read_geo.var["geoacc"] :
                sra_download=True

        jawm.Process.update( var={
            "read1_suffix" : ".read_1.fastq.gz",
            "read2_suffix": ".read_2.fastq.gz",
            "mk.raw_data": os.path.join( var["project_folder"], "raw_data" )
            } )

        if not sra_download :

            ##############
            # local data #
            ##############

            preprocess.execute()
            jawm.Process.wait( preprocess.hash )
        
        else:

            ##########################
            # download data from sra #
            ##########################    
                    
            sra.read_geo.execute()
            jawm.Process.wait( sra.read_geo.hash )

            input_xlsx=os.path.join( sra.read_geo.var["raw_data"], f"{sra.read_geo.var['geoacc']}.samples.xlsx"  )
            df=pd.read_excel( input_xlsx )
            sra_ids=df["Run"].tolist()

            fastq_dump_jobs=[]
            for sraid in sra_ids : 
                
                # download serially to prevent too much io and server query
                prefetch_=sra.prefetch.clone()
                prefetch_.var["sraid"]=sraid
                prefetch_.execute()
                jawm.Process.wait( prefetch_.hash )

                # send extraction jobs to background
                fastq_dump_=sra.fastq_dump.clone()
                fastq_dump_.var["sraid"]=sraid
                fastq_dump_.execute(  )

                fastq_dump_jobs.append( fastq_dump_.hash )
            
            jawm.Process.wait( fastq_dump_jobs )

            sra.relabel_geo.execute( )

            jawm.Process.wait( sra.relabel_geo.hash )

        #################
        # list raw data #
        #################

        # list read1 and read2 files
        read1_files = list(Path( fastqc.fastqc.var["raw_data"] ).glob(f'*{fastqc.fastqc.var["read1_suffix"]}'))
        read1_files=[ str(s) for s in read1_files if not str(s).startswith("tmp.") ]
        read2_files = [ str(s).split( var["read1_suffix"] )[0] +var["read2_suffix"] for s in read1_files  ]

        if not all( os.path.isfile(s) for s in read2_files  ) :
            read2_files = [ ]
            print("Not all Read 1 files have a matching Read 2 file. Treating experiment as single end.")
        else:
            print("Treating experiment as paired end.")

        ##########
        # fastqc #
        ##########

        jawm.Process.update( var={ "mk.fastqc_output":os.path.join( var["project_folder"], "fastqc_output" ) } )

        fastqc_jobs=[]
        for f in read1_files + read2_files :
            fastqc_=fastqc.fastqc.clone()
            fastqc_.var["map.f"]=f
            fastqc_.execute()
            fastqc_jobs.append(fastqc_.hash)

        ############
        # kallisto #
        ############

        gtf_file=os.path.join( 
            var["genomes_folder"], 
            var["organism"], 
            var["release"], 
            f'{var["organism"]}.{var["release"]}.no.rRNA.gtf' 
            )

        kallisto_dic={ 
            "mk.kallisto_index":os.path.join( var["project_folder"], "kallisto_index" ),
            "mk.kallisto_output":os.path.join( var["project_folder"], "kallisto_output" ),
            "gtf_file": gtf_file,
            "fasta_file":os.path.join( 
                var["genomes_folder"], 
                var["organism"], 
                var["release"], 
                f'{var["organism"]}.{var["release"]}.{var["seq_type"]}.dna.fa' 
                )
        }

        jawm.Process.update( var=kallisto_dic )

        # download references
        genomes.download_gtf.var["filename"] = genomes.download_gtf.var["organism"] + "." + str(genomes.download_gtf.var["release"]) + ".gtf.gz"
        genomes.download_gtf.var["mk.output_folder"] = os.path.join( var["project_folder"], var["organism"], var["release"] )
        genomes.download_gtf.execute()

        genomes.download_dna.var["filename"] = genomes.download_dna.var["organism"] + "." + str(genomes.download_dna.var["release"]) + "." + genomes.download_dna.var["seq_type"] + ".dna.fa.gz"
        genomes.download_dna.var["mk.output_folder"] = os.path.join( var["project_folder"], var["organism"], var["release"] )
        genomes.download_dna.execute()
        
        # generate indexes
        kallisto.get_erccs.execute( [ genomes.download_gtf.hash, genomes.download_dna.hash ] )
        kallisto.writecdna.execute( kallisto.get_erccs.hash )
        kallisto.indexer.execute( kallisto.writecdna.hash )

        # check if it is unstranded, 1st strand, or 2nd strand mapping
        kallisto.unstranded_mapping.var["read_1"]=os.path.basename(read1_files[0])
        kallisto.unstranded_mapping.execute( kallisto.indexer.hash )
        
        kallisto.gene_model.execute( kallisto.get_erccs.hash )
        
        kallisto.infer_experiment.var["read_1"]=os.path.basename( read1_files[0] )
        kallisto.infer_experiment.execute( [ kallisto.gene_model.hash, kallisto.unstranded_mapping.hash ] )

        kallisto.strand_checker.execute( kallisto.infer_experiment.hash )

        jawm.Process.wait( kallisto.strand_checker.hash )

        # running mapping jobs
        flagstat_jobs=[]
        mapping_jobs=[]
        for read1 in read1_files:

            # clone the processes required for each file
            mapping_=kallisto.mapping.clone()
            flagstat_=kallisto.flagstat.clone()

            mapping_.var["read_1"]=os.path.basename( read1 )
            mapping_.execute( kallisto.strand_checker.hash )
            mapping_jobs.append( mapping_.hash )

            flagstat_.var["read_1"]=os.path.basename( read1 )
            flagstat_.execute( mapping_.hash ) 

            flagstat_jobs.append( flagstat_.hash )

        jawm.Process.wait( mapping_jobs )

        #################
        # featurecounts #
        #################

        featurecounts_dic={
            "mk.featurecounts_output": os.path.join( var["project_folder"], "featurecounts_output"),
            "map.mapping_folder": os.path.join( var["project_folder"], "kallisto_output"),
            "map.gtf":gtf_file,
            "map.strand_file":os.path.join( var["project_folder"], "kallisto_output", ".strandness.txt")
        }

        jawm.Process.update( var=featurecounts_dic )

        if read2_files:
            featurecounts.geneid.var["paired"]="-p"
            featurecounts.biotype.var["paired"]="-p"

        # run jobs

        id_bams=id_files( os.path.join( var["project_folder"], "kallisto_output") )
        feature_counts_jobs=[]

        for pair_id in id_bams.keys() : 

            geneid_=featurecounts.geneid.clone()
            geneid_.var["pair_id"]=pair_id
            geneid_.var["map.bam"]=id_bams[pair_id][0]
            geneid_.execute( flagstat_jobs )

            biotype_=featurecounts.biotype.clone()
            biotype_.var["pair_id"]=pair_id
            biotype_.var["map.bam"]=id_bams[pair_id][0]
            biotype_.execute( flagstat_jobs )
            
            feature_counts_jobs = feature_counts_jobs + [ geneid_.hash, biotype_.hash ]

        if os.path.basename( list(id_bams.values())[0][0] ) == "pseudoalignments.bam" :

            featurecounts.headers.execute( feature_counts_jobs )

            feature_counts_jobs.append( featurecounts.headers.hash )
        
        ###########
        # multiqc #
        ###########

        multiqc_dic={
            "map.fastqc": os.path.join( var["project_folder"], "fastqc_output" ),
            "map.featurecounts": os.path.join( var["project_folder"], "featurecounts_output" ),
            "map.mapping":  os.path.join( var["project_folder"], "kallisto_output" ),
            "mk.multiqc_output": os.path.join( var["project_folder"], "multiqc_output" )
        }

        jawm.Process.update( var=multiqc_dic )

        multiqc.multiqc.execute( fastqc_jobs + feature_counts_jobs  )

        ##########
        # deseq2 #
        ##########

        deseq2_dic={
            "map.deseq2_output": os.path.join( var["project_folder"], "deseq2_output" ),
            "map.kallisto_output": os.path.join( var["project_folder"], "kallisto_output" ),
            "map.gtf": gtf_file
        }

        jawm.Process.update( var=deseq2_dic )

        deseq2.tx2gene.execute( )
        deseq2.annotations.execute()
        deseq2.parse_submission.execute()
        deseq2.tx2gene_proc.execute( tx2gene.hash )

        jawm.Process.wait([ deseq2.tx2gene_proc.hash, deseq2.parse_submission.hash ])

        deseq2_output = deseq2.deseq2.var["deseq2_output"]
        tests = glob.glob( os.path.join( deseq2_output, "*input.tsv" ) )

        deseq2_jobs=[]

        for file in tests:

            # clone the processes required for each file
            deseq2_=deseq2.deseq2.clone()

            deseq2_.var["input_file"]=os.path.basename(file)
            deseq2_.execute()

            deseq2_jobs.append( deseq2_.hash )
        
        deseq2.mastertable.execute( )
        deseq2_jobs.append( deseq2.mastertable.hash )

        deseq2.annotator.execute( deseq2_jobs )

        # wait for annotator to complete before listing input files and starting david
        # as they are generated by annotator
        jawm.Process.wait( deseq2.annotator.hash )

        deseq2.david.execute( )

        annotated_output = os.path.join(deseq2_output, "annotated")
        tests_results = glob.glob(os.path.join(annotated_output, "*results.tsv"))
        topgo_jobs=[]
        for file in tests_results:
            topgo_=deseq2.topgo.clone()
            topgo_.var["input_file"]=os.path.basename(file)
            topgo_.execute( )
            topgo_jobs.append( topgo_.hash )

        jawm.Process.wait( topgo_jobs + [ david.hash ] )

        david_files=glob.glob(os.path.join(annotated_output, "*.DAVID.tsv"))
        togo_files=glob.glob(os.path.join(annotated_output, "*.topGO.tsv"))

        deseq2.cellplot.var["filetype"]="tsv"
        deseq2.cellplot.var["nterms"]="15"

        cellplot_jobs=[]
        for file in david_files:
            cellplot_=deseq2.cellplot.clone()

            cellplot_.var["map.inFile"]=file
            cellplot_.var["category"]="GOTERM_BP_FAT"

            cellplot__=cellplot_.clone()
            cellplot__.var["category"]="KEGG_PATHWAY"

            cellplot_.execute()
            cellplot__.execute()

            cellplot_jobs.append( cellplot_.hash )
            cellplot_jobs.append( cellplot__.hash )

        for file in togo_files:

            cellplot_=deseq2.cellplot.clone()

            cellplot_.var["map.inFile"]=file
            cellplot_.var["category"]="GOTERM_BP"
            cellplot_.execute()

            cellplot_jobs.append( cellplot_.hash )

        deseq2.qc_plots.execute( )


    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for rnaseq)
        print("Test completed.")


    sys.exit(0)
