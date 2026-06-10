import jawm
import logging

logger = logging.getLogger("jawm_rnaseq")

if __name__ == "__main__":
    import sys
    import os
    from jawm.utils import workflow, load_modules, get_image, id_files
    from pathlib import Path
    import glob
    import string
    import random
    from datetime import datetime
    import csv
    import shutil


    # Start timestamp
    start = datetime.now()

    #####################
    # load dependencies #
    #####################

    # development
    load_modules([
        "submodules",
        "jawm_sra",
        "jawm_fastqc",
        "jawm_genomes",
        "jawm_kallisto",
        "jawm_featurecounts",
        "jawm_multiqc",
        "jawm_deseq2",
        ])

    # production
    # load_modules([
    #     "submodules",
    #     "jawm_kallisto@ceaba6c",
    #     "jawm_deseq2@e2827df",
    #     "jawm_genomes@697ce30",
    #     "jawm_featurecounts@61779d1",
    #     "jawm_fastqc@4b2c63d",
    #     "jawm_multiqc@5a14f42",
    #     "jawm_sra@46317a2",
    # ] )


    ##################
    # read arguments #
    ##################

    workflows, var, args, unknown_args= jawm.utils.parse_arguments( ["main","rnaseq","test", "mastertable"] )
    
    ############
    # set evns #
    ############

    if jawm.utils.docker_available():
        jawm.Process.update(override=False, environment="docker")
    elif jawm.utils.apptainer_available():
        jawm.Process.update(override=False, environment="apptainer")


    ############################
    # download required images #
    ############################

    if str(var.get("skip_prepull", False)).lower() == "false" :
    
        images=get_image( ) 

    ##################
    # process groups #
    ##################

    if var.get("geo", False) :
        logger.warning("""\
Parsing groups out of a GEO accession only works if the sample names finish with replicate infornation. Eg. wt 1, wt rep_1, wt sample 1, wt A, ..\
""")
        groups=_rnaseq.get_gse_groups( var["geo"])
        var["groups"]=groups

    if var.get("groups", False):
        groups, bioproject, organism = _rnaseq.process_groups( var["groups"] )
        jawm.Process.update( var={ "groups" : groups } )
        var["groups"]=groups

    else:
        bioproject=None
        organism=None

    # usage: 
    if workflow( ["main","rnaseq","test", "mastertable" ], workflows ) :

        #########################
        # check & set variables #
        #########################

        # organism
        if organism :
            jawm.Process.update( var={ "organism" : organism } )
            var["organism"]=organism

        # genome release
        if var.get("release") in (None, "last"):
            release, available=genomes.get_latest_ensembl_release(var['organism'])
            jawm.Process.update( var={ "release" : release } )
            var["release"]=release
        else:
            var["release"]=str(var["release"])

        # seq type: primary_assembly vs. toplevel (preference for primary_assembly )
        if "seq_type" not in var :
            seq_type=genomes.check_primary_assembly_exists( var['organism'], var['release'] )
            jawm.Process.update( var={ "seq_type" : seq_type } )
            var["seq_type"]=seq_type

        # output folders ie. genomes and rnaseq
        if not bioproject :
            chars = string.ascii_letters + string.digits
            chars = ''.join(random.choices(chars, k=7))
            default_location=os.path.join( os.getcwd(), "jawm_rnaseq_"+chars )
        else:
            default_location=os.path.join( os.getcwd(), "jawm_rnaseq_"+bioproject )
        
        # output folder
        if "project_folder" not in var :
            logger.info(f"Output folder: {default_location}")
            if os.path.isdir(default_location) :
                logger.warning(f"{default_location} already exists.")
            else:
                os.makedirs(default_location)
            jawm.Process.update( var={ "map.project_folder" : default_location } )
            var["project_folder"]=default_location

        # genomes folder
        if "genomes_folder" not in var:
            jawm.Process.update( var={ "map.genomes_folder" : var["project_folder"] } )
            var["genomes_folder"]= var["project_folder"]

        # BioMart biotypes and go annotation file
        if "biotypes_go" not in var :
            biotypes_go = os.path.join( var["genomes_folder"], var["organism"], var["release"], "biotypes_go.txt" )
            jawm.Process.update( var={ "mk.biotypes_go_dir" : os.path.dirname( biotypes_go ) } )
            jawm.Process.update( var={ "biotypes_go" : biotypes_go } )

        # user's raw data folder
        if ( not bioproject ) and ( "source_folder" not in var ) :
            if "download_link" not in var :
                jawm.Process.update( var={ "mk.source_folder" : "./" } )
            else:
                jawm.Process.update( var={ "mk.source_folder" :  os.path.join( var["project_folder"], "downloads"  ) }) 

        jawm.Process.update( var={
            "read1_suffix" : ".read_1.fastq.gz",
            "read2_suffix": ".read_2.fastq.gz",
            "mk.raw_data": os.path.join( var["project_folder"], "raw_data" )
            } )

        var.update({
            "read1_suffix":".read_1.fastq.gz",
            "read2_suffix":".read_2.fastq.gz"
        })

        # readme file path
        if not var.get("readme_file", False) :
            var["readme_file"]=os.path.join(var["project_folder"], "README.md")

        ##################
        # kallisto index #
        ##################


        # download references
        genomes.download_gtf.var["filename"] = genomes.download_gtf.var["organism"] + "." + str( genomes.download_gtf.var["release"] ) + ".gtf.gz"
        genomes.download_gtf.var["mk.output_folder"] = os.path.join( var["genomes_folder"], var["organism"], var["release"] )
        genomes.download_gtf.execute()

        genomes.download_dna.var["filename"] = genomes.download_dna.var["organism"] + "." + str(genomes.download_dna.var["release"]) + "." + genomes.download_dna.var["seq_type"] + ".dna.fa.gz"
        genomes.download_dna.var["mk.output_folder"] = os.path.join( var["genomes_folder"], var["organism"], var["release"] )
        genomes.download_dna.execute()

        gtf_file=os.path.join( 
            genomes.download_gtf.var["mk.output_folder"],
            f'{var["organism"]}.{var["release"]}.no.rRNA.gtf' 
            )

        kallisto_dic={ 
            "mk.kallisto_index":os.path.join( genomes.download_gtf.var["mk.output_folder"], "kallisto_index" ),
            "mk.kallisto_output":os.path.join( var["project_folder"], "kallisto_output" ),
            "gtf_file": gtf_file,
            "fasta_file":os.path.join( 
                genomes.download_gtf.var["mk.output_folder"], 
                f'{var["organism"]}.{var["release"]}.{var["seq_type"]}.dna.fa' 
                )
        }

        jawm.Process.update( var=kallisto_dic )
        
        # generate indexes
        kallisto.get_genome.execute( [ genomes.download_gtf.hash, genomes.download_dna.hash ] )
        kallisto.writecdna.execute( kallisto.get_genome.hash )
        kallisto.indexer.execute( kallisto.writecdna.hash )

        ####################
        # collect raw data #
        ####################

        if not bioproject :

            ##############
            # local data #
            ##############

            download_jobs=[]
            if "download_link" in var :
                _rnaseq.download_link.execute()
                download_jobs.append( _rnaseq.download_link.hash ) 

                # special handling of samples from the Cologne Center for Genomics
                if ( "fastq.gz" not in var["groups"] ) and ("fq.gz" not in var["groups"]) and ( "delivery.ccg.uni-koeln.de" in var["download_link"] ) :
                    jawm.Process.wait( _rnaseq.download_link.hash )
                    groups=_rnaseq.ccg_samples( groups, source_folder )
                    jawm.Process.update( var={ "groups" : groups } )
                    var["groups"]=groups

            # for some samples there might be more than one file per read
            # which we will concatenate
            n_concats=var["groups"].split('\n')
            n_concats=[ s.split(";")[0] for s in n_concats ]+ [ s.split(";")[1] for s in n_concats ]
            n_concats=[ s for s in n_concats if "," in s ]
            n_concats=len(n_concats)

            if n_concats < 5 :
                _rnaseq.preprocess.var["parallel"]=n_concats
                _rnaseq.preprocess.manager="local"
            elif _rnaseq.preprocess.manager=="slurm" :
                _rnaseq.preprocess.var["parallel"]=n_concats
                _rnaseq.preprocess.manager_slurm={"--mem":"20GB", "-t":"1:00:00", "-c":str(n_concats) }

            _rnaseq.preprocess.execute( download_jobs )

            preprocess_jobs=[ _rnaseq.preprocess.hash ]
        
        else:

            ##########################
            # download data from sra #
            ##########################
                    
            sra.read_acc.execute()
            jawm.Process.wait( sra.read_acc.hash )

            # input_xlsx=os.path.join( sra.read_acc.var["raw_data"], f"{sra.read_acc.var['acc']}.samples.xlsx"  )
            sra_ids = []
            input_tsv=os.path.join( sra.read_acc.var["raw_data"], f"sra.samples.tsv"  )
            sra_ids = sra.read_samples_file( os.path.join( sra.read_acc.var["raw_data"], "sra.samples.tsv"  ) )

            fastq_dump_jobs=[]
            for sraid in sra_ids : 
                
                # download serially to prevent too much io and server query
                prefetch_=sra.prefetch.clone()
                prefetch_.var["sraid"]=sraid
                prefetch_.execute()

                if prefetch_.manager == "local" :

                    jawm.Process.wait( prefetch_.hash )

                # send extraction jobs to background
                fastq_dump_=sra.fastq_dump.clone()
                fastq_dump_.var["sraid"]=sraid
                fastq_dump_.execute( prefetch_.hash )

                fastq_dump_jobs.append( fastq_dump_.hash )
            
            jawm.Process.wait( fastq_dump_jobs )

            # for some samples there might be more than one SRA per sample
            # which we will concatenate
            n_concats=sra.get_nconcatenations( input_tsv, sra.read_acc.var["raw_data"] )

            # define number of parallel concatenations
            if n_concats < 5 :
                sra.relabel_sra.var["parallel"]=n_concats
                sra.relabel_sra.manager="local"
            elif sra.relabel_sra.manager=="slurm" :
                if "max_concats" in var :
                    if n_concats > int(var["max_concats"] ) :
                        n_concats=int(var["max_concats"] )
                else:
                     n_concats=10

                sra.relabel_sra.var["parallel"]=n_concats
                sra.relabel_sra.manager_slurm.update(  {"--mem":"20GB", "-t":"1:00:00", "-c":str(n_concats) } )

            sra.relabel_sra.execute( fastq_dump_jobs )

            preprocess_jobs= [ sra.relabel_sra.hash ] 


        #################
        # list raw data #
        #################

        jawm.Process.wait( preprocess_jobs )

        # list read1 and read2 files
        read1_files = list(Path( fastqc.fastqc.var["raw_data"] ).glob(f'*{fastqc.fastqc.var["read1_suffix"]}'))
        read1_files=[ str(s) for s in read1_files if not str(s).startswith("tmp.") ]
        read2_files = [ str(s).split( var["read1_suffix"] )[0] +var["read2_suffix"] for s in read1_files  ]

        ##########
        # fastqc #
        ##########

        jawm.Process.update( var={ "mk.fastqc_output":os.path.join( var["project_folder"], "fastqc_output" ) } )

        fastqc_jobs=[]
        for f in read1_files + read2_files :
            fastqc_=fastqc.fastqc.clone()
            fastqc_.var["map.f"]=f

            ## debug online failure
            p = Path(f)
            if not p.is_file():
                print(f"ERROR: FastQC input is not a file: {p}", file=sys.stderr)
                print(f"exists={p.exists()} is_dir={p.is_dir()}", file=sys.stderr)
                sys.exit(1)

            fastqc_.execute( preprocess_jobs )
            fastqc_jobs.append(fastqc_.hash)

        ####################
        # kallisto mapping #
        ####################

        # check if it is unstranded, 1st strand, or 2nd strand mapping
        kallisto.unstranded_mapping.var["read_1"]=os.path.basename(read1_files[0])
        kallisto.unstranded_mapping.execute( kallisto.indexer.hash )
        
        kallisto.gene_model.execute( kallisto.get_genome.hash )
        
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
            "map.strand":os.path.join( var["project_folder"], "kallisto_output", ".strandness.featurecounts.txt"),
            }

        jawm.Process.update( var=featurecounts_dic )

        # run jobs

        id_bams=id_files( os.path.join( var["project_folder"], "kallisto_output") )
        id_bams={k: v for k, v in id_bams.items() if not k.startswith("tmp.")}

        feature_counts_jobs=[]

        for pair_id in id_bams.keys() :

            geneid_=featurecounts.geneid.clone()
            geneid_.var["pair_id"]=pair_id
            geneid_.var["map.bam"]=id_bams[pair_id][0]

            biotype_=featurecounts.biotype.clone()
            biotype_.var["pair_id"]=pair_id
            biotype_.var["map.bam"]=id_bams[pair_id][0]

            if os.path.isfile( os.path.join( fastqc.fastqc.var["raw_data"], pair_id+var["read2_suffix"] ) ) :
                geneid_.var["paired"]="-p"
                biotype_.var["paired"]="-p"

            geneid_.execute( flagstat_jobs )
            biotype_.execute( flagstat_jobs )
            
            feature_counts_jobs = feature_counts_jobs + [ geneid_.hash, biotype_.hash ]

        if os.path.basename( list(id_bams.values())[0][0] ) == "pseudoalignments.bam" :

            featurecounts.headers.execute( feature_counts_jobs )

            feature_counts_jobs.append( featurecounts.headers.hash )
        
        ###########
        # multiqc #
        ###########

        multiqc_dic={
            "mk.fastqc": os.path.join( var["project_folder"], "fastqc_output" ),
            "mk.featurecounts": os.path.join( var["project_folder"], "featurecounts_output" ),
            "map.mapping":  os.path.join( var["project_folder"], "kallisto_output" ),
            "mk.multiqc_output": os.path.join( var["project_folder"], "multiqc_output" )
        }

        jawm.Process.update( var=multiqc_dic )

        multiqc.multiqc.execute( fastqc_jobs + feature_counts_jobs  )

        ##########
        # deseq2 #
        ##########

        deseq2_dic={
            "mk.deseq2_output": os.path.join( var["project_folder"], "deseq2_output" ),
            "map.kallisto_output": os.path.join( var["project_folder"], "kallisto_output" ),
            "map.gtf": gtf_file
        }

        jawm.Process.update( var=deseq2_dic )

        deseq2.tx2gene.execute( )
        deseq2.annotations.execute()
        deseq2.parse_submission.var["rep_prefix"]=".rep_"
        deseq2.parse_submission.execute()
        deseq2.tx2gene_proc.execute( deseq2.tx2gene.hash )

        jawm.Process.wait( [ deseq2.tx2gene_proc.hash, deseq2.parse_submission.hash ] )

        deseq2_output = deseq2.deseq2.var["deseq2_output"]
        tests = glob.glob( os.path.join( deseq2_output, "*input.tsv" ) )

        deseq2_jobs=[]

        if not workflow( ["mastertable"], workflows ) :
            # when we do mastertable we do not want to pair-wise analysis
            for file in tests:

                # clone the processes required for each file
                deseq2_=deseq2.deseq2.clone()

                deseq2_.var["input_file"]=os.path.basename(file)
                deseq2_.execute()

                deseq2_jobs.append( deseq2_.hash )

        deseq2.mastertable.execute( deseq2_jobs )

        deseq2.annotator.execute( deseq2.mastertable.hash )

        if workflow( ["mastertable"], workflows ) :
            jawm.Process.wait()
            sys.exit()

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

        jawm.Process.wait( topgo_jobs + [ deseq2.david.hash ] )

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

    # create a reference file for hashing
    _rnaseq.make_padj_hash_file( deseq2_output )

    # write readme file
    readme=_rnaseq.make_readme(var, var["readme_file"] )

    # wait for all jobs to finish
    jawm.Process.wait( )

    # list most relevant files
    if not var.get("report_file", False) :
        var["report_file"]=os.path.join(var["project_folder"], "report_file.txt")
    var["report_file_"]=os.path.join(var["project_folder"], "report_file.txt")


    def add_files(report={}, add={} ):
        for k in add :
            report[k]=report.get(k, []) + add[k]
        return report
    
    report={
        "main" : [ os.path.join( kallisto.unstranded_mapping.var["raw_data"], "samples.mapping.xlsx" ) ,
            var["readme_file"],
            var["report_file_"]
            ]
    }

    report=add_files(report=report, add=fastqc.report_files( fastqc.fastqc.var["fastqc_output"] ) )
    report=add_files(report=report, add=multiqc.report_files( multiqc.multiqc.var["multiqc_output"] ) )
    report=add_files(report=report, add=deseq2.report_files( deseq2.deseq2.var["deseq2_output"] ) )


    with open( f'{var["report_file"]}_', "w") as f :
        for folder, file_list in report.items() :
            for file in file_list :
                f.write(f"{folder} {file}\n")

    if var["report_file_"] != var["report_file"] :
        shutil.copyfile( f'{var["report_file"]}_', var["report_file_"] )
    
    os.rename(f'{var["report_file"]}_', f'{var["report_file"]}')

    if workflow( "test", workflows ) :
        # for the test workflow we also do something more (just for rnaseq)
        var["deseq2_output"]=deseq2.deseq2.var["deseq2_output"]
        deseq2.round_master_table(var)
        logger.info("Test completed.")

    ###################
    # report run time #
    ###################

    # End timestamp
    end = datetime.now()

    # Compute elapsed time
    elapsed = end - start

    # Convert to days, hours, minutes
    days = elapsed.days
    hours, remainder = divmod(elapsed.seconds, 3600)
    minutes, _ = divmod(remainder, 60)

    logger.info(f"Elapsed time: {days}d {hours}h {minutes}m")

    sys.exit(0)
