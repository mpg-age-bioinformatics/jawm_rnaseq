import jawm

rnaseq_p1=jawm.Process( 
    name="rnaseq_p1",
    script="""#!/bin/bash
echo "{{extra_args}} {{my_rnaseq_argument}}" 2>&1 | tee {{output}}/rnaseq.txt
""",

    desc={
        "extra_args": "use this if you want to add not preset arguments",
        "my_rnaseq_argument":"Some text that will be printed to the screen.", 
        "output":"Folder for output", # the prefix "mk." leads to the creation of this folder and volume mapping if you are using containers
    },
    
    container="mpgagebioinformatics/fastqc:0.11.9",
  
)


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow
    from jawm.utils import load_modules

    # load modules from local folders
    load_modules("submodules")

    # load modules from online git repos
    load_modules("jawm_template")

    # it can be used in the form
    # load_modules(["modules","jawm_template@<tag/commit_hash>"])
    
    # or for latest available tag
    # load_modules("jawm_template@latest")

    workflows, var, args, unknown_args= jawm.utils.parse_arguments(["main","rnaseq","test"],)

    # usage: 

    if workflow( ["main","rnaseq","test"], workflows ) :

        # execute process
        rnaseq_p1.execute()

        # execute a process with dependencies
        rnaseq_p2.depends_on=[rnaseq_p1.hash]
        rnaseq_p2.execute()

        # wait for all above processes to complete
        jawm.Process.wait(rnaseq_p2.hash)

        # print the output
        print(rnaseq_p1.get_output())
        print(rnaseq_p2.get_output())

        rnaseq_submodule.rnaseq_submodule_p1.execute()
        template._template_p2.execute()

        jawm.Process.wait([ rnaseq_submodule.rnaseq_submodule_p1.hash, template._template_p2.hash  ])
        print(rnaseq_submodule.rnaseq_submodule_p1.get_output())
        print(template._template_p2.get_output())

    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for rnaseq)
        rnaseq_p3.execute()
        jawm.Process.wait( rnaseq_p3.hash)
        print("Test completed.")


    sys.exit(0)
