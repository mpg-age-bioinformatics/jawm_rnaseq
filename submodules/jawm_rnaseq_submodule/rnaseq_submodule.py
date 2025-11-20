import jawm

rnaseq_submodule_p1=jawm.Process( 
    name="rnaseq_submodule_p1",
    script="""#!/bin/bash
echo "Demo module echo process"
"""  
)