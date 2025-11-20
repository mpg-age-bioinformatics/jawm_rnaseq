# jawm_rnaseq

This is a jawm rnaseq module.

Installing jawm:
```
pip install git+ssh://git@github.com/mpg-age-bioinformatics/jawm.git
```
For more information on jawm please visit jawm's repo on [GitHub.com](https://github.com/mpg-age-bioinformatics/jawm/tree/main).

Example usage:
```
# clone this module
git clone git@github.com:mpg-age-bioinformatics/jawm_rnaseq.git

cd jawm_rnaseq

# download test data
jawm-test -r download

# docker
jawm rnaseq.py rnaseq -p ./yaml/docker.yaml

# slurm & apptainer with multiple yaml files
jawm rnaseq.py rnaseq -p ./yaml/vars.yaml ./yaml/hpc.yaml
```

Additional jawm workflows are available [here (GitHub.com)](https://github.com/mpg-age-bioinformatics?q=jawm_&type=all&language=&sort=).
