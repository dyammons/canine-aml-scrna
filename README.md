# Canine AML analysis
Repository containing code assocaited with the analysis of canine AML dataset.

<br>

## Guidelines
- Create an aptly named branch
- Branches will be named scripts-initals (i.e. scripts-da | scripts-ah) unless this does not work well
- Work on the analysis then submit a pull request when you are ready to merge with main
- Do not commit directly to main (exception would be correcting typos etc)


## File structure:
- [:file\_folder: metaData](/metaData) contains relevant metadata files including cell type annotations and alternate sample names for analysis
- [:file\_folder: main](/) contains the scripts used in the analysis. One file for each major analysis step

<br>

## Setting up Git

```sh
#navigate to working directory
cd /pl/active/CSUClinHeme/users/dylan/proj02_k9_aml_scrna

#clone the repo -- only need to do this once (can use git fetch/pull once cloned)
git clone https://github.com/dyammons/canine-aml-scrna.git

#rename the repository
mv canine-aml-scrna scripts-ah

#enter repo and run/edit scripts
cd scripts-ah/

#create your own branch to work on!
git checkout -b ah
```

## Using the singularity contianer
To run an interactive R session:
```sh
singularity run -B $PWD/../../../ --fakeroot ../software/r4.3.2-seuratv5_v2
```

To submit a job that runs a R script:
```sh
singularity exec -B $PWD/../../../ --fakeroot ../software/r4.3.2-seuratv5_v2 Rscript script.R #change script.R to name of script to run
```


Idea is to clone this repo, pull down singularity container, then run the code!

<br>

## Directory structure
The project directory is located on the `pl` at `/pl/active/CSUClinHeme/users/dylan/proj02_k9_aml_scrna`.
The file structure is as follows:

```sh
proj02_k9_aml_scrna
├── external_data
│   └── any data from published (i.e. Austin’s data)
├── input
│   └── cellranger output count matrices
├── output
│   └── all output files from analysis
├── scripts-initials
│   ├── .git files
│   ├── logFiles
│   ├── metaData
│   └── where is repository should be cloned!
└── software
    └── singularity sandbox (and other containers)
```

