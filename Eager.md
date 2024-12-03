# Day 1
## Eager
Eager is an ancient DNA analysis pipeline. Please find the link to the website: [Eager 2.5.0](https://nf-co.re/eager/2.5.0). It can be used for analysing Human, Animals, Plants, Microbiomes. This pipeline uses a workflow tool called Nextflow (https://www.nextflow.io/).

The aim of this pipeline is to process the data from FASTQ format inputs or preprocessed Bam files inputs to final Bam files outputs. This pipeline also effectates different tasks such as: sequencing quality control (FastQC), removing adapters, alignment, sex determination, damage profil, contamination, ...
## Set-up step
### Working environment in Mjolnir
#### Step 1: Key elements for Mjolnir

For the first connection, you will have to log on with a VPN, then log into Mjolnir with the following command:

```
ssh <KU_ID>@mjolnirhead01fl.unicph.domain
```
Your password is your KU password.

Please find below important paths in Mjolnir to know before starting:

```
/projects/course_sgbb20001/data ==> data directory including shared softwares, databases, etc
/projects/course_sgbb20001/people/USERID ==> your personal folder

```
For the first time Mjolnir user, please follow the following steps before starting any analysis.

1. Copy the standard .bashrc and .bash_profile files to your home directory:

```
nano ~/.bashrc
alias ll="ls -lhtra"
alias xpinfo="squeue -S 'P,t,-p,i' --format '%.20i %.10Q %.20P %.16j %.8u %.10T %.10S %.10M %.10l %.6D %.5C %.10m %20R %80Z %k'
"
alias xpjobs="squeue -S 'P,t,-p,i' --format '%.20i %.16j %.8u %.10T %.10M %.10l %.5C %.10m %R'"
alias xinfo="xpinfo -u $USER"
alias xjobs="xpjobs -u $USER"
```


#### Step 2: Before using Eager, you will need to do the following steps

#### 1. [Tower.nf](https://cloud.tower.nf/)
- Creat account in [Tower.nf](https://cloud.tower.nf/)
- Click on your avatar
- Go to **Your Tokens**
- Press **Add Token**
- Copy the token
- Add it to your .bashrc script in Mjolnir
  ```
  nano ~/.bashrc
  ```
  *Replace the $YOUR_TOKEN to your copied token*
  ```
  export TOWER_ACCESS_TOKEN=$YOUR_TOKEN
  ```
- Save the script by pressing **control + X** and then **Y**
- Load it into environment
  ```
  source ~/.bashrc
  ```
#### 2. Temporary directory 
- Modify your .bashrc or .bash_profile script
    ```
    nano ~/.bashrc

    ```
- Copy and paste the following lines into your .bashrc or .bashprofile, remember to replace all the **USERID** to your KU-ID.
  ```
  export TMPDIR='/projects/course_sgbb20001/USERID/tmp'
  export TMP='/projects/course_sgbb20001/USERID/tmp'
  export TEMP='/projects/course_sgbb20001/USERID/tmp'
  export NXF_SINGULARITY_TMPDIR='/projects/course_sgbb20001/USERID/tmp'
  export NXF_OPTS='-Djava.io.tmpdir=/projects/course_sgbb20001/USERID/tmp -Xms1g -Xmx4g'
  export NXF_TEMP='/projects/course_sgbb20001/USERID/tmp'
  export SINGULARITY_LOCALCACHEDIR='/projects/course_sgbb20001/USERID/tmp'
  export SINGULARITY_TMPDIR='/projects/course_sgbb20001/USERID/tmp'
  export SINGULARITY_LOCALCACHEDIR='/projects/course_sgbb20001/USERID/tmp'
  export SINGULARITY_CACHEDIR='/projects/course_sgbb20001/USERID/tmp'
  export SINGULARITY_TMPDIR='/projects/course_sgbb20001/USERID/tmp'
  mkdir -p $TMPDIR
  ```
- Save the modified .bashrc or .bash_profile script by pressing **control X** and then **Y**. 
- Load the set-up environment
  ```
  source ~/.bashrc
  ### or
  source ~/.bash_profile
  ```
#### 3. Customize Eager config file
A customized Eager config file is provided on the main page. 
- Set running memory, CPUs and time.
  For example, we set up 8 CPUs and total 100 GB memory for running maximum 96 hours on bwa alignment 
  ```
  process {
    withName: bwa {
      time = '96h'
      cpus = 8
      memory = 100.GB
    }
  }
  ```
- Set the temporary directory

  In this setup, we designate our own /scratch/tmp as the working temporary folder instead of the default /tmp folder on the running node. This approach helps avoid filling up the /tmp folder, which can lead to job failures.
  ```
  singularity {
    enabled = true
    autoMounts = true
    runOptions = '-B $SINGULARITY_TMPDIR:/tmp -B $SINGULARITY_TMPDIR:/scratch'
    envWhitelist = ['SINGULARITY_TMPDIR']
  }
  ```
#### 4. Prepare Input.tsv file 
The input .tsv file mainly includes 11 different columns, as the following figure. <img width="1020" alt="Screenshot 2023-12-18 at 17 26 27" src="https://github.com/Schroeder-Group/Eager-Workshop/assets/54803911/e3c30d6a-8345-4650-bb13-bdbc24d54594">

The detail description on those columns would be checked in [Eager website page](https://nf-co.re/eager/2.5.0/docs/usage).

### Step 4: Running Eager on Mjolnir, Step by Step

1. Remember once you are in Mjolnir, in the folder where you want to work, transfer the tsv file

The raw data is located at "/projects/course_sgbb20001/people/hsf378/Raw/RawData.fastq.gz"

2. Open a screen session in the folder where you want to work, by writing in the terminal:
```
screen
```
3. You have to purge the modules you could have used before, by writing:
```
module purge
```
4. Then, load Nextflow version you have. Below an example. !!! I will add a line to check which one you have !!!
```
module load singularity/3.8.7 openjdk/11.0.0 nextflow/22.10.4
```
5. when you are running at your first time, you need to pull the eager:

```
nextflow pull nf-core/eager -r 2.4.6
```

6. create a folder and enter the folder
```
mkdir Eager
cd Eager
cp /projects/course_sgbb20001/people/hsf378/Eager/input.tsv ./
```
7. run Eager with the following code:

```
nextflow run nf-core/eager -r 2.4.6 -profile mjolnir_globe --input input.tsv --fasta /projects/course_sgbb20001/data/databases/Human/hs.build37.1.fasta --fasta_index /projects/course_sgbb20001/data/databases/Human/hs.build37.1.fasta.fai --bwa_index /projects/course_sgbb20001/data/databases/Human --seq_dict /projects/course_sgbb20001/data/databases/Human/hs.build37.1.dict --bam_unmapped_type fastq --bam_mapping_quality_threshold 30 --run_bedtools_coverage --run_bam_filtering -with-tower -c /projects/course_sgbb20001/people/hsf378/Eager_testing.config --complexity_filter_poly_g --anno_file /projects/course_sgbb20001/data/databases/Human/MT.bed -name eager_human 
```

8. To follow up the status, you can check in Nextflow tower

If, an issue happens, and you see written 'failed', you can rerun your commands by adding at the end of the previous command '-resume'. It will run where it stops. You will have also to change the name.
Please find below an example of the commands, if you need to rerun

```
nextflow run nf-core/eager -r 2.4.6 -profile mjolnir_globe --input input.tsv --fasta /projects/course_sgbb20001/data/databases/Human/hs.build37.1.fasta --fasta_index /projects/course_sgbb20001/data/databases/Human/hs.build37.1.fasta.fai --bwa_index /projects/course_sgbb20001/data/databases/Human --seq_dict /projects/course_sgbb20001/data/databases/Human/hs.build37.1.dict --bam_unmapped_type fastq --bam_mapping_quality_threshold 30 --run_bedtools_coverage --run_bam_filtering -with-tower -c /projects/course_sgbb20001/people/hsf378/Eager_testing.config --complexity_filter_poly_g -name YYY -resume
```
9. To close the screen you created before, by pressing **ctrl A** and then **D**; In order to re-enter the screen, you could do this:
    ```
    screen -r 
    ```
10. In Nextflow tower, you should be to see if your run succeeded or not.

If it succeeds, you should be able to see different statistics (e.g., fastqc, ect...).

You should safe the html files in your own folder. Then, you should open it with a browser, go to configure columns, choose all columns and copy/paste it to new google sheet. And keep it :)

You have now the main statistic for your data. Certain colomns are more important than others such as: .....

```
cd results/multiqc
```

open html file with a browser

