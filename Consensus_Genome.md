# Consensus Genome reconstruction
First, we will create a folder called consensus
```
mkdir -p data/consensus
cd data/consensus
```
Then, we create a interactive node to run the following command:
```
screen -S Day4
srun --ntasks=1 --cpus-per-task=3 --mem=30G -t 8:00:00 --reservation=Course_SGBB20001 --account=teaching --pty bash
```
## Extract Reads Mapped to Target Species
To extract reads mapped to a specific target species, we will use `samtools`. The target species is identified using a known accession ID. This process filters the reads from the BAM file that align to the target species.

### Command:
**Remember to replace the accession_id to your target genome id**
```
module load samtools bedtools 
samtools view -b ../filterbam/Lola.dedup.filtered.bam "accession_id" > target_species.bam
samtools index target_species.bam
```

Question1: How many reads in your target_species.bam file?
```
samtools view -c target_species.bam
```
Question2: What is the proportion of genome was covered by those reads?
**Hint: check your filterbam output file**

Question3: How was the coverage evenness score?
**Hint: check your filterbam output file**

### Visualization
#### Genome coverage plot

**Depth of coverage** measures the average number of reads covering each base in a specific region.

$$ Depth\ of\ Coverage = {Total\ Bases\ Covered\ by\ Reads\over Region\ Length} $$

Total Bases Covered by Reads: The sum of the coverage values across all positions in the region.
Region Length: The total number of bases in the region.

**Breadth of coverage** measures the proportion of a region that is covered by at least one read.

$$ Breadth\ of\ coverage = {Number\ of\ Covered\ Bases\over Regio\ Length} $$

Number of Covered Bases: The number of bases in the region with at least one read mapped.

1. calculate the depth of coverage at each position, and custom the format of the depth file.
```
samtools depth -a target_species.bam > depth.txt
awk '{print $1, $2, $2, $3}' OFS="\t" depth.txt > depth_as_bed.txt
```
2. calculate the average depth of coverage in 20000 base pairs window size
**replace the accession id to your target species id **
```
samtools view -H target_species.bam | grep "accession_id" |awk '$1 == "@SQ" {print substr($2, 4) "\t" substr($3, 4)}' > genome.txt
bedtools makewindows -g genome.txt -w 20000 > windows.bed
bedtools map -a windows.bed -b depth_as_bed.txt -c 4 -o mean > window_depth.bed
```
3. redo the previous steps, but we filtered out reads with mapping quality lower than 20
```
samtools view -q 20 target_species.bam -Sb > target_species.mq20.bam
```
Question: How many reads were filtered out?
```
samtools view -c target_species.mq20.bam
```
```
samtools depth -a target_species.mq20.bam > depth_mq20.txt
awk '{print $1, $2, $2, $3}' OFS="\t" depth_mq20.txt > depth_as_bed_mq20.txt
bedtools map -a windows.bed -b depth_as_bed_mq20.txt -c 4 -o mean > window_depth_mq20.bed
```
4. visulize the depth of coverage plot
```
module load conda
conda activate /projects/course_sgbb20001/data/envs/Pathopipe
cp /projects/course_sgbb20001/people/hsf378/scripts/coverage_plot.R ./
Rscript coverage_plot.R
```

#### Damage plot
Take a screen shot or save the pdf file for the metaDMG results.

## Variant calling
1. Download the reference genome from NCBI
```
module load datasets
datasets download genome accession XXXX
cp /projects/course_sgbb20001/people/hsf378/genomes/ncbi_dataset/data/GCA_003665785.1* ./
cp /projects/course_sgbb20001/people/hsf378/consensus/GCA_000195855.1.fasta ./
```
**Replace the XXX to your target species genome accession id**

2. Calling the variants
```
module load gsl/2.5 perl bcftools
samtools faidx target_species_genome.fna
bcftools mpileup -f target_species_genome.fna target_species.bam | bcftools call -m -Ov -o Lola.variant.vcf --ploidy 1
bgzip Lola.variant.vcf
tabix Lola.variant.vcf.gz
```
**Replace the target_species_genome.fna to the genome you downloaded before**
Question 1: How many variants were called?
```
grep -v "^#" Lola.variant.vcf | wc -l 
```
Question 2: How many SNPs were called?
```
bcftools view -v snps Lola.variant.vcf  | grep -vc "^#" 
```
## Reconstruct consensus genome 
```
bcftools consensus -f target_species_genome.fna Lola.variant.vcf.gz > target_species.consensus.fa
```
**Replace the target_species_genome.fna to the genome you downloaded before**

# Phylogenetic inference

To reconstruct the evolutionary relationships among species, genes, or other biological entities, Phylogenetic inference analysis provides a framework to understand the shared ancestry and divergence of species/genes/traits. 

Commonly used tools for phylogenetic tree building include MEGA, RAxML, BEAST, IQ-TREE, MrBayes, and Phylip, among others. These tools employ different theories and algorithms for tree construction, including:

Distance-based methods: Tools like MEGA use methods such as Neighbor-Joining (NJ) and Unweighted Pair Group Method with Arithmetic Mean (UPGMA). These approaches compute pairwise genetic distances and cluster taxa based on minimizing the overall tree length.

Maximum Likelihood (ML): Tools like RAxML and IQ-TREE use ML methods that evaluate the likelihood of a tree given a specific model of sequence evolution. The tree that maximizes this likelihood is selected as the best representation of evolutionary relationships.

Bayesian Inference (BI): Tools like MrBayes and BEAST apply Bayesian frameworks to estimate tree topologies and branch lengths, integrating prior knowledge and using Markov Chain Monte Carlo (MCMC) sampling to approximate posterior probabilities of trees.

Parsimony-based methods: Tools like PAUP* use the principle of parsimony, aiming to find the tree that minimizes the total number of evolutionary changes.

## Building the tree based on SNPs data and ML method

```
mkdir -p ~/data/Phylogeny
cd ~/data/Phylogeny
cp /projects/course_sgbb20001/people/hsf378/Phylogeny/SNPs_306_lola.leprae.fa ./
cp /projects/course_sgbb20001/people/hsf378/Phylogeny/iqtree.sh ./
module load iqtree
sbatch iqtree.sh
```
