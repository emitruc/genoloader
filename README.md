# GENOLOADER

This tutorial is for the Genomic Load activity at the **Workshop on Population and Speciation Genomics 2025**

First, log into your AWS instance and move to the folder activity `~/workshop_materials/30_genetic_load/genotype_counts`

Activate the conda environment

```
conda activate evomics_ml_workshop
```

You'll need to execute some scripts in the bash terminal and then visualize pdf files, either via guacamole or downloading the files with `scp`

## HOW TO ESTIMATE GENETIC LOAD IN INDIVIDUALS AND POPULATIONS FROM GENOMIC DATA

Genetic load is estimated by summing up the fitness deficit caused by all deleterious alleles present in the genome of an individual.

Rarely we have accurate estimates of the fitness effects of the segregating variants in out population when analysing non-model species, as it is often the case in conservation biology. There are ways to estimate the distribution of fitness effects (DFE, see [here](https://academic.oup.com/mbe/article/41/5/msae070/7641109) for a very recent approach) but, in principle, we need to know the actual fitness effect, or the selection coefficient, of each deleterious variant in our population.

So our first question is:  

**What are the selection coefficients of the variants in our target individuals/populations?**  

**How can we get a good guess of the selection coefficient for the variants in our genomic dataset?**  

If any direct estimate of fitness effects is not feasible, we need to resort to model-based predictions. 

The two most common approaches are: 

1) Predictions based on coding sequence annotation (relying on [SNPEff/SIFT](https://pcingola.github.io/SnpEff/) software):  
---It relies on the annotation of the genes in the reference genome assembly sequence it was used to align and call the variant positions in the target populations/individuals:  
---The rationale is that changes in the coding sequence can be ranked by their effect on the protein product (nonsynonyomus, synonymous, Loss-Of-Function variants, splice sites modifications, etc);  
---Main limitation: the inference is restricted to the coding sequence while non-coding variants are considered as neutral (not assessing the changes in e.g. promoters).

2) Predictions based on evolutionary conservation scores (eg. [GERP](https://doi.org/10.1371/journal.pcbi.1001025); [PhyloP](http://compgen.cshl.edu/phast/phyloP-tutorial.php)):  
---It relies on the alignment of the target reference genome assembly sequence with other genome assmblies from species which are more or less phylogenetically distant;  
---The rationale is that sites which are more conserved across evolutionary time are under stronger selection so that any change would have more negative impact on fitness;  
---Main limitation: the inference is restricted to the portion of the genomes that can be accurately aligned across distant species. Computationally very intensive. Integration with General Language Model, as in [GPN-MSA](https://www.nature.com/articles/s41587-024-02511-w), holds promise for big improvments.  
---Genome alignments can be done with software like [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), but check also the method outlined in [Dehasque et al 2024](https://www.cell.com/cell/fulltext/S0092-8674(24)00577-4).  

In this tutorial, we will focus on the first approach using SNPEff to categorize the fitness impact of a set of variants in HIGH, MODERATE, LOW and MODIFIER impact. Please check what each of these categories represent at this [page](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)

**How do you think the genome assembly and annotation will impact your prediction of the fitness effects of the target variants?**  

## SNPEff variants annotation

We will actually skip SNPEff actual annotation here but you can see how it works below and check an annotated vcf file in the `~/workshop_material/30_genetic_load/genotypes_counts/penguins` directory

### EXTRA 1: Building a database and annotating a vcf with SNPEff --- SKIP THIS STEP

To annotate a set genetic variants with SNPEff, the first step is to build a database for the species of interest. SNPEff databases get frequently updated so you are likely to find your species of interest. Otherwise, you can build your own database with a reference genome, an annotation of the genes, and, likely, a protein dataset. See below an example of this step:

Inside the SNPEff installation folder: create a folder called `data and, inside this a folder called `genomes` and a folder called `yourSpecies`.  

Copy in `genomes` folder  your species' reference genome fasta file giving its scientifc name (e.g., Pterois_miles.fa)

Copy your species' genome annotation file (`genes.gff`) into `yourSpecies` folder 

Look for the `snpEff.config` in the SNPEff installation folder and add at the end:

`Your_species.genome : Your_species`
`Your_species.reference : Your_species`

Here is an `awk` command that can help:
```
awk '1; END { print "Your_species.genome : Your_species" }' snpEff.config > temp && mv temp snpEff.config
awk '1; END { print "Your_species.reference : Your_species" }' snpEff.config > temp && mv temp snpEff.config
```

Then build the dataset:
```
java -Xmx4g -jar /path/to/snpEff/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v Your_species
```

Now, you can annotate the variants in your vcf file using SNPEff and the genomic database for your species of interest as follows:
```
java -Xmx4g -jar /path/to/snpEff/snpEff.jar -v Your_speciess YOUR_VCF_FILE.vcf > /path/to/outputfolder/Annotated.vcf
```

## Summarizing genotype counts per variant effect categories

As we skipped the actual annotation with SNPEff, answer the following questions to get acquainted with SNPEff annotated files.  

If you have not done it yet, move to the directory for this part of the activity:  

`~/workshop_materials/30_genetic_load/genotypes_counts/penguins`

Check SNPEff annotation details at this [page](https://pcingola.github.io/SnpEff/snpeff/inputoutput/) 

**What is the predicted effect of the variant at position 3978434 of the `penguins_scaf1.vcf`? Why is it predicted to be so deleterious**  

**Why there are sometimes more than one annotation per variant?**

We will now parse the annotated vcf to summarize fitness impacts, individual genotypes and types of variant while taking into account missing data.  
MORE IMPORTANTLY: we need to correctly polarize each variant, identifying ancestral and derived alleles. Remember about the slides we just saw on polarization of ancestral-derived alleles.

We have put together a little python script called GENOLOADER that can polarize ancestral derived alleles while parsing the vcf, summarizing the counts of the different genotypes per category of SNPeff fitness impact, and make a plot of the counts per genotype (homozygote for the ancestral allele, heterozygote,homozygote for the derived allele).

Type the following command to see the options available while running GENOLOADER:

```
./GENOLOADER -h
```

GENOLOADER summarize variants in vcf files taking two ingroup populations (or one target population and one close outgroup) and one more distant outgroup. We need to set the minimum number of individuals per each population/outgroup. In this example, we will actually set no missing data in the ingroup populations and require at least three individuals in the outgroup population.

**Why allowing for missing data in the ingroups while summarizing genotype counts to estimate genetic load is not a good idea?**  

We then choose how GENOLOADER polarizes ancestral/derived alleles with the option `-r`: 

POP_OUT: use the outgroup allele as ancestral  

POP1: use the major allele in population 1 as ancestral  

POP2: use the major allele in population 2 as ancestral  

POP1POP2: use the major allele in a joint population 1&2 sample as ancestral  

Finally, we need to tell GENOLOADER to counts the genotypes only using the sites we believe are best polarized given the option used for polarization using the option `-pol`:  
unfolded, in1Fold, in2Fold, unfoldOutMiss.

Refer to the slides on polarization to decide which one of these flags is the most suitable given the polarization method.

We will now do some tests on a sample of penguinss. In the folder `penguins`, there is a vcf file called `penguins_scaf1.vcf`, and three text files called `pop1` (Emperor penguin), `pop2` (King penguin), `pop_out` (Adelie + Gentoo penguin).

Check these files.

**How many individuals per populations do you see?**  

**How many variants are in total in the vcf?**  

Here is the command to run GENOLOADER using the outgroup population to polarize ancestral/derived alleles:

```
./GENOLOADER -p1 pop1 -p2 pop2 -p0 pop_out -f penguins_scaf1.vcf -m1 24 -m2 24 -m0 3 -r POP_OUT -pol unfolded
```

**How many variants have been repolarized (inverted)?**  

Now check the output files: a file with the same name as the vcf but with the .gt extension, one with the .counts extension and a pdf file.

**How many HIGH impact variants can you count in the gt file? And how many MODIFIER ones?**  

**Which one between the Emperor and the King penguins has the higher masked load and realized load**  

We will discuss the results further during the final wrap-up.  


The step above will take about 15 minutes on the AWS instance. Launch it, then open another terminal and navigate to the same folder `/home/wpsg/workshop_materials/30_genetic_load/genotype_counts/penguins`.  

Here make a subset of 100000 lines of the `penguins_scaf1.vcf` with this command:  
```
head -n 100000 penguins_scaf1.vcf > test.vcf
```

First, re-run the command line above with this reduced vcf to make it faster.  
Then, change the command line above to run GENOLOADER again with the same test.vcf file but:

1) use major allele in POP1 as ancestral;
   
2) use major allele in POP2 as ancestral;

3) use major allele in both POP1 and POP2 as ancestral.

Make sure you pair each of this option with the right polarization flag.


**Compare the output plots (full vcf vs. reduced vcf, and among different polarization options) and get ready for discussion altogether**

**How polarization is impacting our estimates of masked and realized load?**


### EXTRA 2: More fish quiz to test for purging of deleterious variants

We are back to the same invasive fish populations as in the previous activity (ROH estimates). We sampled three populations from the invasive range as well as one population from the source range. If you have completed the previous activity, you should have seen the effect of the invasion bottleneck on the genetic diversity in terms of ROHs. Let's now check what happened to the genetic load during the invasion.

Move to the directory called `~/workshop_material/30_genetic_load/genotype_counts/fish/`

For more fun, we shuffled the labels of the populations, so that you do not know which one is from the source range. We already ran the polarization using the major allele in the source population as ancestral. Our aim here is to track segregating sites during the highetened drift experienced during the invasion and the following range expasion.

**Do you agree with this choice?**  

Our main questions are:

**What is the fate of medium or high segregating variants during a range expansion?**

**Can we see the effects of purging or of relaxation of purifying selection?**

First run the counting and plotting steps of the modified GENOLOADER script, now called GENOLOADER_countsPlots, using the file `fish.gt` and
assigning the `popA` to `popD` files to the pop1 to pop4 options in the script.

```
./GENOLOADER_countsPlots -f fish.gt -p1 PopA -p2 PopB -p3 PopC -p4 PopD -pol in1Fold
```

**Which one is the source population?**  

**Can you say what is the likely direction of the expansion?**  

**How are the homozygous genotypes counts for the MD or HI changing along the expansion?**

Try to run the following script to plot individual homozygote derived counts normalized by the average counts in the source population.

```
./normCountsPlotting.py -f fish.counts -ref {choose one among A, B, C, D}
```

**What is this plot showing?**  

**Is there any signal of purging?**



We will discuss the results further during the final wrap-up.
