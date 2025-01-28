# GENOLOADER

This tutorial is for the Genomic Load activity at the **Workshop on Population and Speciation Genomics 2025**

First, log into your AWS instance and move to the folder activity `~/workshop_material/30_genetic_load_genotypes_counts`

You'll need to execute some scripts in the bash terminal and then visualize pdf files, either via guacamole or downloading the files with `scp`

## HOW TO ESTIMATE GENETIC LOAD IN INDIVIDUALS AND POPULATIONS FROM GENOMIC DATA

Genetic load is estimated by summing up the fitness deficit caused by all deleterious alleles present in the genome of an individual.

Rarely we have accurate estimates of the fitness effects of the segregating variants in out population when analysing non-model species, as it is often the case in Conservtion Biology. There are ways to estimate the distribution of fitness effects (DFE, see [here](https://academic.oup.com/mbe/article/41/5/msae070/7641109) for a very recent approach) but, in principle, we need to know the actual fitness effect, or the selection coefficient, of each of the deleterious variants in our population.

So our first question is:  

**What are the selection coefficients of the variants in our target individuals/populations?**  

**How can we get a good guess of the selection coefficient for the variants in our genomic dataset?**  

If any direct estimate of fitness effects is not feasible, we need to resort to model-based predictions. 

Two main approaches are currently popular: 

1) Predictions based on coding sequence annotation (relying on [SNPEff/SIFT](https://pcingola.github.io/SnpEff/) software):  
---It relies on the annotation of the genes in the reference genome assembly sequence we used to align and call the variant positions in the target individuals:  
---The rationale is that changes in the coding sequence can be ranked by their effect on the protein product (nonsynonyomus, synonymous, Loss-Of-Function variants, splice sites modifications, etc);  
---Main limitation: the inference is restricted to the coding sequence while non-coding variants are considered as neutral (not assessing the changes in e.g. promoters).

2) Predictions based on evolutionary conservation scores (eg. [GERP](https://doi.org/10.1371/journal.pcbi.1001025); [PhyloP](http://compgen.cshl.edu/phast/phyloP-tutorial.php)):  
---It relies on the alignment of the target reference genome assembly sequence with other genome assmblies from species which are more or less distant phylogenetically speaking;  
---The rationale is that the more conserved a site across evolutionary time the more intense the negative selection coefficient on not common alleles;  
---Main limitation: the inference is restricted to the portion of the genomes that can be accurately aligned across distant species. Computationally very intensive. Integration with General Language Model integration as in [GPN-MSA](https://www.nature.com/articles/s41587-024-02511-w) holds promise for big improvments.  
---Genome alignments are achieved with approaches like the one in the software [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), but check also the method outlined in [Dehasque et al 2024](https://www.cell.com/cell/fulltext/S0092-8674(24)00577-4).  

In this tutorial we will focus on the first approach using SNPEff to categorize the fitness impact of our target variants in HIGH, MODERATE, LOW and MODIFIER impact. Please check what each of these categories represent at this [page](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)

**How do you think the genome assembly and annotation will impact your prediction of the fitness effects of the target variants?**  

## SNPEff variants annotation

We can actually skip SNPEff actual annotation in this tutorial but you can see how it works below and check an annotated vcf file in the `~/workshop_material/30_genetic_load_genotypes_counts/penguins` directory

### EXTRA 1: Running SNPEff and check the outputs

To annotate your genetic variants with SNPEff you a database for the species of interest. SNPEff database gets frequently updated so you are likely to find your species of interest. Otherwise, you can build your own database with a reference genome, an annotation of the genes, and, likely, a protein dataset. See below an example of this step:

>>>>>>>>>>>>>>>>>>>
Move to the snpEff folder:

cd ~/path/to/snpEff/

Create a folder called data inside the snpEff folder:
mkdir data

Move inside the data folder:
cd data/

Create a folder called "genomes" inside your data folder:
mkdir genomes

Create a folder called as your species' scientific name (e.g. Pterois miles) inside the data folder:
mkdir Pterois miles

Now copy your species' reference genome to the "genomes" folder and call the file as your species' scientific name (e.g. Pterois_miles.fa):
cp ~/path/to/your/reference/Pterois_miles.fa ./genomes/

Now copy your species' genome annotation file into the folder called as your species' scientific name (e.g. Pterois_miles) with the following command:
cp ~/path/to/your/reference/genes.gff ./Pterois_miles/

Call the annotation file "genes.gff":
mv ./Pterois_miles/Annotation.gff ./Pterois_miles/genes.gff

Now move back to the snpEff main folder: 
cd ../

And edit the file "snpEff.config" adding the following entries
"Pterois_miles.genome : Pterois_miles"
"Pterois_miles.reference : Pterois_miles"
with the following command:
awk '1; END { print "Pterois_miles.genome : Pterois_miles" }' snpEff.config > temp && mv temp snpEff.config
awk '1; END { print "Pterois_miles.reference : Pterois_miles" }' snpEff.config > temp && mv temp snpEff.config

Now run the following command to build the dataset

java -Xmx4g -jar ~/path/to/snpEff/snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v Pterois_miles

Then, you can annotate the variants in your vcf file using SNPEff and the genomic database for your species of interest using a command line like this one 

java -Xmx4g -jar ~/path/to/snpEff/snpEff.jar -v Pterois_miles YOUR_VCF_FILE.vcf > ~/path/to/outputfolder/Annotated.vcf


## Summarizing genotype counts per variant effect categories

Before starting and as we skipped the actual annotation with SNPEff, answer the following questions to get acquainted with SNPEff annotated files.  

If you have not done it yet, move to this directory:  
`~/workshop_material/30_genetic_load_genotypes_counts/penguins`

Check SNPEff annotation details at this [page](https://pcingola.github.io/SnpEff/snpeff/inputoutput/) 

**What is the predicted effect of the variant at position 3978434 of the `penguins_scaf1.vcf`? Why it is predicted to be so deleterious**  

**Why there could be more than one annotation?**

After getting SNPEff fitness effect predictions, we will parse the annotated vcf to summarize fitness impact, individual genotypes, type of variant, missing data, allele frequency per population, etc. MORE IMPORTANTLY: we need to correctly identify ancestral and derived allele at each of the segregating sites. You should have seen before the slides about polarization of ancestral-derived alleles.

We have put together a little python script called GENOLOADER that we can use to polarize ancestral derived alleles while parsing the vcf, summarizing the counts of the different genotypes per categories of SNPeff fitness impact, and making a plot.

Run the following command to see the options available to run GENOLOADER:

```
./GENOLOADER -h
```

GENOLOADER summarize variants in vcf files taking two ingroup populations (or one target population and one close outgroup) and one more distant outgroup population. We need to set the minimum number of individuals per population. In this example, we will actually set no missing data in the ingroup populations and require at least three individuals in the outgroup population.

**Why allowing for missing data in the ingroups while summarizing genotype counts to estimate genetic load is not a good idea?**  

Then, we will set the way GENOLOADER is polarizing ancestral/derived alleles (`-r`).   
The options are:  

POP_OUT: use the outgroup allele as ancestral  

POP1: use the major allele in population 1 as ancestral  

POP2: use the major allele in population 2 as ancestral  

POP1POP2: use the major allele in a joint population 1&2 sample as ancestral  

Finally, we need to tell GENOLOADER to counts the genotypes only using the sites we believe are best polarized given the option used for polarization (`-pol`).  The options are:
unfolded, in1Fold, in2Fold, unfoldOutMiss

Refer to the slides on polarization that we saw before to decide which one of these options is the best given the polarization method.

We will now do some tests on a sample of different penguin individuals. In the folder XXX, there is a vcf file called XXX, and three text files called pop1 (Emperor penguin), pop2 (King penguin), pop_out (Adelie + Gentoo penguin).

Check these files.

**How many individuals per populations do you see?**  

**How many variants are in total in the vcf?**  

Here is the command line to run GENOLOADER using the outgroup population to polarize ancestral/derived alleles:

```
./GENOLOADER -p1 pop1 -p2 pop2 -p0 pop_out -f penguins_scaf1.vcf -m1 24 -m2 24 -m0 3 -r POP_OUT -pol unfolded
```

**How many variants have been repolarized?**  

Now check the output files: gt file, counts file and plots.

**How many HIGH impact variants can you see in the gt file? And how many MODIFIER ones?**  

**Which of the two penguins has the higher masked load and realized load**  

We will discuss the results further during the final wrap-up.

Now change the command line above to run GENOLOADER again but:

1) using major allele in POP1 as ancestral;
   
2) using major allele in POP2 as ancestral;

3) using major allele in both POP1 and POP2 as ancestral.
 

**Compare the output plots and get ready for discussion altogether**


**How polarization is impacting our estimates of masked and realized load?**


We will discuss the results further during the final wrap-up.


### EXTRA 2: More fish quiz to test for purging of deleterious variants

We are back to the same invasive fish populations as in the previous activity (ROH estimates). We sampled three populations from the invasive range as well as one population from the source range. If you have completed the previous activity, you should have seen the effect of the invasion bottleneck on the genetic diversity in terms of ROHs. Let's now check what happened those variants with medium-high deleteriousness.

For more fun, we shuffled the labels of the populations, so that you do not know which one is the source range. But we already ran the polarization using the major allele in the source population as ancestral. Our aim here is to track segregating sites during the highetened drift experienced during the invasion and the following range expasion.

Our main questions are:

**What is the fate of medium or high segregating variants during a range expansion?**

**Are we able to detect different trajectories which we can identify as puryfing selection relaxation or, on the contrary, with purging?**

First run the counting and plotting steps of the modified GENOLOADER script, called GENOLOADER_countsPlots, using the file `fish.complete.gt` and
assigning the popA to popD to the pop1 to pop4 options in the script

```
./GENOLOADER_countsPlots -f fish.gt -p1 PopA -p2 PopB -p3 PopC -p4 PopD -pol in1Fold
```

**Which one is the source population?**  

**Can you say what is the likely direction of the expansion?**  

**How are the homozygous genotypes counts for the MD or HI changing along the expansion?**

Try to run the following script to plot individual HOMO_DER counts normalized by the average counts in the source population.

```
./normCountsPlotting.py -f fish.counts -ref {choose one among A, B, C, D}
```

**What is this plot showing?**  

**Is there any signal of purging?**

We will discuss the results further during the final wrap-up.
