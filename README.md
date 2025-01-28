# GENOLOADER

## HOW TO ESTIMATE GENETIC LOAD IN INDIVIDUALS AND POPULATIONS FROM GENOMIC DATA

Genetic load is estimated by summing up the fitness deficit caused by all deleterious alleles present in the genome of an individual.

Rarely we have accurate estimates of the fitness effects of the segregating variants in out population when analysing non-model species, as it is often the case in Conservtion Biology. There are ways to estimate the distribution of fitness effects (DFE, see [here](https://academic.oup.com/mbe/article/41/5/msae070/7641109) for a very recent approach) but, in principle, we need to know the actual fitness effect, or the selection coefficient, of each of the deleterious variants in our population.

So our first question is:  

**What are the selection coefficients of all the variant positions in the target individuals/populations?**  

**How can get a good guess of the selection coefficient for the variants in our genomic dataset?**  

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

**What can you do to improve your estimates?**  

## SNPEff variants annotation

We can actually skip SNPEff actual annotation in this tutorial but you can see how it works below and check the output (the summary, an annotated vcf file)

### EXTRA 1: Running SNPEff and check the outputs

To annotate your genetic variants with SNPEff you a database for the species of interest. SNPEff database gets frequently updated so you are likely to find your species of interest. Otherwise, you can build your own database with a reference genome, an annotation of the genes, and, likely, a protein dataset. See below an example of this step:

##FRANCESCO/SEBA help me with this

Then, you can annotate the variants in your vcf file using SNPEff and the genomic database for your species of interest using a command line like this one 

#Annotating the vcf files

#Check what SNPEff added to your vcf. Understand the info added (Refer to this page: https://pcingola.github.io/SnpEff/snpeff/inputoutput/) 

#?#Question: ask something about the annotation like Why there could be more than one annotation? ADD HERE

## Summarizing genotype counts per variant effect categories

After getting SNPEff fitness effect predictions, we will parse the annotated vcf to summarize fitness impact, individual genotypes, type of variant, missing data, allele frequency per population, etc. MORE IMPORTANTLY: we need to correctly identify ancestral and derived allele at each of the segregating sites. You should have seen before the slides about polarization of ancestral-derived alleles.

We have put together a little python script called GENOLOADER that we can use to polarize ancestral derived alleles while parsing the vcf, summarizing the counts of the different genotypes per categories of SNPeff fitness impact, and making a plot.

Run the following command to see the options available to run GENOLOADER:

```
./GENOLOADER -h
```

GENOLOADER summarize variants in vcf files taking two ingroup populations (or one target population and one close outgroup) and one more distant outgroup population. We need to set the minimum number of individuals per population. In this example, we will actually set no missing data in the ingroup populations and require at least three individuals in the outgroup population.

**Why allowing for missing data in the ingroups while summarizing genotype counts to estimate genetic load is not a good idea?**  

Then, we will set the way GENOLOADER is polarizing ancestral/derived alleles (`-r`).The options are:  
POP_OUT: use the outgroup allele as ancestral  
POP1: use the major allele in population 1 as ancestral  
POP2: use the major allele in population 2 as ancestral  
POP1POP2: use the major allele in a joint population 1&2 sample as ancestral  

Finally, we need to tell GENOLOADER to counts the genotypes only using the sites we believe are best polarized given the option used for polarization (`-pol`).  The options are:
unfolded, in1Fold, in2Fold, unfoldOutMiss

Refer to the slides on polarization that we saw before to decide which one of these options is the best given the polarization method.

We will now do some tests on a sample of different penguin individuals. In the folder XXX, there is a vcf file called XXX, and three text files called pop1 (Emperor penguin), pop2 (King penguin), pop_out (Adelie + Gentoo penguin).

Check these files.

**How many individuals per populations do we have?**  

**How many variants are in the vcf?**  

Here is the command line to run GENOLOADER using the outgroup population to polarize ancestral/derived alleles:

```
./GENOLOADER -p1 pop1 -p2 pop2 -p0 pop_out -f penguins_scaf1.vcf -m1 24 -m2 24 -m0 3 -r POP_OUT -pol unfolded
```

**How many variants have been repolarized?**  

Now check the output: gt file, counts file and plots.

**How many HIGH impact variants can you see in the gt file? And how many MODIFIER ones?**  

**Which of the two penguins has the higher masked load and realized load**  

We will discuss the output more altogether.

Now change the command line to run again GENOLOADER:

1) using major allele in POP1 as ancestral;
   
2)  using major allele in POP2 as ancestral;

3) using major allele in both POP1 and POP2 as ancestral.

**Compare the output plots and get ready for discussion altogether**

**How polarization is impacting our estimates of masked and realized load?**

### EXTRA 2: More fish quiz to test for purging of deleterious variants

#Describe the case study and the material provided

#Polarize on the reference source population and follow the fate of segregating variants during a range expansion

#Count the genotypes per individual and plot per population.

#Normalize the counts of Homo_der and make a plot of normalized increase of homozygosity along an expansion route







