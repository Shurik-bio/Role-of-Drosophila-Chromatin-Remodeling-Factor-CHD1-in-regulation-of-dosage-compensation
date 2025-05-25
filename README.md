# Role-of-Drosophila-Chromatin-Remodeling-Factor-CHD1-in-regulation-of-dosage-compensation

**Students:**
- Egor Kulikov, 
Head of the Laboratory of Applied Genetics and Molecular Diagnostics, FSC ARRTPI

- Aleksandra Mikhailova 
Scientific Researcher, Lab of Postgenomic Investigation, VIR

**Supervisor:**

- Aleksandr Konev
PhD, Senior Researcher, Group of Genetic Studies Chromatin and Reparation,
Petersburg Nuclear Physics Institute named by B.P.Konstantinov of NRC «Kurchatov Institute»


### Introduction

Gene dosage compensation (DC) is an epigenetic mechanism that allows equalization of the expression level of sex-linked genes in males and females of species in which sex is determined by sex chromosomes [1]. In the fruit fly *Drosophila melanogaster*, DC is achieved by doubly increasing the level of expression of genes on the single X chromosome in males. It was reported earlier that X chromosome in the *Chd1* mutant males becomes shortened, decondensed and thickened, while in females it characterizes a normal structure [2]. In this case, the *Chd1* gene product of maternal origin is concentrated in the X chromosome, causing its specific staining. Thus discovered the role of CHD1 in regulating the structure of the X chromosome of males, associated with the phenomenon of DC. The CHD1 protein is considered to be one of the key factors in the replication-independent assembly of chromatin containing variant histone H3.3. In *Drosophila* absence of CHD1 protein in embryos caused misassembly of H3.3 in paternal pronuclear chromatin, while loss of CHD1 in the adult brain resulted in decreased H3.3 incorporation into chromatin, global chromatin disruption, transcriptional dysregulation, and metabolic reprogramming [3].

In this study we aimed on estimation of potential influence of the CHD1 factor on gene expression profiles in *Drosophila*. 

### Aim and Tasks

Study the specific roles of the *Drosophila* CHD1 in dosage compensation

The tasks include:
- Analyze the raw RNA-sequenced data:
  - Perform quality control
  - Align reads on rRNA reference genome sequence with Bowtie2
  - Align reads without rRNA on reference genome with STAR
  - Gene reading counting with FeatureCounts

- Visualize and analyze mapped reads in RStudio: 
  - Find differentially expressed genes using DESeq2
  - Compare different groups: male-female and mutant-control using MA-plots 
  - Reveal statistically overrepresented GO-terms by the GO  enrichment analysis in FlyBase.org database
  - Study the role of CHD1 in dosage compensation in Drosophila

### Data

The head of four-day-old imago of fruit fly *Drosophila melanogaster* in two technical replicates was used for RNA isolation. The wild type is represented by the Oregon R line, while the mutants are characterized by depletion in the *Chd1* gene. Samples are presented in two biological replicates. To prepare the libraries, we used the MGIEasy RNA Directional Library Prep Set kit.ver. 2.1 reagent kit. These libraries were sequenced using an MGI platform.

### Workflow

Our workflow includes following steps:

![workflow](Figures/Workflow.png)

### Technical notes




### Results and Discussion

[Текст ссылки](https://www.example.com)

### References

1. Shevelyov YY, Ulianov SV, Gelfand MS, Belyakin SN, Razin SV. Dosage Compensation in Drosophila: Its Canonical and Non-Canonical Mechanisms. Int J Mol Sci. 2022 Sep 19;23(18):10976. doi: 10.3390/ijms231810976. PMID: 36142884; PMCID: PMC9506574.
2. Konev A.Y., Tiutiunnik A.A., Baranovskaya I.L. The influence of the Chd1 chromatin assembly and remodeling factor mutations on Dro-sophila polytene chromosome organization. Citologija. 2016;58(4):281-284. (in Russian)
3. Schoberleitner I, Bauer I, Huang A, Andreyeva EN, Sebald J, Pascher K, Rieder D, Brunner M, Podhraski V, Oemer G, Cázarez-García D, Rieder L, Keller MA, Winkler R, Fyodorov DV, Lusser A. CHD1 controls H3.3 incorporation in adult brain chromatin to maintain metabolic homeostasis and normal lifespan. Cell Rep. 2021 Oct 5;37(1):109769. doi: 10.1016/j.celrep.2021.109769


