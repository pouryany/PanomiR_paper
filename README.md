## Introduction

This document contains the companion code to the paper: “PanomiR:
inferring miRNA targeting of coordinated multi-pathway disease
transcriptional programs.”

For information about PanomiR R-package please refer to its [Github
repository](https://github.com/pouryany/PanomiR) or the [Bioconductor
directory of PanomiR](https://bioconductor.org/packages/PanomiR)

For questions, comments, and other queries, contact pnaderiy \[at\]
bidmc \[dot\] harvard \[dot\] edu or whide \[at\] bidmc \[dot\] harvard
\[dot\] edu

## Installation

PanomiR can be accessed via Bioconductor. To install, start R (version
\> 4.2.0) and run the following code.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("PanomiR")
```

You can also install the latest development version of PanomiR using
GitHub.

``` r
devtools::install_github("pouryany/PanomiR")
```

## Overview

PanomiR is a pipeline to prioritize disease-associated miRNAs based on
activity of disease-associated pathways. The input datasets for PanomiR
are (a) a gene expression disease dataset along with covariates, (b) a
background collection of pathways/genesets, and (c) a collection of
miRNAs containing gene targets.

The general workflow of PanomiR is (a) generation of pathway summary
statistics from gene expression data, (b) detection of differentially
activated pathways, (c) finding coherent groups, or clusters, of
differentially activated pathways, and (d) detecting miRNAs targeting
each group of pathways.

Individual steps of the workflow can be used in isolation to carry out
different analyses. The following sections outline each step and
material needed to execute PanomiR.

## 1. Pathway summarization

PanomiR can generate pathway activity profiles given a gene expression
dataset and a list of pathways.Pathway summaries are numbers that
represent the overall activity of genes that belong to each pathway.
These numbers are calculated based on a methodology previously described
in part in (Altschuler et al. 2013; Joachim et al. 2018). Briefly, genes
in each sample are ranked by their expression values and then pathway
summaries are calculated as the average rank-squared of genes within a
pathway. The summaries are then center and scaled (zNormalized) across
samples.

The default list of background pathways in PanomiR is formatted into a
table (`data("path_gene_table")`). The table is based on canonical
pathways collection of Molecular Signatures Database (MSigDB) V6.2 and
it contains annotated pathways from a variety of sources (Liberzon et
al. 2011).

\*\* Users interested in using other pathway/geneset backgrounds, such
as newer versions of MSigDB or KEGG, should refer to the
[appendix](#geneset) of this manual.

This section uses a reduced example dataset from The Cancer Genome Atlas
(TCGA) Liver Hepatocellular Carcinoma (LIHC) dataset to generate Pathway
summary statistics (Ally et al. 2017). **Note:** Make sure that you
select gene representation type that matches the rownames of your
expression data. The type can be modified using the `id` argument in the
function below. The default value for this argument is `ENSEMBL`.

``` r
library(PanomiR)

# Pathway reference from the PanomiR package
data("path_gene_table")
data("miniTestsPanomiR")
# Generating pathway summary statistics 

summaries <- pathwaySummary(miniTestsPanomiR$mini_LIHC_Exp,
                            path_gene_table, method = "x2",
                            zNormalize = TRUE, id = "ENSEMBL")

head(summaries)[,1:2]
```

## 2. Differential Pathway activation

Once you generate the pathway activity profiles, as discussed in the
last section, there are several analysis that you can perform. We have
bundled some of the most important ones into standalone functions. Here,
we describe differential pathway activation profiling, which is
examining differences in pathway activity profiles in user-determined
conditions.

At this stage you need to provide a pathway-gene association table, an
expression dataset, and a covariates table. You need to specity what
covariates you would like to contrast. You also need to provide a
contrast, as formatted in limma. If the contrast is not provided, the
function assumes the first two levels of the provided contrast
covariate. **Note:** make sure the contrast covariate is formatted as
factor.

``` r
output0 <- differentialPathwayAnalysis(
                        geneCounts = miniTestsPanomiR$mini_LIHC_Exp,
                        pathways =  path_gene_table,
                        covariates = miniTestsPanomiR$mini_LIHC_Cov,
                        condition = 'shortLetterCode')

de.paths <- output0$DEP

head(de.paths,3)
```

## 3. Finding clusters of pathways

PanomiR provides a function to find groups coordinated differentially
activated pathways based on a pathway co-expression network (PCxN)
previously described in (Pita-Juárez et al. 2018). Briefly, PCxN is a
network where nodes are pathways and links are co-expression between the
nodes. It is formatted into a table were rows represent edges. The edges
of PCxN are marked by two numbers, 1- a correlation co-efficient and 2-
a significance adjusted p-value. Cut-offs for both of these numbers can
be manually set using PanomiR functions. See function manuals for more
info.

PCxN and its associated genesets are already released and can be
accessed through following Bioconductor packages:
[pcxn](http://bioconductor.org/packages/release/bioc/html/pcxn.html) and
[pcxnData](http://bioconductor.org/packages/release/data/experiment/html/pcxnData.html).

Here we have provided a small version of PCxN for tutorial purposes. A
more recent version of PCxN based on MSigDB V6.2 is available through
the data repository accompanying PanomiR manuscript, which can be found
[here](https://github.com/pouryany/PanomiR_paper).

``` r
# using an updated version of pcxn 

set.seed(2)
pathwayClustsLIHC <- mappingPathwaysClusters(
                            pcxn = miniTestsPanomiR$miniPCXN, 
                            dePathways = de.paths[1:300,],
                            topPathways = 200,
                            outDir=".",
                            plot = FALSE,
                            subplot = FALSE,
                            prefix='',
                            clusteringFunction = "cluster_louvain",
                            correlationCutOff = 0.1)


head(pathwayClustsLIHC$Clustering)
```

## 4. Prioritizing miRNAs per cluster of pathways.

PanomiR identifies miRNAs that target clusters of pathways, as defined
in the last section. In order to this, you would need a reference table
of miRNA-Pathway association score (enrichment). We recommend using a
customized miRNA-Pathway association table, tailored to your
experimental data. This section provides an overview of prioritization
process. Readers interested in knowing more about the technical details
of PanomiR are refered to accompaniying publication (Work under
preparation).

### Enrichment reference

Here, we provide a preprocessed small example table of miRNA-pathway
enrichment in `miniTestsPanomiR$miniEnrich` object. This table contains
enrichment analysis results using Fisher’s Exact Test between MSigDB
pathways and TargetScan miRNA targets. The individual components are
accessible via `data(msigdb_c2)` and `data(targetScan_03)` (Agarwal et
al. 2015; Liberzon et al. 2011). This example table is contains only a
full subset of the full pairwise enrichment. You can refer to [section
5](#geneset) of this manual on how to create full tables and how to
customize them to your specific gene expression data.

### Generating targeting scores

PanomiR generates a score for individual miRNAs targeting a group of
pathways. These scores are generated based on the reference enrichment
table. We are interested in knowing to what extent each miRNA targets
pathway clusters identified in the last step (see previous section).
PanomiR constructs a null distribution of this targeting score for each
miRNA. The significance of observed scores from a given group of
pathways (clusters in this case) is contrasted against the null
distribution to generate a targeting p-value. These p-values are used to
rank miRNAs per cluster.

### Sampling parameter

The above described process requires repeated sampling to empirically
obtain the null distribution. The argument `sampRate` denotes the number
of repeats in the process. Note that in the example below, we use a
sampling rate of 50, the recommended rate is between 500-1000. Also, we
set the saveSampling argument to FALSE. This argument, if set TRUE,
ensures that the null distribution is obtain only once. This argument
should be set to TRUE if you wish to save your sampling and check for
different outputs from the clustering algorithms or pathway thresholds.

``` r
set.seed(1)
output2 <- prioritizeMicroRNA(enriches0 = miniTestsPanomiR$miniEnrich,
                    pathClust = miniTestsPanomiR$miniPathClusts$Clustering,
                    topClust = 1,
                    sampRate = 50, 
                    method = c("aggInv"),
                    outDir = "Output/",
                    dataDir = "outData/",
                    saveSampling = FALSE,
                    runJackKnife = FALSE,
                    numCores = 1,
                    prefix = "outmiR",
                    saveCSV = FALSE)

head(output2$Cluster1)
```

## 5. miRNA-Pathway enrichment tables

PanomiR best performs on tissue/experiment-customized datasets. In order
to do this, you need to create a customized enrichment table. You can
simply do so by using the pathway and miRNA list that we have provided
as a part of the package. simply, plug in the name of the genes present
(expressed) in your experiment in the following code

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")


customeTableEnrich <- miRNAPathwayEnrichment(mirSets = targetScan_03,
                                              pathwaySets = msigdb_c2,
                                              geneSelection = yourGenes,
                                              mirSelection = yourMicroRNAs,
                                              fromID = "ENSEMBL",
                                              toID = "ENTREZID",
                                              minPathSize = 9,
                                              numCores = 1,
                                              outDir = ".",
                                              saveOutName = NULL)
```

In the above section, the field `fromID` denotes the gene representation
format of your input list. Here is a quick example that runs fast. Note
that the `miRNAPathwayEnrichment()` function creates a detailed report
with parameters that are used internally. To get a smaller table that is
suitable for publication purposes, use `reportEnrichment()` function.

``` r
# using an updated version of pcxn 
data("msigdb_c2")
data("targetScan_03")

tempEnrich <-miRNAPathwayEnrichment(targetScan_03[1:30],msigdb_c2[1:30])

head(reportEnrichment(tempEnrich))
```

## 6. Customized genesets and recommendations

PanomiR can integrate genesets and pathways from external sources
including those annotated in MSigDB. In order to do so, you need to
provide a `GeneSetCollection` object as defined in the `GSEABase`
package.

The example below illustrates how to use external sources to create your
own customized pathway-gene association table. This customized can then
replaced `path_gene_table` input in functions described in sections 1,2,
and 5 of this manual.

``` r
data("gscExample")

newPathGeneTable <-tableFromGSC(gscExample)

head(newPathGeneTable)
```

The the pathway correlation network in section 3 is build upon an MSigDB
V6.2, canonical pathways (cp) collection dataset that includes KEGG
Pathways. KEGG prohibits distribution of its pathways by third parties.
However, use can access desired versions of MSigDB in gmt format via
[this link](https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)
(Subramanian et al. 2005).

The library `msigdb` provides an programmatic interface to download
different geneset collections. Including how to add KEGG pathways or
download mouse genesets. Use the this [MSigDB
tutorial](https://bioconductor.org/packages/release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html)
to create your desired gene sets.

You can also use the following code chunk to create pathway-gene
association tables from gmt files.

``` r
library(GSEABase)

yourGeneSetCollection <- getGmt("YOUR GMT FILE")
newPathGeneTable      <- tableFromGSC(yourGeneSetCollection)
```

## Session info

``` r
sessionInfo()
```

    ## R Under development (unstable) (2021-12-03 r81290)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.2.0  magrittr_2.0.1  fastmap_1.1.0   tools_4.2.0    
    ##  [5] htmltools_0.5.2 yaml_2.2.1      stringi_1.7.6   rmarkdown_2.11 
    ##  [9] knitr_1.37      stringr_1.4.0   xfun_0.29       digest_0.6.29  
    ## [13] rlang_0.4.12    evaluate_0.14

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-agarwal2015predicting" class="csl-entry">

Agarwal, Vikram, George W Bell, Jin-Wu Nam, and David P Bartel. 2015.
“Predicting Effective microRNA Target Sites in Mammalian mRNAs.” *Elife*
4: e05005.

</div>

<div id="ref-ally2017comprehensive" class="csl-entry">

Ally, Adrian, Miruna Balasundaram, Rebecca Carlsen, Eric Chuah, Amanda
Clarke, Noreen Dhalla, Robert A Holt, et al. 2017. “Comprehensive and
Integrative Genomic Characterization of Hepatocellular Carcinoma.”
*Cell* 169 (7): 1327–41.

</div>

<div id="ref-altschuler2013pathprinting" class="csl-entry">

Altschuler, Gabriel M, Oliver Hofmann, Irina Kalatskaya, Rebecca Payne,
Shannan J Ho Sui, Uma Saxena, Andrei V Krivtsov, et al. 2013.
“Pathprinting: An Integrative Approach to Understand the Functional
Basis of Disease.” *Genome Medicine* 5 (7): 1–13.

</div>

<div id="ref-joachim2018relative" class="csl-entry">

Joachim, Rose B, Gabriel M Altschuler, John N Hutchinson, Hector R Wong,
Winston A Hide, and Lester Kobzik. 2018. “The Relative Resistance of
Children to Sepsis Mortality: From Pathways to Drug Candidates.”
*Molecular Systems Biology* 14 (5): e7998.

</div>

<div id="ref-liberzon2011molecular" class="csl-entry">

Liberzon, Arthur, Aravind Subramanian, Reid Pinchback, Helga
Thorvaldsdóttir, Pablo Tamayo, and Jill P Mesirov. 2011. “Molecular
Signatures Database (MSigDB) 3.0.” *Bioinformatics* 27 (12): 1739–40.

</div>

<div id="ref-pita2018pathway" class="csl-entry">

Pita-Juárez, Yered, Gabriel Altschuler, Sokratis Kariotis, Wenbin Wei,
Katjuša Koler, Claire Green, Rudolph E Tanzi, and Winston Hide. 2018.
“The Pathway Coexpression Network: Revealing Pathway Relationships.”
*PLoS Computational Biology* 14 (3): e1006042.

</div>

<div id="ref-subramanian2005gene" class="csl-entry">

Subramanian, Aravind, Pablo Tamayo, Vamsi K Mootha, Sayan Mukherjee,
Benjamin L Ebert, Michael A Gillette, Amanda Paulovich, et al. 2005.
“Gene Set Enrichment Analysis: A Knowledge-Based Approach for
Interpreting Genome-Wide Expression Profiles.” *Proceedings of the
National Academy of Sciences* 102 (43): 15545–50.

</div>

</div>
