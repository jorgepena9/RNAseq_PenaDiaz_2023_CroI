RNA-seq analysis - Pathway enrichment analysis
================
Jorge Peña Díaz
2023-05-21

------------------------------------------------------------------------

Analysis used to compare the gene expression profile of the ∆*croI* and
wild-type (WT) *C. rodentium* DB1000 strains as published in:

> Peña‑Díaz, J., Woodward, S.E., Creus‑Cuadros, A., Serapio‑Palacios,
> A., Deng, W., and Finlay, B.B. (2023). Quorum Sensing Modulates
> Bacterial Virulence and Colonization Dynamics During an Enteric
> Infection. bioRxiv,
> [10.1101/2023.03.14.532656](https://www.biorxiv.org/content/10.1101/2023.03.14.532656v1).

------------------------------------------------------------------------

### 1. Setup

Load the required packages for the analysis.

``` r
library(tidyverse)
library(EnrichmentBrowser)
library(clusterProfiler)
library(ggsci)
```

### 2. Import the required data

``` r
# Import differentially regulated genes
croI_res_sig <- read.delim("data/DEGs_dcroI.txt")

# Import gene expression
croI_res_all <- read.delim("data/ExpressionData_dcroI.txt")
```

### 3. Gene Set Enrichment Analysis (GSEA)

#### 3.1 Data preparation

Create the required files before analysis

``` r
# Create a ranked list of genes arranged according to the Wald-test statistic
croI_res_all %>%
  select(old_locus_tag, stat) %>%
  na.omit() %>% 
  filter(is.finite(stat)) %>%
  arrange(stat)  %>%
  write.table("data/croI.rnk",sep="\t",col.names = FALSE, row.names=FALSE,quote=FALSE)


# Create a gene set matrix (GMT) by using KEGG
cro <- getGenesets(org = "cro", db = "kegg", cache = TRUE)
```

    ## NCBI Gene ID mapping is not available for: cro

    ## Returning KEGG native gene IDs ...

``` r
# Write GMT file
writeGMT(cro, gmt.file = "data/Citrobacter_rodentium.gmt")
```

#### 3.2 GSEA

Run the gene set enrichment analysis in the terminal

    GSEA_4.2.3/gsea-cli.sh GSEAPreRanked \
    -gmx data/Citrobacter_rodentium.gmt \
    -rnk data/croI.rnk \
    -collapse false \
    -nperm 100000 \
    -scoring_scheme weighted \
    -rpt_label croI \
    -plot_top_x 20 \
    -rnd_seed 12345 \
    -set_max 500 \
    -set_min 5 \
    -zip_report false \
    -out croI_GSEA > gsea_output.txt

#### 3.3 Import results for data vizualization

``` r
# Import results
Upregulated <- read.delim("data/GSEA/gsea_report_for_na_pos.tsv") %>%
  mutate(Direction = "Upregulated")

Downregulated <- read.delim("data/GSEA/gsea_report_for_na_neg.tsv") %>%
   mutate(Direction = "Downregulated")

# Clean data
GSEA <- bind_rows(Upregulated, Downregulated) %>%
  mutate(Pathway = str_sub(NAME, end = 8)) %>%
  mutate(NAME = str_sub(NAME, start = 10)) %>%
  mutate(NAME = str_replace_all(NAME, "_", " ")) %>%
  mutate(NAME = str_to_sentence(NAME)) %>%
  select(NAME, Pathway, SIZE, ES, NES, NOM.p.val, q.val = "FDR.q.val", p.adj = "FWER.p.val",
         RANK.AT.MAX, LEADING.EDGE, Direction) %>%
  filter(q.val <= 0.05) 
```

### 4. Over-representation analysis (ORA)

#### 4.1 Data preparation

Create the required files before running the over-representation
analysis

``` r
#  List of all background genes
universe <- croI_res_all %>%
  as.data.frame() %>%
  select(old_locus_tag, stat) %>%
  na.omit() %>% 
  filter(is.finite(stat)) %>%
  arrange(desc(stat))

# Format the data
universe_geneList = universe[,2]
names(universe_geneList) = as.character(universe[,1])
universe_geneList <- names(sort(universe_geneList, decreasing = TRUE))

# List of differentially regulated genes
d <- croI_res_sig %>%
  select(old_locus_tag, stat) %>%
  na.omit() %>% 
  filter(is.finite(stat)) %>%
  arrange(desc(stat)) %>%
  as.data.frame()

# Format the data
geneList = d[,2]
names(geneList) = as.character(d[,1])

# Subset up- and down-regulated genes
geneList_Up <- geneList[geneList > 0]
geneList_Down <- geneList[geneList < 0]

# Sort in decreasing order
geneList_Up <- names(sort(abs(geneList_Up), decreasing = TRUE))
geneList_Down <- names(sort(abs(geneList_Down), decreasing = TRUE))
```

#### 3.2 ORA (KEGG pathways)

KEGG pathway over-representation analysis

``` r
KEGG_up <- enrichKEGG(gene = geneList_Up,
                      organism = 'cro',
                      minGSSize = 5,
                      universe = universe_geneList,
                      pvalueCutoff = 0.05)
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/cro/pathway"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/cro"...

``` r
KEGG_down <- enrichKEGG(gene = geneList_Down,
                        organism = 'cro',
                        minGSSize = 5,
                        universe = universe_geneList,
                        pvalueCutoff = 0.05)

# Merge data
KEGG_up <- KEGG_up %>%
  mutate(Direction = "Upregulated") %>%
  as.data.frame()

KEGG_down <- KEGG_down %>%
  mutate(Direction = "Downregulated") %>%
  as.data.frame()

KEGG <- bind_rows(KEGG_up, KEGG_down)
```

#### 3.2 ORA (KEGG modules)

KEGG module over-representation analysis

``` r
# KEGG Module over-representation
Module_up <- enrichMKEGG(gene = geneList_Up,
                         organism = 'cro',
                         universe = universe_geneList,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05)
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/cro/module"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/module"...

``` r
Module_down <- enrichMKEGG(gene = geneList_Down,
                           organism = 'cro',
                           universe = universe_geneList,
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

# Merge data
Module_up <- Module_up %>%
  mutate(Direction = "Upregulated") %>%
  as.data.frame()

Module_down <- Module_down %>%
  mutate(Direction = "Downregulated") %>%
  as.data.frame()

Module <- bind_rows(Module_up, Module_down)
```

### 5. Data visualization

We will combine the results from GSEA and module over-representation
analysis and then plot them in a lollipop plot.

``` r
# Clean GSEA data
GSEA_clean <- GSEA %>%
  mutate(p.adj = ifelse(p.adj == 0, 0.00001, p.adj)) %>%  # Replace since we cannot calculate actual p value
  mutate(q.val = ifelse(q.val == 0, 0.0000000001, q.val)) %>% # Replace since we cannot calculate actual q value
  select(NAME, Pathway, p.adj, q.val, Direction, SIZE) %>%
  mutate(Analysis = "GSEA")

# Clean KEGG Module over-representation data
Module_clean <- Module %>%
  rename(NAME = "Description", Pathway = "ID", p.adj = "p.adjust", q.val = "qvalue") %>%
  select(NAME, Pathway, p.adj, q.val, Direction, Count) %>%
  rename(SIZE = "Count") %>%
  mutate(Analysis = "Module")

# Merge the data and do final cleaning before plotting
Data_clean <- rbind(GSEA_clean, Module_clean) %>%
  mutate(NAME = str_replace_all(NAME, " ", "\n")) %>%
  mutate(NAME = str_replace(NAME, "EHEC/EPEC\npathogenicity\nsignature,\nT3SS\nand\neffectors", 
                                   "EHEC/EPEC\npathogenicity signature\n(T3SS and effectors)")) %>%
  mutate(NAME = str_replace(NAME, "Starch\nand\nsucrose\nmetabolism",
                            "Starch and\nsucrose\nmetabolism")) %>%
  mutate(qscore = -log10(p.adj)) %>% 
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) 
  
# Plot the data
Data_clean %>%
  ggplot(aes(x = qscore, reorder(NAME, q.val))) +
  geom_linerange(aes(xmin = 0, xmax = qscore, colour = Analysis, y = reorder(NAME, q.val, decreasing = TRUE)), linewidth = 2.5) +
  geom_point(aes(size = SIZE, colour = Analysis)) +
  xlab("-log10(p-value)") +
  ylab("") +
  facet_grid(Direction~., scales = "free", space = "free") +
  labs(colour = "Analysis") +
  theme_bw() +
  scale_y_discrete() +
  scale_colour_d3(labels = c("GSEA (KEGG Pathways)", "ORA (KEGG Modules)")) +
  scale_x_continuous(limits = c(0,25)) +
  scale_size(range = c(1.5, 9.5), name="Size", limits = c(5, 60)) +
  theme(plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.8)),
        axis.text.y = element_text(size = rel(1.8), hjust = 1),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.15)),
        strip.text = element_text(size = rel(1.125)))  +
  guides(colour = guide_legend(order = 1))
```

![](Pathway_enrichment_files/figure-gfm/GSEA_Module%20data-1.png)<!-- -->

------------------------------------------------------------------------

### 6. Session information

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggsci_3.0.0                 clusterProfiler_4.6.2      
    ##  [3] EnrichmentBrowser_2.28.2    graph_1.76.0               
    ##  [5] SummarizedExperiment_1.28.0 Biobase_2.58.0             
    ##  [7] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
    ##  [9] IRanges_2.32.0              S4Vectors_0.36.2           
    ## [11] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
    ## [13] matrixStats_0.63.0          lubridate_1.9.2            
    ## [15] forcats_1.0.0               stringr_1.5.0              
    ## [17] dplyr_1.1.2                 purrr_1.0.1                
    ## [19] readr_2.1.4                 tidyr_1.3.0                
    ## [21] tibble_3.2.1                ggplot2_3.4.2              
    ## [23] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] shadowtext_0.1.2       fastmatch_1.1-3        BiocFileCache_2.6.1   
    ##   [4] plyr_1.8.8             igraph_1.4.2           lazyeval_0.2.2        
    ##   [7] GSEABase_1.60.0        splines_4.2.3          BiocParallel_1.32.6   
    ##  [10] digest_0.6.31          yulab.utils_0.0.6      htmltools_0.5.5       
    ##  [13] GOSemSim_2.24.0        viridis_0.6.3          GO.db_3.16.0          
    ##  [16] fansi_1.0.4            magrittr_2.0.3         memoise_2.0.1         
    ##  [19] tzdb_0.4.0             Biostrings_2.66.0      annotate_1.76.0       
    ##  [22] graphlayouts_1.0.0     timechange_0.2.0       enrichplot_1.18.4     
    ##  [25] colorspace_2.1-0       blob_1.2.4             rappdirs_0.3.3        
    ##  [28] ggrepel_0.9.3          xfun_0.39              crayon_1.5.2          
    ##  [31] RCurl_1.98-1.12        jsonlite_1.8.4         scatterpie_0.1.9      
    ##  [34] ape_5.7-1              glue_1.6.2             polyclip_1.10-4       
    ##  [37] gtable_0.3.3           zlibbioc_1.44.0        XVector_0.38.0        
    ##  [40] DelayedArray_0.24.0    Rgraphviz_2.42.0       scales_1.2.1          
    ##  [43] DOSE_3.24.2            DBI_1.1.3              Rcpp_1.0.10           
    ##  [46] viridisLite_0.4.2      xtable_1.8-4           gridGraphics_0.5-1    
    ##  [49] tidytree_0.4.2         bit_4.0.5              httr_1.4.6            
    ##  [52] fgsea_1.24.0           RColorBrewer_1.1-3     pkgconfig_2.0.3       
    ##  [55] XML_3.99-0.14          farver_2.1.1           dbplyr_2.3.2          
    ##  [58] utf8_1.2.3             ggplotify_0.1.0        tidyselect_1.2.0      
    ##  [61] labeling_0.4.2         rlang_1.1.1            reshape2_1.4.4        
    ##  [64] AnnotationDbi_1.60.2   munsell_0.5.0          tools_4.2.3           
    ##  [67] cachem_1.0.8           downloader_0.4         cli_3.6.1             
    ##  [70] generics_0.1.3         RSQLite_2.3.1          gson_0.1.0            
    ##  [73] evaluate_0.21          fastmap_1.1.1          yaml_2.3.7            
    ##  [76] ggtree_3.6.2           knitr_1.42             bit64_4.0.5           
    ##  [79] tidygraph_1.2.3        KEGGREST_1.38.0        ggraph_2.1.0          
    ##  [82] nlme_3.1-162           KEGGgraph_1.58.3       aplot_0.1.10          
    ##  [85] compiler_4.2.3         rstudioapi_0.14        filelock_1.0.2        
    ##  [88] curl_5.0.0             png_0.1-8              treeio_1.22.0         
    ##  [91] tweenr_2.0.2           stringi_1.7.12         highr_0.10            
    ##  [94] lattice_0.21-8         Matrix_1.5-4           vctrs_0.6.2           
    ##  [97] pillar_1.9.0           lifecycle_1.0.3        data.table_1.14.8     
    ## [100] cowplot_1.1.1          bitops_1.0-7           patchwork_1.1.2       
    ## [103] qvalue_2.30.0          R6_2.5.1               gridExtra_2.3         
    ## [106] codetools_0.2-19       MASS_7.3-60            withr_2.5.0           
    ## [109] GenomeInfoDbData_1.2.9 parallel_4.2.3         hms_1.1.3             
    ## [112] grid_4.2.3             ggfun_0.0.9            HDO.db_0.99.1         
    ## [115] rmarkdown_2.21         ggforce_0.4.1

------------------------------------------------------------------------

### 7. References

``` r
packages <- c("tidyverse", "EnrichmentBrowser", "clusterProfiler", "ggsci")

packages %>%
  map(citation) %>%
  print(style = "text")
```

    ## [[1]]
    ## Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
    ## Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E,
    ## Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi
    ## K, Vaughan D, Wilke C, Woo K, Yutani H (2019). "Welcome to the
    ## tidyverse." _Journal of Open Source Software_, *4*(43), 1686.
    ## doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
    ## 
    ## [[2]]
    ## Geistlinger L, Csaba G, Zimmer R (2016). "Bioconductor's
    ## EnrichmentBrowser: seamless navigation through combined results of set-
    ## & network-based enrichment analysis." _BMC Bioinformatics_, *17*, 45.
    ## doi:10.1186/s12859-016-0884-1
    ## <https://doi.org/10.1186/s12859-016-0884-1>.
    ## 
    ## [[3]]
    ## Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L,
    ## Fu x, Liu S, Bo X, Yu G (2021). "clusterProfiler 4.0: A universal
    ## enrichment tool for interpreting omics data." _The Innovation_, *2*(3),
    ## 100141. doi:10.1016/j.xinn.2021.100141
    ## <https://doi.org/10.1016/j.xinn.2021.100141>.
    ## 
    ## Yu G, Wang L, Han Y, He Q (2012). "clusterProfiler: an R package for
    ## comparing biological themes among gene clusters." _OMICS: A Journal of
    ## Integrative Biology_, *16*(5), 284-287. doi:10.1089/omi.2011.0118
    ## <https://doi.org/10.1089/omi.2011.0118>.
    ## 
    ## [[4]]
    ## Xiao N (2023). _ggsci: Scientific Journal and Sci-Fi Themed Color
    ## Palettes for 'ggplot2'_. R package version 3.0.0,
    ## <https://CRAN.R-project.org/package=ggsci>.
