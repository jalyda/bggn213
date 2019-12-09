Cancer Genomics
================

``` r
library(GenomicDataCommons)
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'GenomicDataCommons'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(TCGAbiolinks)
library(maftools)
```

``` r
status()
```

    ## $commit
    ## [1] "955a105f3f2ba797e1d9d8de013226a495feae56"
    ## 
    ## $data_release
    ## [1] "Data Release 20.0 - November 11, 2019"
    ## 
    ## $status
    ## [1] "OK"
    ## 
    ## $tag
    ## [1] "1.23.0"
    ## 
    ## $version
    ## [1] 1

``` r
projects <- getGDCprojects()
head(projects)
```

    ##   dbgap_accession_number
    ## 1                   <NA>
    ## 2                   <NA>
    ## 3              phs000466
    ## 4                   <NA>
    ## 5              phs000467
    ## 6              phs000465
    ##                                                                                                                                           disease_type
    ## 1                                                                                                                                Mesothelial Neoplasms
    ## 2                                                                                                                         Adenomas and Adenocarcinomas
    ## 3                                                                                                                  Complex Mixed and Stromal Neoplasms
    ## 4 Myomatous Neoplasms, Soft Tissue Tumors and Sarcomas, NOS, Fibromatous Neoplasms, Lipomatous Neoplasms, Nerve Sheath Tumors, Synovial-like Neoplasms
    ## 5                                                                                                       Neuroepitheliomatous Neoplasms, Not Applicable
    ## 6                                                                                                                    Myeloid Leukemias, Not Applicable
    ##   releasable released state
    ## 1      FALSE     TRUE  open
    ## 2      FALSE     TRUE  open
    ## 3      FALSE     TRUE  open
    ## 4      FALSE     TRUE  open
    ## 5       TRUE     TRUE  open
    ## 6       TRUE     TRUE  open
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   primary_site
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Heart, mediastinum, and pleura, Bronchus and lung
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Adrenal gland
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Kidney
    ## 4                                                                                                                                                                                                                                                                                        Corpus uteri, Stomach, Other and unspecified parts of tongue, Meninges, Other and unspecified male genital organs, Colon, Connective, subcutaneous and other soft tissues, Bones, joints and articular cartilage of limbs, Ovary, Retroperitoneum and peritoneum, Peripheral nerves and autonomic nervous system, Uterus, NOS, Kidney
    ## 5 Heart, mediastinum, and pleura, Stomach, Bones, joints and articular cartilage of other and unspecified sites, Lymph nodes, Liver and intrahepatic bile ducts, Unknown, Uterus, NOS, Skin, Other endocrine glands and related structures, Adrenal gland, Renal pelvis, Connective, subcutaneous and other soft tissues, Bones, joints and articular cartilage of limbs, Other and ill-defined sites, Meninges, Spinal cord, cranial nerves, and other parts of central nervous system, Retroperitoneum and peritoneum, Peripheral nerves and autonomic nervous system, Hematopoietic and reticuloendothelial systems, Kidney
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Unknown, Hematopoietic and reticuloendothelial systems
    ##    project_id          id                             name tumor
    ## 1   TCGA-MESO   TCGA-MESO                     Mesothelioma  MESO
    ## 2    TCGA-ACC    TCGA-ACC         Adrenocortical Carcinoma   ACC
    ## 3 TARGET-CCSK TARGET-CCSK Clear Cell Sarcoma of the Kidney  CCSK
    ## 4   TCGA-SARC   TCGA-SARC                          Sarcoma  SARC
    ## 5  TARGET-NBL  TARGET-NBL                    Neuroblastoma   NBL
    ## 6  TARGET-AML  TARGET-AML           Acute Myeloid Leukemia   AML

``` r
cases_by_project <- cases() %>%
  facet("project.project_id") %>%
  aggregations()
head(cases_by_project)
```

    ## $project.project_id
    ##                      key doc_count
    ## 1                  FM-AD     18004
    ## 2             TARGET-AML      2146
    ## 3          TARGET-ALL-P2      1587
    ## 4             TARGET-NBL      1132
    ## 5              TCGA-BRCA      1098
    ## 6          MMRF-COMMPASS       995
    ## 7              TARGET-WT       652
    ## 8               TCGA-GBM       617
    ## 9                TCGA-OV       608
    ## 10             TCGA-LUAD       585
    ## 11     BEATAML1.0-COHORT       583
    ## 12             TCGA-UCEC       560
    ## 13             TCGA-KIRC       537
    ## 14             TCGA-HNSC       528
    ## 15              TCGA-LGG       516
    ## 16             TCGA-THCA       507
    ## 17             TCGA-LUSC       504
    ## 18             TCGA-PRAD       500
    ## 19          NCICCR-DLBCL       489
    ## 20             TCGA-SKCM       470
    ## 21             TCGA-COAD       461
    ## 22             TCGA-STAD       443
    ## 23             TCGA-BLCA       412
    ## 24             TARGET-OS       383
    ## 25             TCGA-LIHC       377
    ## 26               CPTAC-2       342
    ## 27               CPTAC-3       322
    ## 28             TCGA-CESC       307
    ## 29             TCGA-KIRP       291
    ## 30             TCGA-SARC       261
    ## 31             TCGA-LAML       200
    ## 32         TARGET-ALL-P3       191
    ## 33             TCGA-ESCA       185
    ## 34             TCGA-PAAD       185
    ## 35             TCGA-PCPG       179
    ## 36              OHSU-CNL       176
    ## 37             TCGA-READ       172
    ## 38             TCGA-TGCT       150
    ## 39             TCGA-THYM       124
    ## 40            CGCI-BLGSP       120
    ## 41             TCGA-KICH       113
    ## 42              TCGA-ACC        92
    ## 43             TCGA-MESO        87
    ## 44              TCGA-UVM        80
    ## 45   ORGANOID-PANCREATIC        70
    ## 46             TARGET-RT        69
    ## 47             TCGA-DLBC        58
    ## 48              TCGA-UCS        57
    ## 49 BEATAML1.0-CRENOLANIB        56
    ## 50             TCGA-CHOL        51
    ## 51           CTSP-DLBCL1        45
    ## 52         TARGET-ALL-P1        24
    ## 53           TARGET-CCSK        13
    ## 54             HCMI-CMDC         7
    ## 55        VAREPOP-APOLLO         7

``` r
library(bio3d)
library()
```

``` r
seqs <- read.fasta("lecture18_sequences.fa")
seqs
```

    ##              1        .         .         .         .         .         60 
    ## P53_wt       MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    ## P53_mutant   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
    ##              **************************************** ******************* 
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## P53_wt       DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ## P53_mutant   DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ##              **** ******************************************************* 
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## P53_wt       SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ## P53_mutant   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ##              ************************************************************ 
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .         .         .         240 
    ## P53_wt       RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ## P53_mutant   RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ##              ******************************** *************************** 
    ##            181        .         .         .         .         .         240 
    ## 
    ##            241        .         .         .         .         .         300 
    ## P53_wt       SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    ## P53_mutant   SCMGGMNRRPILTIITLEV-----------------------------------------
    ##              ******************                                           
    ##            241        .         .         .         .         .         300 
    ## 
    ##            301        .         .         .         .         .         360 
    ## P53_wt       PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    ## P53_mutant   ------------------------------------------------------------
    ##                                                                           
    ##            301        .         .         .         .         .         360 
    ## 
    ##            361        .         .         .  393 
    ## P53_wt       GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    ## P53_mutant   ---------------------------------
    ##                                                
    ##            361        .         .         .  393 
    ## 
    ## Call:
    ##   read.fasta(file = "lecture18_sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 393 position columns (259 non-gap, 134 gap) 
    ## 
    ## + attr: id, ali, call

``` r
id <- conserv(seqs, method = "identity")
which(id != 1)
```

    ##   [1]  41  65 213 259 260 261 262 263 264 265 266 267 268 269 270 271 272
    ##  [18] 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289
    ##  [35] 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306
    ##  [52] 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323
    ##  [69] 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340
    ##  [86] 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357
    ## [103] 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374
    ## [120] 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391
    ## [137] 392 393

``` r
pos <- which(id != 1)[1]
```
