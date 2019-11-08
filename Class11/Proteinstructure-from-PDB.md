Untitled
================

Q1a. Determine the percentage of structures solved by X-Ray and Electron
Microscopy

``` r
data <- read.csv("Data Export Summary.csv")
sum <- sum(data$Total)
total <- data$Total
total[1]/sum * 100
```

    ## [1] 89.06177

``` r
total[3]/sum * 100
```

    ## [1] 2.514442

``` r
ans <- data$Total/sum(data$Total) * 100
names(ans) <- data$Experimental.Method
round(ans, 2)
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##               89.06                8.13                2.51 
    ##               Other        Multi Method 
    ##                0.19                0.10

Q1b. Also can you determine what proportion of structures are protein?

``` r
round(sum(data$Proteins)/sum(data$Total) * 100, 2)
```

    ## [1] 92.71

# Working with biomolecular data in R

``` r
library(bio3d)
pdb <- read.pdb("1hsg.pdb")
summary(pdb)
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
#Change all residue names to single letters
aa321(pdb$atom[pdb$calpha,"resid"])
```

    ##   [1] "P" "Q" "I" "T" "L" "W" "Q" "R" "P" "L" "V" "T" "I" "K" "I" "G" "G"
    ##  [18] "Q" "L" "K" "E" "A" "L" "L" "D" "T" "G" "A" "D" "D" "T" "V" "L" "E"
    ##  [35] "E" "M" "S" "L" "P" "G" "R" "W" "K" "P" "K" "M" "I" "G" "G" "I" "G"
    ##  [52] "G" "F" "I" "K" "V" "R" "Q" "Y" "D" "Q" "I" "L" "I" "E" "I" "C" "G"
    ##  [69] "H" "K" "A" "I" "G" "T" "V" "L" "V" "G" "P" "T" "P" "V" "N" "I" "I"
    ##  [86] "G" "R" "N" "L" "L" "T" "Q" "I" "G" "C" "T" "L" "N" "F" "P" "Q" "I"
    ## [103] "T" "L" "W" "Q" "R" "P" "L" "V" "T" "I" "K" "I" "G" "G" "Q" "L" "K"
    ## [120] "E" "A" "L" "L" "D" "T" "G" "A" "D" "D" "T" "V" "L" "E" "E" "M" "S"
    ## [137] "L" "P" "G" "R" "W" "K" "P" "K" "M" "I" "G" "G" "I" "G" "G" "F" "I"
    ## [154] "K" "V" "R" "Q" "Y" "D" "Q" "I" "L" "I" "E" "I" "C" "G" "H" "K" "A"
    ## [171] "I" "G" "T" "V" "L" "V" "G" "P" "T" "P" "V" "N" "I" "I" "G" "R" "N"
    ## [188] "L" "L" "T" "Q" "I" "G" "C" "T" "L" "N" "F"

``` r
#Call atoms to get rid of drug residues
ca.inds <- atom.select(pdb, "calpha")
ca.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
    ## 
    ##    Atom Indices#: 198  ($atom)
    ##    XYZ  Indices#: 594  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
lig <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(lig, file="1hsg_ligand.pdb")
```

``` r
prot <- atom.select(pdb, "protein", value = TRUE)
write.pdb(prot, file = "1hsg_prot.pdb")
```

``` r
library(bio3d.view)
view(lig)
```

    ## Computing connectivity from coordinates...

``` r
view(pdb, "overview", col = "sse")
```

    ## Computing connectivity from coordinates...
