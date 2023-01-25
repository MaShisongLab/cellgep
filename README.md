# CellGEP

**An algorithm for transferring cell-type labels across single-cell datasets**

<b>CellGEP</b> (<u>Cel</u>l <u>l</u>abeling via <u>G</u>ene <u>E</u>xpression <u>P</u>rograms) is an algorithm for transferring cell-type label across single-cell datasets. Cell-type specific marker GEPs are first identified via single-cell gene co-expression network analysis via [scGGM](https://github.com/MaShisongLab/scggm) in a reference dataset(s). The identified marker GEPs can then be used to label and annotate a novel dataset for their cell types. For more details please refer to [Xu *et al*, 2023](#References). 

## Table of Contents
- [Install](#Install)
- [Usage](#Usage)
- [References](#References)

## Install
This algorithm works in R. It depends on the R packages [`Seurat`](https://satijalab.org/seurat/), [`graph`](https://www.bioconductor.org/packages/release/bioc/html/graph.html), [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html), and [`RBGL`](https://www.bioconductor.org/packages/release/bioc/html/RBGL.html). Copy the file `cellgep.R` to your R working folder, and source the file in R and start using it.


## Usage

The CellGEP algorithm uses marker GEPs identified in a reference dataset(s) to annotate cells in a novel single-cell dataset. We provide a list of marker GEPs as an example, which we identified from a mouse MCA dataset via scGGM ([Han *et al*, 2018](#References); [Xu *et al*, 2023](#References)). The GEPs are saved in the file 'data/mca.marker.geps.txt'. Below list the steps to use CellGEP to annotate a novel dataset.

### 1. Identify neighboring cells from the novel dataset.

Use the typical Seurat procedure to creat a Seurat single-cell analysis object for the novel dataset, e.g., named 'sc_object'; conduct PCA and UMAP analysis on the object; and then identify cell neighbors for the cells within the object via the following command in R.
```R
#R code
require(Seurat)
sc_object = FindNeighbors(sc_object, return.neighbor = T)
sc_neighbor = sc_object@neighbors$RNA.nn@nn.idx
sc_neighbor[1:5,]
```
sc_neighbor stores the neighbor information for the cells. Here is an example of the first five rows:
```shell
     [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20]
[1,]    1  1057   282   794   594   475   451   714   592   388  1302  1360   557   369   645   686   529  1485   948   353
[2,]    2 60882 15798 19477 14832 34639 35975 39103 24483 24851 35868 15407 15143 15612 34907 24539 61170 15680 23513 15155
[3,]    3  1302   388   217   461  1230  1327  1353   607  1339   948  1468   147  1014   339   895   784   451   936    68
[4,]    4  1372    64  1486   639   651   869   560   514  1411     6    23   408  1370  1460  1031  1485  1313   440   451
[5,]    5   204   708   657   127   961   372   262   800    65  1181  1048  1025   838  1497   352  1229   701   313  1287
```

### 2. Calculate marker GEPs' average expression values in the novel dataset.

The function `average_by_gep`calculates the marker GEP's average expression values in the novel dataset. </br>

<B>gep_exp = average_by_gep (`expression`, `geps`, `gep_origin`, `expression_origin`)</B></br>
`expression` - gene count matrix of the novel single-cell dataset, the row.names of matrix should be gene name. The format of gene names should macth those of the GEPs.</br>
`geps` - a table storing GEP infomration, with the 1st column being gene name, and the 2nd colum being GEP names. An example is store in the file "data/mca.marker.geps.txt".</br>
`gep_origin` - the origin of GEPs, e.g., "mca"</br>
`expression_origin` - the origin of expression, e.g., "novel_sc_dataset"
```R
# R code
require(edgeR)
source('cellgep.R')

# If gene IDs are in Ensembl ID format, the file 'mca.marker.geps.txt' should be used.
# If gene IDs are in gene symbl format, the file 'mca.marker.symbl.geps.txt' should be used instead.
geps = read.table('data/mca.marker.geps.txt')
geps_expression = average_by_gep( expression, geps, 'mca', 'novel_sc_dataset' )
```

### 3. Annotate the cells in the novel dataset via CellGEP. 

The `cellgep`function will using the GEPs to annotate the cells in the novel dataset.</br>

<B>anno = cellgep(`gep_expression`,`neighbor`,`geps_to_use`,`dynamic.cutoff=T`)</B></br>
`gep_expression` - the gep expression matrix</br>
`neighbor` - the cell neigbores, as calculated in step 1 </br>
`geps_to_use` - an integer array specifying which geps to be used for annotation</br>
`dyamic.cutoff` - whether to use a dynamic cutoff process to determine cutoff values of GEPs to be expressed

```R
# R code
require(graph)
require(RBGL)
source('cellgep.R')
gep_num = dim(geps_expression$expression)[1]
anno = cellgep(geps_expression$expression, sc_neighbor, 1:gep_num, dynamic.cutoff=T)
anno$annotation[1:5,1:13]
```
The table below lists the annotation results, with the first column being the maker GEPs annotation for the cells.
```shell
    GEP GEP_id is_original Expression GEP_mean GEP_sd GEPs_expressed GEP_Ori  N1  N2  N3  N4  N5  ...
1  M002      2           Y      1.014    0.107  0.121      M002-1.01       2   2   2   2   2   2  ...
2  M008      8           Y      0.957    0.078  0.141      M008-0.96       8  na   8  na   8  na  ...
3  M002      2           Y      0.793    0.107  0.121      M002-0.79       2   2   2   2  na   2  ...
4  M002      2           Y      0.929    0.107  0.121      M002-0.93       2   2   2   2   2   2  ...
5  M047     47           Y      0.875    0.058  0.065      M047-0.88      47  47  47  47  47  na  ...
... ...    ...          ...     ...      ...    ...           ...        ... ... ... ... ... ...  ...  
```
### 4. Visualize the results.

```R
# R code
sc_object@meta.data$marker_gep = anno$annotation[,1]
DimPlot(sc_object, group.by = 'marker_gep', label = T)
```

## References

Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, Saadatpour A, Zhou Z, Chen H, Ye F, Huang D, Xu Y, Huang W, Jiang M, Jiang X, Mao J, Chen Y, Lu C, Xie J, Fang Q, Wang Y, Yue R, Li T, Huang H, Orkin SH, Yuan GC, Chen M, and Guo G. 2018. Mapping the Mouse Cell Atlas by Microwell-Seq. *Cell* 172: 1091-1107.

Xu Y, Wang Y, and Ma S. 2023. scGGM identifies gene expression programs from single-cell transcriptomes and facilitates universal cell label transfer. *submitted*
