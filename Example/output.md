

## 🔍 Output
### cns calling
| Directory               | File suffix        | Describe                            |
| ----------------------- | ------------------ | ----------------------------------- |
| {Workdir}/03-phastCons/ | {species}.all.bw   | PhastCons Conservative Scoring File |
| {Workdir}/03-phastCons/ | {species}.CNSs.bed | CNS file of {species}               |

### panCNS
| Directory             | File suffix        | Describe                                              |
| --------------------- | ------------------ | ----------------------------------------------------- |
| {Workdir}/Ref\_{ref}_ | .panGene.final.csv | The output panCNS file, each line represents an index |

### pangene
| Directory                     | File suffix  | Describe           |
| ----------------------------- | ------------ | ------------------ |
| {Workdir}/Ref\_{ref}_IndexDir | .panGene.csv | The result pangene |



The Group column is the homology group identified by orthofinder. 

| Group column        | Describe                                                     |
| ------------------- | ------------------------------------------------------------ |
| OGxxxxxxx.x         | Indicates the gene index subdivided in the homology group    |
| OGxxxxxxx.x.Un      | The .Un suffix indicates a set of genes that still exist independently in a single species after CPM. |
| OGxxxxxxx.x.tree_x  | Indicates the gene index subdivided by gene evolution relationship based on the gene index |
| OGxxxxxxx.x.tree_Un | The gene set ending with .tree_Un is a gene set that is not classified using evolutionary relationships. |
| UnMapOGXXXXXXX.x    | UnMap prefix is the gene that orthofinder has no clustering  |


