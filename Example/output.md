

## 🔍 Output
### callCns
| Directory               | File suffix           | Describe                                                     |
| ----------------------- | --------------------- | ------------------------------------------------------------ |
| {Workdir}/02-model/     | nonconserved_4d.mod   | The non-conservative model obtained by phyloFit              |
| {Workdir}/02-model/     | allCEs.sort.merge.bed | Coordinates of all conserved sequences in the genome, including CDS. |
| {Workdir}/03-phastCons/ | {species}.all.bw      | PhastCons Conservative Scoring File                          |
| {Workdir}/03-phastCons/ | {species}.CNSs.bed    | CNS file of {species}                                        |
### pangene
| Directory                     | File suffix                    | Describe             |
| ----------------------------- | ------------------------------ | -------------------- |
| {Workdir}/diamondDir          | .diamond_blastp.filter.anchors | blastp anchors       |
| {Workdir}/JCVIDir             | .graph.pkl                     | Gene synteny network |
| {Workdir}/Ref\_{ref}_IndexDir | .panGene.csv                   | The result pangene   |

### panCNS
| Directory                | File suffix                              | Describe                                                     |
| ------------------------ | ---------------------------------------- | ------------------------------------------------------------ |
| {Workdir}/halLiftoverDir | .halLiftover.anchors                     | Homologous CNS combining evidence from multiple sequence alignment. |
| {Workdir}/blastnDir      | .blastn.halLiftoverFilter.anchors        | Homologous CNS combining evidence from multiple sequence alignment, similarity. |
| {Workdir}/JCVIDir        | .halLiftoverFilter.rescue.lifted.anchors | Homologous CNS combining evidence from multiple sequence alignment, similarity, and collinearity |
| {Workdir}/Ref\_{ref}_    | .pancns.sort.csv                         | The output panCNS file, each line represents an index        |

The Group column is the homology group identified by orthofinder. 

| Group column        | Describe                                                     |
| ------------------- | ------------------------------------------------------------ |
| OGxxxxxxx.x         | Indicates the gene index subdivided in the homology group    |
| OGxxxxxxx.x.Un      | The .Un suffix indicates a set of genes that still exist independently in a single species after CPM. |
| OGxxxxxxx.x.tree_x  | Indicates the gene index subdivided by gene evolution relationship based on the gene index |
| OGxxxxxxx.x.tree_Un | The gene set ending with .tree_Un is a gene set that is not classified using evolutionary relationships. |
| UnMapOGXXXXXXX.x    | UnMap prefix is the gene that orthofinder has no clustering  |

