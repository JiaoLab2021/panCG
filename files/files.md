### panCNS

| Directory                          | File suffix                                | Describe                                                     |
| ---------------------------------- | ------------------------------------------ | ------------------------------------------------------------ |
| halLiftoverDir                     | {que}.{ref}.bed                            | que’s CNS position on the ref map <ref_map_chrID> <ref_map_start> <ref_map_end> <que_cnsID> |
| halLiftoverDir                     | {que}.{ref}.merge.bed                      | The position of que’s CNS in the ref map (merge if the distance is less than the threshold) <ref_map_chrID> <ref_map_start> <ref_map_end> <que_cnsID> |
| halLiftoverDir                     | {que}.{ref}.halLiftover.anchors            | The correspondence between que's CNS in the hal multiple sequence alignment file and ref's CNS <que_cnsID> <ref_cnsID> |
| halLiftoverDir                     | {que}.{ref}.merge.bw.bed                   | `{que}.{ref}.merge.bed` add averageBwScore and effecve_len   |
| blastnDir                          | .blastn.fmt6.txt                           | Original blastn alignment file                               |
| blastnDir                          | .blastn.halLiftoverFilter.anchors          | blastn anchors after halLiftover filtering                   |
| blastnDir                          | .blastn.halLiftoverFilter.anchors.fmt6.txt | fmt6 format of `.blastn.halLiftoverFilter.anchors` file      |
| JCVIDir                            | .lifted.anchors                            | Recruits additional anchors file output by JCVI              |
| JCVIDir                            | .anchors                                   | High quality anchors file output by JCVI                     |
| JCVIDir                            | .halLiftoverFilter.lifted.anchors          | JCVI recruits additional anchors after halLiftover filtering |
| JCVIDir                            | .halLiftoverFilter.anchors                 | JCVI high quality anchors after halLiftover filtering        |
| {Workdir}/Index/                   | cnsCluster.csv                             | Clustering information of CNS of all species                 |
| {Workdir}/Index/{species}/         | {species}.csv                              | Clustering information of all {species} CNS before           |
| {Workdir}/CEsDir/                  | {species}.recall_cds.bed                   | recall_cds coordinates for each species                      |
| {Workdir}/Ref\_{species}_IndexDir/ | .cnsIndexAssign.csv                        | Results of cnsIndexAssign                                    |
| {Workdir}/Ref\_{species}_IndexDir/ | .cnsIndexMerge.csv                         | Results of cnsIndexMerge                                     |
| {Workdir}/Ref\_{species}_IndexDir/ | .recallCEs.csv                             | The results of recallCEs, which records the results of CE obtained by recall |
| {Workdir}/Ref\_{species}_IndexDir/ | .ReCnsIndexMerge.csv                       | The result after merging `.recallCEs.csv`                    |
| {Workdir}/Ref\_{species}_IndexDir/ | .TripleCnsIndexMerge.csv                   | Classify the CE in `.ReCnsIndexMerge.csv` (cns, cds) and then merge the results |
| {Workdir}/Ref\_{species}_IndexDir/ | .recall.csv                                | The result of merging `.TripleCnsIndexMerge.csv` with cell   |
| {Workdir}/Ref\_{species}_IndexDir/ | .sort.csv                                  | The result of sort `.recall.csv`                             |

### 