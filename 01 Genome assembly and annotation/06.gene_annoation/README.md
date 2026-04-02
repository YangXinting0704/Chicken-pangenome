ab initio prediction use augustus
```
sh 01.ab initio prediction.sh

```
protein homology prediction use protein sequences from the GRC7b genome and Gallus protein sequences

```
sh 02.homology-based prediction.sh
```

transcriptome-based prediction use RNA-Seq data from 11 tissues and Iso-Seq data

```
sh 03.transcriptome-based prediction.sh

```

Intergrate ab initio prediction, transcriptome-based prediction, and homology-based prediction

```
sh 04.EVM.sh

```

Functional annotations use eggmapper

```
sh 05.eggmapper.sh

```