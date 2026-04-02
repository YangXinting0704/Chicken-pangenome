eliminate redundant sequences
```
sh purge_haplotigs

```
validate and anchor contigs onto chromosomes

```
sh AssemblyMapper.sh
```
Illumina read pairs are mapped to W chromosome contigs,calculated the read depth and considered which contigs originate from the W chromosome

```
sh mosdepth.sh
```

place the scaffolds onto chromosomes based on their interactions

```
juicer_3D-DNA.sh
```



