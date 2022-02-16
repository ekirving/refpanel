## Pipeline steps

The pipeline can also be broken down into separate steps, for distribution across multiple nodes in a cluster.

### Downloading data

To download all the `CRAM` files for the four IGSR data sources, run:

```bash
for source in 1000g hgdp sgdp ggvp; do
  ./refpanel data/source/${source}/cram/download.done
done
```

To download all the `gVCF` files for 1000G, run:

```bash
./refpanel data/source/1000g/gVCF/download.done
```

To download all the `fastq` files for the other data sources, run:

```bash
for source in appg PRJEB26721 PRJEB28198; do
  ./refpanel data/source/${source}/fastq/download.done
done
```

To avoid issues with rate limits on the EBI FTP service, downloads are limited to 15 concurrent FTP connections. If your
cluster nodes have unique IP addresses, then you may wish to run each of these download commands on separate nodes.

### Aligning CRAMs

To align all samples in a data source, run:

```bash
./refpanel data/source/${source}/cram/align.done
```

### QC metrics

To calculate alignment QC metrics for all samples in a data source, run:

```bash
./refpanel data/source/${source}/cram/metrics.done
```

### Calling gVCFs

To call genotype likelihoods in all samples in a data source, run:

```bash
./refpanel data/source/${source}/gVCF/call.done
```

### Merging gVCFs

To merge all `gVCF` files in a data source, run:

```bash
./refpanel data/source/${source}/gVCF/merge.done
```

This is one of the most time-consuming steps, as `gatk CombineGVCFs` is extremely slow when there are thousands of
samples.

### Joint-calling VCFs

To joint-call a reference panel, run:

```bash
./refpanel data/panel/${panel}/vcf/joint-call.done
```

### Read-based phasing

To perform read-based phasing of all samples in a reference panel, run:

```bash
./refpanel data/panel/${panel}/vcf/whatshap.done
```

### Pedigree phasing

To perform pedigree based phasing of all trios in a reference panel, run:

```bash
./refpanel data/panel/${panel}/vcf/trios.done
```

### Statistical phasing

To perform statistical phasing of all samples a reference panel, run:

```bash
./refpanel data/panel/${panel}/vcf/phase.done
```

