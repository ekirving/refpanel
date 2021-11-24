## Configuration

If you wish to build a customised joint-callset (e.g., one that contains additional samples), `refpanel`
requires a minimal amount of configuration.

There are four tiers to this configuration:

* `panels`, which define a joint-callset containing `samples` from multiple `sources`
* `sources`, which contain a list of `samples` from a single datasource (
  e.g., [HGDP](https://www.internationalgenome.org/data-portal/data-collection/hgdp))
* `samples`, which aggregate one or more sequencing `accessions` (
  e.g., [HGDP00094](https://www.internationalgenome.org/data-portal/sample/HGDP00094))
* `accessions`, which uniquely identify an individual sequencing library (
  e.g., [ERR1373257](https://www.ebi.ac.uk/ena/browser/view/ERR1373257))

### Adding a new source

To add a new data `source` you will need to create two `tsv` metadata files:

* `<source>-samples.tsv`; required columns `(sample, sex)`
* `<source>-accessions.tsv`; required columns `(sample, accession, <readgroup-cols>, fastq_r1, fastq_r2)`

For examples, see [test-source-samples.tsv](../test/test-source-samples.tsv) and [test-source-accessions.tsv](
../test/test-source-accessions.tsv)

The mapping between `@RG` tags and `<readgroup-cols>` is defined in [config.yaml](../config.yaml)

By default, the following mapping is used:

* `ID:<accession>`
* `SM:<sample>`
* `CN:<center>`
* `PL:<platform>`
* `LB:<library>`
* `DS:<description>`

### Adding a new panel

To add a new reference `panel` you will need to create one additional `tsv` metadata file:

* `<panel>-samples.tsv` metadata sheet; required columns `(source, sample, superpopulation)`

For an example, see [test-panel-samples.tsv](../test/test-panel-samples.tsv)

