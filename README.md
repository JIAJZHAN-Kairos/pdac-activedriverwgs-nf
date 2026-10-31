# ActiveDriverWGS scatter/gather

Nextflow orchestration for official ActiveDriverWGS 1.2.1, GRCh38 and the
BSgenome.Hsapiens.UCSC.hg38 1.4.4 reference package. This repository contains
code only. Private mutation inputs, patient identifiers and results belong in
the authorized private S3 workspace and must not be committed here.

Independent elements are split into blocks. Each task receives the complete
mutation input for its chromosome, preserving all local background mutations.
Patients must be filtered for hypermutation using their **whole-genome**
SNV+indel counts before constructing these chromosome shards. Do not filter
patients separately within blocks. Each task invokes the official
`ADWGS_test` with the complete segmented element and a 50,000 bp background.
The container must be pinned by digest and provide the validated runtime at
`/opt/adwgs/Rlib` plus an AWS CLI for execution without Fusion.

The single JSON `--input` manifest has `schema_version: 1`, `max_vcpus`, a
`blocks` array (`block_id`, `mutations`, `elements`, `ids` and their SHA256s),
and `universe`, `membership`, `seed_results`, `expected_ids`, and
`independent_counts` paths. These files are
required for `--mode full`. Block ID files list unique, nonzero element IDs.
Seed results contain verified existing tests for identical inputs and elements.
Block results and seed IDs must cover the exact expected nonzero test set
once; the collector rejects duplicate, missing, or substituted test IDs and
compares observed burdens against independent unique-donor counts. All unmutated
elements are retained with P=1. BH is performed separately for each prespecified
complete analysis family, across all its chromosomes and region classes.
There is no per-block significance selection.

The pinned engine can count two indels from one donor when one indel crosses
an element/background boundary. The wrapper preserves native observed, expected
and P values. For this case it reconstructs the native interval-counting rules,
records `native_indel_boundary_excess`, and validates
`element_muts_obs = mutated_donors + native_indel_boundary_excess`.
Unexplained discrepancies still terminate execution. Report `mutated_donors`
as the independent donor support and flag candidates with native count excess;
the native observed count is not always a count of distinct donors.

Use `--mode pilot` first, comparing raw P, expected burden and observed donor
counts against independently generated native reference tests before starting
the complete cohort. Computational completion is not scientific validation.

```bash
nextflow run . -profile local \
  --input manifest.json --outdir output --mode full \
  --expected_fits 123 --max_parallel 64 --container registry/image@sha256:DIGEST
```

On Seqera, select the existing authorized AWS Batch compute environment,
provide private S3 paths, and set an independent S3 work directory. Concurrency
is limited to 64 single-vCPU fit tasks. The gather step starts after fitting.
Spot host reclamation is retried up to three times. Other failures terminate
and require inspection; no statistical parameter is changed automatically.
