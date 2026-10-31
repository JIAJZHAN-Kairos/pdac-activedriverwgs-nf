nextflow.enable.dsl=2
import groovy.json.JsonSlurper

process FIT_BLOCK {
    tag "${block_id}"
    cpus 1
    memory '4 GB'
    time '2h'
    maxForks params.max_parallel
    publishDir "${params.outdir}/blocks", mode: 'copy'
    input:
    tuple val(block_id), path(mutations), path(elements), path(ids), val(mutations_sha), val(elements_sha), val(ids_sha)
    output:
    path "${block_id}.tsv", emit: tests
    path "${block_id}.qc.json", emit: qc
    script:
    """
    echo '${mutations_sha}  ${mutations}' | sha256sum -c -
    echo '${elements_sha}  ${elements}' | sha256sum -c -
    echo '${ids_sha}  ${ids}' | sha256sum -c -
    export R_LIBS=/opt/adwgs/Rlib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    Rscript ${projectDir}/bin/fit_block.R '${mutations}' '${elements}' '${ids}' '${block_id}'
    """
    stub:
    """
    printf 'id\\tpp_element\\n' > ${block_id}.tsv
    echo '{}' > ${block_id}.qc.json
    """
}

process COLLECT_TESTS {
    cpus 1
    memory '8 GB'
    time '1h'
    publishDir params.outdir, mode: 'copy'
    input:
    path tests, stageAs: 'chunk_tests/*'
    path universe
    path membership
    path seed
    path expected_ids
    output:
    path 'results', emit: results
    script:
    """
    export R_LIBS=/opt/adwgs/Rlib OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    Rscript ${projectDir}/bin/collect_tests.R chunk_tests '${universe}' '${membership}' '${seed}' '${params.expected_fits}' '${expected_ids}'
    """
    stub:
    """
    mkdir results
    echo 'stub only' > results/STUB_ONLY
    """
}

workflow {
    if(!params.input || !params.outdir || !params.container) error 'Required: --input, --outdir, --container'
    def manifest = new JsonSlurper().parseText(file(params.input).text)
    if(manifest.schema_version != 1) error 'Unsupported manifest schema'
    if(manifest.max_vcpus > 64 || params.max_parallel > 64) error 'Authorized ceiling is 64 vCPU'
    def blocks = Channel.fromList(manifest.blocks).map { b ->
        tuple(b.block_id,file(b.mutations),file(b.elements),file(b.ids),b.mutations_sha256,b.elements_sha256,b.ids_sha256)
    }
    FIT_BLOCK(blocks)
    if(params.mode == 'full') {
        COLLECT_TESTS(FIT_BLOCK.out.tests.collect(),file(manifest.universe),file(manifest.membership),file(manifest.seed_results),file(manifest.expected_ids))
    } else {
        FIT_BLOCK.out.tests.subscribe { log.info "Pilot result: ${it}" }
    }
}
