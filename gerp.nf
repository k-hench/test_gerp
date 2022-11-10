// ----- Config section -----
params.base = "${baseDir}"
params.outdir = "results"
params.indexBase = "reference"
params.splitFastaN = 1
params.lastParams = "-m 10 -j 3 -u 1 -p HOXD70"
params.minimap2Params = "-a -cx asm20"
params.gsalignParams = "-sen -no_vcf"
params.roastParams = "+ X=2 E="


log.info"""
GERP TEST
===================================
base_dir        : ${params.base}
outdir          : ${params.outdir}
refence name    : ${params.indexBase}
"""

nextflow.enable.dsl = 2

// ----- workflow components -----
Channel
  .fromPath("data/${params.indexBase}.fa")
  .set{ in_ch }

process build_last_db {
  publishDir "data/", mode: 'copy', pattern: "*_db.*"
  publishDir "${params.outdir}", mode: 'copy', pattern: "*.tsv"
  label "last_db"
  conda "${HOME}/miniconda3/envs/msa_align"

  input:
  file( in )

  output:
  tuple file ( "${params.indexBase}_db.*" ), file( "*_sizes.tsv" )

  script:
  """
  faSize -detailed ${in} > ${params.indexBase}_sizes.tsv
  lastdb -R 10 -u YASS -c ${params.indexBase}_db ${in}
  """
}

process gsalign_index {
  publishDir "data/", mode: 'copy', pattern: "*_idx*"
  label "gs_align_index"
  conda "${HOME}/miniconda3/envs/msa_align"

  input:
  file( in )

  output:
  file ( "${params.indexBase}_idx.*" )

  script:
  """
  faSize -detailed ${in} > ${params.indexBase}_sizes.tsv
  GSAlign index ${in} ${params.indexBase}_idx
  """
}

process build_index {
  publishDir "${params.outdir}/genome/nib", mode: 'copy'
  label "build_index"
  conda "${HOME}/miniconda3/envs/msa_ucsc"

  input:
  file( in )

  output:
  file ( "${params.indexBase}.2bit" )

  script:
  """
  faToTwoBit ${in} ${params.indexBase}.2bit
  """
}

process split_fasta {
  publishDir "data/", mode: 'copy'
  label "split_index"
  conda "${HOME}/miniconda3/envs/msa_split"

  input:
  file( in )

  output:
  tuple file( "*.split.fa" ), file( "*.fa.flat" ), file( "*.fa.gdx" )

  script:
  """
  pyfasta split -n ${params.splitFastaN} ${in}
  """
}

Channel
    .fromPath("data/query*")
    .set{ species_ch }

process align_single_last {
  publishDir "${params.outdir}/alignments", mode: 'copy'
  label "last_single"
  conda "${HOME}/miniconda3/envs/msa_align"

  input:
  tuple file( ref ), file( spec ), val( spec_name )

  output:
  file( "${spec_name}.psl" )

  script:
  """
  lastal \
    ${params.lastParams} \
    ${params.indexBase}_db \
    ${spec} | \
    maf-convert psl /dev/stdin > ${spec_name}.psl
  """
}

// ----- run workflow -----
workflow {
  main:
  last_ch = in_ch | build_last_db
  in_ch | gsalign_index
  in_ch | build_index
  in_ch | split_fasta
  last_ch.combine( species_ch ) | \
    map{ [ it[0], it[2], it[2].getName() - ".fa"] } | \
    align_single_last
}