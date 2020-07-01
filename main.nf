/*
Copyright (c) 2020 Bisrat J. Debebe

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// TODO: survivor can be specified as `params.both.asInteger()`
nextflow.preview.dsl = 2


def baseOutdir = params.name ? "${params.outdir}/${params.name}" : "${params.outdir}"
println baseOutdir
println params.ref
println params.reads

process alignReads {
    publishDir  "${baseOutdir}/alignment", mode: 'copy'
    input:
        path reads
        path ref
    output:
        path 'alignment.sort.bam', emit: alignment

    "minimap2 -t ${task.cpus} -ax map-ont -MD ${ref} ${reads} | samtools sort -o alignment.sort.bam - "

}

process svimCall {
    publishDir "${baseOutdir}", mode: 'copy'
    input:
        path bam
        path ref
    output:
        path 'svim_calls'
    script:
    def sample_name = params.name ? "--sample ${params.name}" : ""

    "svim alignment --read_names ${sample_name} svim_calls ${bam} ${ref}"

}

process snifflesCall {
    publishDir  "${baseOutdir}/sniffles_calls", mode: 'copy'
    input:
        path bam
        path ref
    output:
        path '*sniffles.raw*'
    script:
    def sample_name = params.name ? "-v ${params.name}.sniffles.raw.vcf" : "sniffles.raw.vcf"
    def cluster = params.snifflesCluster ? "--cluster" : ""
    def genotype = params.snifflesGenotype ? "--genotype" : ""


    """ 
        sniffles -t ${task.cpus} -m ${bam} -s ${params.snifflesMinSupport} \
        --cluster_support ${params.snifflesClusterSupport}  \
        -l ${params.snifflesMinLength} -r ${params.snifflesMinSeqSize} \
        ${cluster} ${genotype} ${sample_name} ${params.snifflesAdvanced}
        bcftools sort -Oz -o ${sample_name}.gz ${sample_name} && tabix ${sample_name}.gz
    """

}

workflow {
    reads = Channel.fromPath(params.reads).ifEmpty { error "Please specify FASTQ files to use for variant calling" } 
    ref = Channel.fromPath(params.ref).ifEmpty { error "Please specify a reference file to use for alignment (can be remote file https/ftp)" } 
    alignReads(reads, ref)
    svimCall(alignReads.out.alignment, ref)
    snifflesCall(alignReads.out.alignment, ref)
}
