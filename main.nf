/*
Copyright (c) 2020 Bisrat J. Debebe

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// TODO: survivor can be specified as `params.both.asInteger()`
nextflow.preview.dsl = 2


def baseOutdir = params.name ? "${params.outdir}/${params.name}" : "${params.outdir}"

process alignReads {
    publishDir  "${baseOutdir}/alignment", mode: 'copy'
    input:
        path reads
        path ref
    output:
        path '*alignment.sort.bam', emit: alignment
        path 'lowmq*'
    script:
    def bamName = params.name ? "${params.name}.alignment.sort.bam" :"alignment.sort.bam" 
    """
        minimap2 --cs --MD -t ${task.cpus} -ax map-ont ${ref} ${reads} | samtools sort -o ${bamName} - 
        samtools view -h ${bamName} | awk '\$1 ~ /^@/ || \$5<5 '  | samtools sort -o lowmq.sort.bam -
        samtools depth lowmq.sort.bam > lowq.depth
        SURVIVOR bincov ${params.windowSize} ${params.lowQualSupport} 2 > lowmq_regions.bed
    """

}

process svimCalls {
    publishDir "${baseOutdir}/sv/svim_calls", mode: 'copy'
    input:
        path bam
        path ref
        file regions
    output:
        path 'svim.results.tar.gz', emit: svim_dir
        path "**.vcf.gz*", emit: svim_vcf
        path "*.log", emit: svim_log
        path "*.png", emit: svim_images
    script:
        def (sampleCmd, outRawVcf, outBndVcf, outFilteredVcf) = params.name ? 
            [ "--sample ${params.name}", "${params.name}.svim.raw.vcf.gz",
            "${params.name}.svim.bnd.vcf.gz", "${params.name}.svim.filtered.vcf.gz",
            ] 
            : 
            [ "", "svim.raw.vcf.gz", "svim.bnd.vcf.gz", "svim.filtered.vcf.gz"]


        """
            svim alignment --read_names ${sampleCmd} svim_calls ${bam} ${ref}
            tar czf svim.results.tar.gz svim_calls 
            mv svim_calls/*{log,png} .
            bcftools sort -T ${params.tmpDir} -Oz -o ${outRawVcf} svim_calls/variants.vcf 
            bcftools filter -i 'FILTER=="PASS" && SVTYPE=="BND"' ${outRawVcf} | bcftools sort -Oz  -o ${outBndVcf} - 
            bcftools filter -i "${params.svimFilter}" ${outRawVcf} | bcftools sort -Oz -o ${outFilteredVcf}  -
        """

}

        
process snifflesCalls {
    publishDir  "${baseOutdir}/sv/sniffles_calls", mode: 'copy'
    input:
        path bam
        path ref
    output:
        path '*sniffles.raw*', emit: sniffles_vcf
    script:
    def vcfName = params.name ? "${params.name}.sniffles.raw.vcf" : "sniffles.raw.vcf"
    def cluster = params.snifflesCluster ? "--cluster" : ""
    def genotype = params.snifflesGenotype ? "--genotype" : ""


    """ 
        sniffles -t ${task.cpus} -m ${bam} -s ${params.snifflesMinSupport} \
        --cluster_support ${params.snifflesClusterSupport} -v ${vcfName}  \
        -l ${params.snifflesMinLength} -r ${params.snifflesMinSeqSize} \
        ${cluster} ${genotype} ${params.snifflesAdvanced}
        bcftools sort -T ${params.tmpDir} -Oz -o ${vcfName}.gz ${vcfName} && tabix ${vcfName}.gz
    """

}

process highConfCalls {
    publishDir  "${baseOutdir}/sv/high_conf", mode: 'copy'
    input:
        path svim_calls
        path sniffles_calls
    output:
        path '*highconf.vcf.gz', emit: highconf_vcf
    script:
    def vcfName = params.name ? "${params.name}.highconf.vcf" :"highconf.vcf" 
    def sameStrand = params.sameStrand ? 1 : 0
    def sameType = params.sameType ? 1 : 0
    def estDist = params.estDist ? 1 : 0

    """ 
    SURVIVOR merge <(ls ${svim_calls} ${sniffles_calls}) ${params.isecDist}\
    ${params.callerSupport} ${sameStrand} ${sameType} ${estDist}\
    ${params.isecMinLength} ${vcfName}
    bcftools sort -Oz  -o ${vcfName}.gz ${vcfName} &&  bcftools index ${vcfName}.gz
    """
}

process extractRegions {
    publishDir "${baseOutdir}/isec_regions", mode: 'copy'
    
    input:
        path vcf
        each regions

    output:
        path "*isec*", emit: isec_vcfs
    script:
        def outBed = vcf.name.replace("vcf.gz", "isec.bed")
        def outVcf = vcf.name.replace("vcf.gz", "isec.vcf.gz")

        """
            intersectBed -a ${regions} -b ${vcf} -wb  > ${outBed}
            cut -f7 ${outBed} > id.lst
            bcftools filter -i 'ID=@id.lst' ${vcf} | bcftools sort -T ${params.tmpDir} -Oz -o ${outVcf}
        """
}

process goldCompare {
    publishDir "${baseOutdir}/gold_compare", mode: 'copy'
    
    input:
        path calls
        path gold_set

    output:
        path "*gold.shared.vcf.gz", emit: gold_shared_vcfs
    script:
        def outVcf = calls.name.replaceAll("highconf.*.vcf.gz", "gold.shared.vcf")

        """
            SURVIVOR merge <(ls ${calls} ${gold_set}) ${params.isecDist}\
            ${params.callerSupport} ${sameStrand} ${sameType} ${estDist}\
            ${params.isecMinLength} ${outVcf}
            bcftools sort -Oz  -o ${outVcf}.gz ${outVcf} &&  bcftools index ${outVcf}.gz
        """
}



workflow {
    if (params.sra) {
        reads = Channel.fromSRA(params.sra)
    } else {
        reads = Channel.fromPath(params.reads).ifEmpty { error "Please specify FASTQ files to use for variant calling" } 
    }

    ref = Channel.fromPath(params.ref).ifEmpty { error "Please specify a reference file to use for alignment (can be remote file https/ftp)" } 
    regions = params.regions ? Channel.fromPath(params.regions) : Channel.from("NO_FILE")
    alignReads(reads, ref)
    svimCalls(alignReads.out.alignment, ref, regions)
    snifflesCalls(alignReads.out.alignment, ref)
    svimCalls.out.svim_vcf | map { it.findAll { it =~/filtered.vcf.gz$/ }} | set { svim_filtered } 
    snifflesCalls.out.sniffles_vcf | map { it.findAll { it =~/raw.vcf.gz$/ }} | set { sniffles_calls } 
    highConfCalls(svim_filtered, sniffles_calls)

    if (params.regions) {
        svimCalls.out.svim_vcf | merge(snifflesCalls.out.sniffles_vcf) | merge(highConfCalls.out.highconf_vcf) | flatten |  filter { it =~ /(highconf|bnd|filtered|raw).vcf.gz$/ } | set { vcfs }
        extractRegions(vcfs, Channel.fromPath(params.regions))
    }

    if (params.goldSet) {
        if (params.regions) { 
            extractRegions.out.isec_vcfs | flatten() | filter { ~/highconf.isec.vcf.gz$/ } | set { gold_comp_vcf }
        } else {
            highConfCalls.out.highconf_vcf | set { gold_comp_vcf }
        }
        goldSetShared(gold_comp_vcf, Channel.fromPath(params.goldSet))
        
    }
}
