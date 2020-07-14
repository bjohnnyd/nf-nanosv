/*
Copyright (c) 2020 Bisrat J. Debebe

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// TODO: survivor can be specified as `params.both.asInteger()`
nextflow.preview.dsl = 2


def baseOutdir = params.name ? "${params.outdir}/${params.name}" : "${params.outdir}"
String.metaClass.isEmpty = { delegate.allWhitespace }

process alignReads {
    publishDir  "${baseOutdir}/alignment", mode: 'copy'
    input:
        path reads
        path ref
    output:
        path '*alignment.sort.bam', emit: alignment
    script:
    def bamName = params.name ? "${params.name}.alignment.sort.bam" :"alignment.sort.bam" 
    """
        minimap2 --cs --MD -t ${task.cpus} -ax map-ont ${ref} ${reads} | samtools sort -o ${bamName} - 
    """

}

process getLowMapQ {
    publishDir  "${baseOutdir}/alignment/lowmq", mode: 'copy'
    input:
        path alignment
    output:
        path '*lowmq*'
    script:
    def outName = params.name ? "${params.name}.lowmq" :"lowmq" 
    """
        samtools view -h ${alignment} | awk '\$1 ~ /^@/ || \$5<5 '  | samtools sort -o ${outName}.sort.bam -
        samtools depth ${outName}.sort.bam > lowq.depth
        SURVIVOR bincov lowq.depth ${params.windowSize} ${params.lowQualSupport} > ${outName}.bed
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
        path "*bnd*", emit: bnd_results
    script:
        def (sampleCmd, outRawVcf, outBndVcf, outFilteredVcf, outBndInfo, outBndBam, outBndBedPe) = params.name ? 
            [ "--sample ${params.name}", "${params.name}.svim.raw.vcf.gz",
            "${params.name}.svim.bnd.vcf.gz", "${params.name}.svim.filtered.vcf.gz",
            "${params.name}.bnd.lst", "${params.name}.bnd.sort.bam", "${params.name}.bnd.sort.bedpe"
            ] 
            : 
            [ "", "svim.raw.vcf.gz", "svim.bnd.vcf.gz", "svim.filtered.vcf.gz",
            "bnd.lst", "bnd.sort.bam",  "bnd.sort.bedpe"
            ]


        """
            svim alignment --read_names ${sampleCmd} svim_calls ${bam} ${ref}
            tar czf svim.results.tar.gz svim_calls 
            mv svim_calls/*{log,png} .
            bcftools sort -Oz -o ${outRawVcf} svim_calls/variants.vcf 
            bcftools filter -i 'FILTER=="PASS" && SVTYPE=="BND" && SUPPORT >= 1 && QUAL >= 1' ${outRawVcf} | bcftools sort -Oz  -o ${outBndVcf} - 
            bcftools filter -i "${params.svimFilter}" ${outRawVcf} | bcftools sort -Oz -o ${outFilteredVcf}  -
            bcftools query -f '%ID\\t%CHROM\\t%POS\\t%ALT\\t%READS\\t%QUAL\n' ${outBndVcf} | sort -k6nr > ${outBndInfo}
            cut -f5 ${outBndInfo} | tr ',' '\\n' > reads.lst
            samtools view -H ${bam} > header.sam
            samtools view ${bam} | fgrep -w -f reads.lst > alignments.sam
            cat header.sam alignments.sam | samtools sort -o ${outBndBam} -
            bamToBed -bedpe -cigar -i ${outBndBam} > ${outBndBedPe}
        """

}

        
process snifflesCalls {
    publishDir  "${baseOutdir}/sv/sniffles_calls", mode: 'copy'
    input:
        path bam
        path ref
    output:
        path '*sniffles.raw.vcf.gz*', emit: sniffles_vcf
    script:
    def vcfName = params.name ? "${params.name}.sniffles.raw.vcf" : "sniffles.raw.vcf"
    def cluster = params.snifflesCluster ? "--cluster" : ""
    def genotype = params.snifflesGenotype ? "--genotype" : ""


    """ 
        sniffles -t ${task.cpus} -m ${bam} -s ${params.snifflesMinSupport} \
        --cluster_support ${params.snifflesClusterSupport} -v ${vcfName}  \
        -l ${params.snifflesMinLength} -r ${params.snifflesMinSeqSize} \
        ${cluster} ${genotype} ${params.snifflesAdvanced}
        bcftools sort -Oz -o ${vcfName}.gz ${vcfName} && tabix ${vcfName}.gz
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
        bcftools view -Ov ${svim_calls} > svim_calls.vcf
        bcftools view -Ov ${sniffles_calls} > sniffles_calls.vcf
        SURVIVOR merge <(ls svim_calls.vcf sniffles_calls.vcf) ${params.isecDist}\
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
            bcftools filter -i 'ID=@id.lst' ${vcf} | bcftools sort -Oz -o ${outVcf}
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
        def sameStrand = params.sameStrand ? 1 : 0
        def sameType = params.sameType ? 1 : 0
        def estDist = params.estDist ? 1 : 0

        """
            bcftools view -Ov ${calls} > calls.vcf
            bcftools view -Ov ${gold_set} > gold_set.vcf
            SURVIVOR merge <(ls calls.vcf gold_set.vcf) ${params.isecDist}\
            ${params.callerSupport} ${sameStrand} ${sameType} ${estDist}\
            ${params.isecMinLength} ${outVcf}
            bcftools sort -Oz  -o ${outVcf}.gz ${outVcf} &&  bcftools index ${outVcf}.gz
        """
}

process plotTopQual {
    publishDir "${baseOutdir}/plots/top_qual", mode: 'copy'
    
    input:
        path calls
        path alignment
        val topN

    output:
        path "*.png", emit: topQualPlots
    script:
        def runName = calls.name.replaceAll("highconf.vcf.gz", "")
        def outName = runName.trim().length() == 0 ? "" : "${runName}"

        """
            samtools index ${alignment}
            bcftools query -f '%ID\\t%CHROM\\t%POS\\t%SVLEN\\t%QUAL\\t%ALT\\n' ${calls} | sort -k5nr | head -${topN} | awk 'BEGIN {OFS="\\t";} {gsub("-","",\$4);gsub("(<|>)","",\$6);print \$1,\$2,\$3,\$3+\$4,\$4,\$5,\$6,FNR}' > variants.lst
            parallel -j ${task.cpus} -C '\\t' "samplot plot -n '{1}-Q{5}' -b ${alignment} -o ${outName}{7}.png -c {2} -s {3} -e {4} -t {6}"  :::: <(cut -f1-4,6-8 variants.lst)
        """
}
    

//TODO: Svim BND reads --> filter alignment bam --> bamtobedPe


workflow {
    if (params.sra) {
        reads = Channel.fromSRA(params.sra)
    } else {
        if (!params.reads) {
            error "Please provide path to FASTQ with --reads <path/to/fastq>"
        } else {
            reads = Channel.fromPath(params.reads).ifEmpty { error "Please specify FASTQ files to use for variant calling" } 
        }
    }

    if (!params.name) {
        println "You have not provided a name for this run please make sure the output is unique as the previous run will be overwritten"
    }

    ref = Channel.fromPath(params.ref).ifEmpty { error "Please specify a reference file to use for alignment (can be remote file https/ftp)" } 
    regions = params.regions ? Channel.fromPath(params.regions) : Channel.from("NO_FILE")
    alignReads(reads, ref)
    getLowMapQ(alignReads.out.alignment)
    svimCalls(alignReads.out.alignment, ref, regions)
    snifflesCalls(alignReads.out.alignment, ref)
    svimCalls.out.svim_vcf | map { it.findAll { it =~/filtered.vcf.gz$/ }} | set { svim_filtered } 
    snifflesCalls.out.sniffles_vcf | map { it.findAll { it =~/raw.vcf.gz$/ }} | set { sniffles_calls } 
    highConfCalls(svim_filtered, sniffles_calls)
    plotTopQual(highConfCalls.out.highconf_vcf, alignReads.out.alignment, params.topN)

    if (params.regions) {
        svimCalls.out.svim_vcf | merge(snifflesCalls.out.sniffles_vcf) | merge(highConfCalls.out.highconf_vcf) | flatten |  filter { it =~ /(highconf|bnd|filtered|raw).vcf.gz$/ } | set { vcfs }
        extractRegions(vcfs, Channel.fromPath(params.regions))
    }

    if (params.goldSet) {
        if (params.regions) { 
            extractRegions.out.isec_vcfs | flatten() | filter { it =~/highconf.isec.vcf.gz$/ } | set { gold_comp_vcf }
        } else {
            highConfCalls.out.highconf_vcf | set { gold_comp_vcf }
        }
        goldCompare(gold_comp_vcf, Channel.fromPath(params.goldSet))
        
    }
}
