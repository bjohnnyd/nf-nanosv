/*
Copyright (c) 2020 Bisrat J. Debebe

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

nextflow.preview.dsl = 2

def baseOutdir = params.name ? "${params.outdir}/${params.name}" : "${params.outdir}"
String.metaClass.isEmpty = { delegate.allWhitespace }

process alignReads {
    publishDir  "${baseOutdir}/alignment", mode: 'copy'
    input:
        path reads
        path ref
    output:
        tuple path('*alignment.sort.bam'), path('*alignment.sort.bam.bai'), emit: alignment
    script:
    def bamName = params.name ? "${params.name}.alignment.sort.bam" :"alignment.sort.bam" 
    """
        minimap2 --cs --MD -t ${task.cpus} -ax map-ont ${ref} ${reads} | samtools sort -o ${bamName} - 
        samtools index ${bamName}
    """

}

process getLowMapQ {
    publishDir  "${baseOutdir}/alignment/lowmq", mode: 'copy'
    input:
        tuple path(alignment), path(alignment_index)
    output:
        path '*lowmq*', emit: lowmq
        path '*lowmq.bed', emit: lowmq_bed
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
        tuple path(bam), path(bam_index), val(dist)
        path ref
        file regions
    output:
        path 'svim.results.tar.gz', emit: svim_dir
        path "**.vcf.gz*", emit: svim_vcf
        path "*.log", emit: svim_log
        path "*.png", emit: svim_images
        path "*bnd*", emit: bnd_results
        tuple val(dist), path("*svim.solo${dist}.filtered.vcf.gz"),  emit: svim_filtered_vcf
    script:
        def (sampleCmd, fullPrefix, bndPrefix) = params.name ?  [ "--sample ${params.name}", "${params.name}.svim",  "${params.name}.bnd"] : ["", "svim", "bnd"]
        def clusterDist = "${dist}"

        if(clusterDist == "0" || clusterDist == "None") {
            def filterCmd = """
                bcftools sort -Oz -o ${fullPrefix}.solo${clusterDist}.vcf.gz ${fullPrefix}.filtered.vcf.gz
                bcftools view -h -Oz -o ${fullPrefix}.clustered${clusterDist}.vcf.gz  ${fullPrefix}.filtered.vcf.gz
                tabix ${fullPrefix}.solo${clusterDist}.vcf.gz && tabix ${fullPrefix}.clustered${clusterDist}.vcf.gz

            """
        } else {
            def filterCmd = """
                bcftools sort -Ov variants.vcf > calls.vcf
                for i in {0..1};do echo "calls.vcf" >> vcf.lst
                SURVIVOR merge vcf.lst ${clusterDist} 2 0 0 0 ${params.soloVarSizeMin} clustered_vars.vcf
                bedtools intersect -v -wa -a calls.vcf -b clustered_vars.vcf > solo_vars.vcf
                cat <(bcftools view -h calls.vcf) solo_vars.vcf | bcftools sort -Oz -o ${fullPrefix}.solo${clusterDist}.vcf.gz -
                bcftools filter -i "${params.svimFilter}" ${fullPrefix}.solo${clusterDist}.vcf.gz | bcftools sort -Oz -o ${fullPrefix}.solo${clusterDist}.filtered.vcf.gz  -
                bcftools sort -Oz -o ${fullPrefix}.clustered${clusterDist}.vcf.gz clustered_vars.vcf
                tabix ${fullPrefix}.solo${clusterDist}.vcf.gz && tabix ${fullPrefix}.clustered${clusterDist}.vcf.gz
            """
        }

        """
            svim alignment --read_names ${sampleCmd} svim_calls ${bam} ${ref}
            tar czf svim.results.tar.gz svim_calls 
            mv svim_calls/*{log,png} .
            bcftools view -Oz -o ${fullPrefix}.raw.vcf.gz svim_calls/variants.vcf 
            bcftools filter -i 'FILTER=="PASS" && SVTYPE=="BND" && SUPPORT >= ${params.bndSupport} && QUAL >= ${params.bndQual}' ${fullPrefix}.raw.vcf.gz | bcftools sort -Oz  -o ${fullPrefix}.bnd.vcf.gz - 
            bcftools filter -i "${params.svimFilter}" ${fullPrefix}.raw.vcf.gz | bcftools sort -Oz -o ${fullPrefix}.filtered.vcf.gz  -
            bcftools query -f '%ID\\t%CHROM\\t%POS\\t%ALT\\t%READS\\t%QUAL\\n' ${fullPrefix}.bnd.vcf.gz | sort -k6nr > ${bndPrefix}.lst
            cut -f5 ${bndPrefix}.lst | tr ',' '\\n' > reads.lst
            samtools view -H ${bam} > header.sam
            samtools view ${bam} | fgrep -w -f reads.lst > alignments.sam || :
            cat header.sam alignments.sam | samtools sort -o ${bndPrefix}.sort.bam -
            bamToBed -bedpe -cigar -i ${bndPrefix}.sort.bam > ${bndPrefix}.sort.bedpe
            ${filterCmd}
        """

}

        
process snifflesCalls {
    label 'cpu'

    publishDir  "${baseOutdir}/sv/sniffles_calls", mode: 'copy'
    input:
        tuple path(bam), path(bam_index), val(dist)
        path ref
    output:
        path "*vcf.gz"
        path '*sniffles.raw.vcf.gz', emit: sniffles_raw_vcf
        tuple val(dist), path("*sniffles.solo${dist}.filtered.vcf.gz"),  emit: sniffles_filtered_vcf
    script:
    def vcfName = params.name ? "${params.name}.sniffles" : "sniffles"
    def cluster = params.snifflesCluster ? "--cluster" : ""
    def genotype = params.snifflesGenotype ? "--genotype" : ""
    def clusterDist = "${dist}"

    if(clusterDist == "0" || clusterDist == "None") {
        def filterCmd = 
        """
            bcftools sort -Oz -o ${vcfName}.solo${clusterDist}.filtered.vcf.gz ${vcfName}.raw.vcf.gz
            bcftools view -h -Oz -o ${vcfName}.clustered${clusterDist}.vcf.gz  ${vcfName}.filtered.vcf.gz
            tabix ${vcfName}.solo${clusterDist}.filtered.vcf.gz && tabix ${vcfName}.clustered${clusterDist}.vcf.gz

        """
    } else {
        def filterCmd = 
        """
            bcftools view -Ov ${vcfName}.raw.vcf.gz > calls.vcf
            for i in {0..1};do echo "calls.vcf" >> vcf.lst
            SURVIVOR merge vcf.lst ${clusterDist} 2 0 0 0 ${params.soloVarSizeMin} clustered_vars.vcf
            bedtools intersect -v -wa -a calls.vcf -b clustered_vars.vcf > solo_vars.vcf
            cat <(bcftools view -h calls.vcf) solo_vars.vcf | bcftools sort -Oz -o ${vcfName}.solo${clusterDist}.vcf.gz -
            bcftools filter -i "${params.svimFilter}" ${vcfName}.solo${clusterDist}.vcf.gz | bcftools sort -Oz -o ${vcfName}.solo${clusterDist}.filtered.vcf.gz  -
            bcftools sort -Oz -o ${vcfName}.clustered${clusterDist}.vcf.gz clustered_vars.vcf
            tabix ${vcfName}.solo${clusterDist}.filtered.vcf.gz && tabix ${vcfName}.clustered${clusterDist}.vcf.gz
        """
    }


    """ 
        sniffles -t ${task.cpus} -m ${bam} -s ${params.snifflesMinSupport} \
        --cluster_support ${params.snifflesClusterSupport} -v ${vcfName}.raw.vcf  \
        -l ${params.snifflesMinLength} -r ${params.snifflesMinSeqSize} \
        ${cluster} ${genotype} ${params.snifflesAdvanced}
        bcftools sort -Oz -o ${vcfName}.raw.vcf.gz ${vcfName}.raw.vcf
    """

}

process highConfCalls {
    publishDir  "${baseOutdir}/sv/high_conf", mode: 'copy'
    input:
        tuple path(svim_calls), path(sniffles_calls)
        each path(lowmq_bed)
    output:
        path '*highconf*vcf.gz*', emit: highconf
    script:
    def vcfName = params.name ? "${params.name}.highconf" :"highconf" 
    def sameStrand = params.sameStrand ? 1 : 0
    def sameType = params.sameType ? 1 : 0
    def estDist = params.estDist ? 1 : 0

    """ 
        bcftools view -Ov ${svim_calls} > svim_calls.vcf
        bcftools view -Ov ${sniffles_calls} > sniffles_calls.vcf
        SURVIVOR merge <(ls svim_calls.vcf sniffles_calls.vcf) ${params.isecDist}\
        ${params.callerSupport} ${sameStrand} ${sameType} ${estDist}\
        ${params.isecMinLength} ${vcfName}.vcf

        bcftools sort -Oz  -o ${vcfName}.vcf.gz ${vcfName}.vcf
        bcftools index ${vcfName}.vcf.gz && bcftools tabix ${vcfName}.vcf.gz
        intersectBed -v -a ${vcfName}.vcf.gz -b ${lowmq_bed} > highconf.nolowmq.body
        cat <(bcftools view -h ${vcfName}.vcf.gz) highconf.nolowmq.body | bcftools sort -Oz -o ${vcfName}.nolowmq.vcf.gz -

        bcftools filter -i 'SVTYPE == "DEL" || SVTYPE == "INS"' ${vcfName}.vcf.gz | bcftools sort -Oz -o ${vcfName}.indel.vcf.gz -
        bcftools filter -i 'SVTYPE == "DEL" || SVTYPE == "INS"' ${vcfName}.nolowmq.vcf.gz | bcftools sort -Oz -o ${vcfName}.nolowmq.indel.vcf.gz -
        bcftools index ${vcfName}.indel.vcf.gz && tabix ${vcfName}.indel.vcf.gz
        bcftools index ${vcfName}.nolowmq.indel.vcf.gz && tabix ${vcfName}.nolowmq.indel.vcf.gz
    """
}

process extractRegions {
    label 'cpu'
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
            bcftools index ${outVcf} && tabix ${outVcf}
        """
}

process calculateCoverage {
    publishDir  "${baseOutdir}/alignment/coverage", mode: 'copy'
    
    input:
        tuple path(bam), path(bam_index)
        each regions

    output:
        path "./*{mosdepth,per-base,regions,quantized,thresholds}*.{txt,bed.gz}" optional true
    script:
        def quantizeCmd = params.quantCoverage ? "-q ${params.quantCoverage}" : "" 
        def regionsCmd = regions.name == "NO_FILE" ? "" : "-b ${regions}"
        def prefix = params.name ? "${params.name}" : "" 

        """
            mosdepth -t ${task.cpus} ${quantizeCmd} ${regionsCmd} ${prefix} ${bam}
        """
}

process goldCompare {
    publishDir "${baseOutdir}/gold_compare", mode: 'copy'
    
    input:
        path calls
        each path(gold_set)

    output:
        path "*gold.shared.vcf.gz", emit: gold_shared_vcfs
    script:
        def outVcf = calls.name.replaceAll("vcf.gz", "gold.shared.vcf")
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

process lowMapQualFilter {
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

process extractPlotRegions {
    input:
        path calls
        val topN

    output:
        path "variants.lst", emit: plot_regions
    shell:
        '''
            bcftools filter -e 'SVTYPE == "INS"' !{calls} | bcftools query -f '%ID\\t%CHROM\\t%POS\\t%SVLEN\\t%QUAL\\t%ALT\\n' | sort -k5nr | awk -v topN=!{topN}  'BEGIN {OFS="\\t";} FNR <= topN {gsub("-","",$4);gsub("(<|>)","",$6);print $1,$2,$3,$3+$4,$4,$5,$6,FNR}' > variants.lst
        '''
}
    
process plotTopQual {
    publishDir "${baseOutdir}/plots/top_qual", mode: 'copy'
    
    input:
        tuple path(alignment), path(alignment_index), val(variant_row)

    output:
        path "*.png", emit: topQualPlots
    script:
        def runName = alignment.name.replaceAll("alignment.sort.bam", "")
        def outName = runName.trim().length() == 0 ? "" : "${runName}"
        def(id, chrom, start, end, qual, type, rank) = variant_row.tokenize(',')

        """
            samplot plot -n '${id}-Q${qual}' -b ${alignment} -o ${outName}${rank}.png -c ${chrom} -s ${start} -e ${end} -t ${type}
        """
}

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
    dists = Channel.from(params.soloDist.split(","))

    // Main
    alignReads(reads, ref)
    getLowMapQ(alignReads.out.alignment)
    calculateCoverage(alignReads.out.alignment, regions) 
    alignReads.out.alignment | combine(dists) | set { align_dist_ch} 
    svimCalls(align_dist_ch, ref, regions)
    snifflesCalls(align_dist_ch, ref)
    svimCalls.out.svim_vcf | map { it.findAll { it =~/filtered.vcf.gz$/ }} | set { svim_filtered } 
    svimCalls.out.svim_filtered_vcf | join(snifflesCalls.out.sniffles_filtered_vcf) | map { tuple(it[1], it[2])} |set { var_dist_ch} 
    highConfCalls(var_dist_ch, getLowMapQ.out.lowmq_bed)

    // Optional
    if (params.regions) {
        svimCalls.out.svim_vcf | merge(snifflesCalls.out.sniffles_raw_vcf) | merge(highConfCalls.out.highconf) | flatten |  filter { it =~ /(highconf|bnd|filtered|sniffles|indel|nolowmq).vcf.gz$/ } | set { vcfs } 
        extractRegions(vcfs, Channel.fromPath(params.regions))
    }

    if (params.goldSet) {
        if (params.regions) { 
            extractRegions.out.isec_vcfs | flatten() | filter { it =~/isec.vcf.gz$/ } | set { gold_comp_vcf }
        } else {
            highConfCalls.out.highconf | set { gold_comp_vcf }
        }
        goldCompare(gold_comp_vcf, Channel.fromPath(params.goldSet))
        
    }

    // Plotting
    highConfCalls.out.highconf | flatten | filter { it =~ /highconf.vcf.gz$/ } | set {to_plot_vcf}
    extractPlotRegions(to_plot_vcf, params.topN)
    extractPlotRegions.out.plot_regions | splitCsv(sep: '\t') | map { it.removeAt(4); it.join(",") }  | set { plot_info }
    plotTopQual(alignReads.out.alignment.combine(plot_info))

}
