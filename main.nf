nextflow.enable.dsl = 2

include { extractVcf } from './modules/extract_vcf'
include { extractHap } from './modules/extract_hap'
include { runDashM } from './modules/run_dashm'

workflow {
    Channel.fromPath(params.sample_sheet, checkIfExists: true)
        | splitCsv(header: true, sep: '\t')
        | map { row ->
            [
                row.familyID,
                row.father,
                row.mother,
                file(row.jointvcf, checkIfExists: true),
                file(row.phasedvcf, checkIfExists: true),
                row.chrom,
                row.position,
                row.sample,
                row.CMtype,
                row.parent,
                row.sex
            ]
          }
        | set { samples }
    
    println "baseDir = $baseDir"

    dashm_scripts_ch = Channel.fromPath("${baseDir}/bin/*.{py,R}", checkIfExists: true)
    common_vcf_ch = Channel.fromPath("${baseDir}/data/00-common_all.vcf.txt.gz", checkIfExists: true)
    recom_prob_ch = Channel.fromPath("${baseDir}/Recom_Prob", type: 'dir') 



    extractVcf(samples)
    extractHap(samples)


    vcf_keyed = extractVcf.out.map { familyID, chrom, vcf_file ->
        def key = "${familyID}_${chrom}" // 
        [key, [familyID, chrom, vcf_file]]
    }

    hap_keyed = extractHap.out.map { familyID, chrom, hap_file ->
        def key = "${familyID}_${chrom}"
        [key, [familyID, chrom, hap_file]]
    }

    samples_keyed = samples.map { familyID, father, mother, jointvcf, phasedvcf, chrom, position, sample, CMtype, parent, sex ->
        def key = "${familyID}_${chrom}"
        [key, [familyID, father, mother, position, sample, CMtype, parent, sex]]
    }


    joined_data = vcf_keyed
        .join(hap_keyed)
        .join(samples_keyed)
        .map { key, vcf_data, hap_data, sample_data ->
            def familyID = vcf_data[0]
            def chrom = vcf_data[1]
            def vcf_file = vcf_data[2]
            def hap_file = hap_data[2]
            def father = sample_data[1]
            def mother = sample_data[2]
            def position = sample_data[3]
            def sample = sample_data[4]
            def CMtype = sample_data[5]
            def parent = sample_data[6]
            def sex = sample_data[7]

            [
                familyID,    // familyID
                chrom,       // chrom
                vcf_file,    // vcf_txt
                hap_file,    // hap_txt
                position,    // position
                father,      // father
                mother,      // mother
                sample,      // sample
                CMtype,      // CMtype
                parent,      // parent
                sex          // sex
            ]
        }


    joined_data.ifEmpty { 
        println "WARNING: joined_data is empty! Check if all channels have matching keys."
    }

    runDashM(joined_data,common_vcf_ch.collect(),recom_prob_ch.collect())

}
