process runDashM {

    tag "${familyID}_${chrom}"

    publishDir "${params.outdir ?: 'results'}/${familyID}",
        pattern: 'temp_output/*_results.xlsx', 
        mode: 'copy',
        saveAs: { it.split('/')[-1] }

    publishDir "${params.outdir ?: 'results'}/${familyID}",
        pattern: 'temp_output/*.svg',
        mode: 'copy',
        saveAs: { it.split('/')[-1] }

    input:
    tuple val(familyID), val(chrom), path(vcf_txt), path(hap_txt), val(position), val(father), val(mother), val(sample), val(CMtype), val(parent), val(sex)
    path common_vcf  
    path recom_prob 

    output:
    path "temp_output/*_results.xlsx", emit: results, optional: true
    path "temp_output/*.svg", emit: plots, optional: true

    script:
    """
    ~/miniconda3/ATAC2/bin/python ${baseDir}/bin/DashM_Main.2.py \\
        --raw_data $vcf_txt \\
        --phased_data $hap_txt \\
        --chrom $chrom \\
        --position $position \\
        --father $father \\
        --mother $mother \\
        --sample $sample \\
        --cm_type $CMtype \\
        --parent $parent \\
        --sex $sex \\
        --output_dir temp_output
    """
}
