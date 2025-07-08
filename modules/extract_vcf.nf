nextflow.enable.dsl = 2

process extractVcf {
    tag "${familyID}_${chrom}"
    publishDir "results/vcf", mode: 'copy', pattern: "*.{vcf.gz,vcf.gz.tbi,txt}" , enabled: false
    
    input:
    tuple val(familyID), val(father), val(mother), path(jointvcf), val(phasedvcf), 
          val(chrom), val(position), val(sample), val(CMtype), val(parent), val(sex)
    
    output:
    tuple val(familyID), val(chrom), path("${familyID}_chr*.vcf.txt")
    
    script:
    """
    output_chrom="chr\$(echo $chrom | sed 's/^chr//')"

    tabix -p vcf $jointvcf
    bcftools view -r \$output_chrom $jointvcf -O z -o ${familyID}_\$output_chrom.vcf.gz
    tabix -p vcf ${familyID}_\$output_chrom.vcf.gz
    
    Rscript ${projectDir}/bin/s1.prepare_vcf.R ${familyID}_\$output_chrom.vcf.gz
    """
}
