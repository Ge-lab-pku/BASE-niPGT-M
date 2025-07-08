nextflow.enable.dsl = 2

process extractHap {
    tag "${familyID}_${chrom}"
    publishDir "results/hap", mode: 'copy' , enabled: false
    
    input:
    tuple val(familyID), val(father), val(mother), path(jointvcf), path(phasedvcf), 
          val(chrom), val(position), val(sample), val(CMtype), val(parent), val(sex)
    
    output:
    tuple val(familyID), val(chrom), path("${familyID}_chr*.hap.txt")
    
    script:
    """
    output_chrom="chr\$(echo $chrom | sed 's/^chr//')"

    tabix -p vcf $phasedvcf
    echo -e "CHROM\\tPOS\\tREF\\tALT\\t${father}\\t$mother" > ${familyID}_\${output_chrom}.hap.txt
    bcftools view -r \$output_chrom -s $father,$mother $phasedvcf | \
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\n' \
    >> ${familyID}_\$output_chrom.hap.txt
    """
}
