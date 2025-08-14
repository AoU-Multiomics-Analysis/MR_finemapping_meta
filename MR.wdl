task MR {
    File MungedSumstats 
    File SusieFinemapping
    String OutputPrefix
    Int Memory
    Int QTLSampleSize
    String QTLGroup

command{
    
    Rscript /mendelian_randomization.R \
        --MungedSumstats ${MungedSumstats} \
        --SusieFinemapping ${SusieFinemapping} \
        --OutputPrefix ${OutputPrefix} \
        --QTLSampleSize ${QTLSampleSize} \
        --QTLGroup ${QTLGroup}
        } 

runtime {
        docker: 'ghcr.io/aou-multiomics-analysis/mr_finemapping_meta:main'        
        memory: "${Memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "1"
        zones: ["us-central1-c"]
    }

output {
    File MR_output = "${OutputPrefix}_MR.tsv" 
    #File MR_enrichr = "${OutputPrefix}_MR_enrich.tsv"
    #File MR_LOO = "{OutputPrefix}_MR_LOO.tsv"

    }
}

workflow mendelian_randomization_workflow {
    call MR
}
