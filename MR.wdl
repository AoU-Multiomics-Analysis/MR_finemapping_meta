task MR {
    File MungedSumstats 
    File SusieFinemapping
    File MR_script
    String OutputPrefix
    Int Memory
    Int QTLSampleSize
    String QTLGroup

command{
    
    Rscript ${MR_script} \
        --MungedSumstats ${MungedSumstats} \
        --SusieFinemapping ${SusieFinemapping} \
        --OutputPrefix ${OutputPrefix} \
        --QTLSampleSize ${QTLSampleSize} \
        --QTLGroup ${QTLGroup}
        } 

runtime {
        docker: 'quay.io/kfkf33/susier:v24.01.1'        
        memory: "${Memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        cpu: "1"
        zones: ["us-central1-c"]
    }

output {
    File MR_output = "${OutputPrefix}_MR.tsv" 
    File MR_enrichr = "${OutputPrefix}_MR_enrich.tsv"

    }
}

workflow mendelian_randomization_workflow {
    call MR
}
