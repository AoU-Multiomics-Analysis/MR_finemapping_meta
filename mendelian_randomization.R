message('MR start')

library(R.utils)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(enrichR)
library(optparse)

message('Libraries loaded')

#install.packages('R.utils')
theme_set(theme_classic())
theme_update(panel.border = element_rect(fill = NA,linewidth = .9), 
            axis.line = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 20))

# set dbs for enrichr
dbs <- c('GO_Biological_Process_2025',
         'GO_Cellular_Component_2025',
         'GO_Molecular_Function_2025',
         'MSigDB_Hallmark_2020',
         'KEGG_2021_Human',
         'Reactome_Pathways_2024'
        )



############ FUNCTIONS ###########
# helper function to load susie parquet files 
load_finemapping_data <- function(path){
split_name <- str_split(basename(path),'_|\\.') %>% unlist()
group <- split_name[3]
    
fm_data <- arrow::read_parquet(path) %>% 
    #separate(variant_id,into = c('chrom','pos','alt')) %>% 
    #extract(pos, into = c("pos", "ref"), regex = "([0-9]+)([A-Za-z]+)") %>%
    mutate(group = group)
fm_data
    
}

# helper function to load munged sumstats 
# and format for MR 
load_gwas_data <- function(GWAS_path){
GWAS_dat <- fread(GWAS_path)
GWAS_dat_cols <- colnames(GWAS_dat)

message(paste0('GWAS columns: ',GWAS_dat_cols))

if (!'SE' %in% GWAS_dat_cols & 'OR' %in% GWAS_dat_cols){
message('SE measurement is missing, computing from OR and P value')
GWAS_dat$SE <- get_se(GWAS_dat$OR,GWAS_dat$P)
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'OR')

if (!'FRQ' %in% GWAS_dat_cols){
GWAS_dat$FRQ <- NA

}

    
} else if (!'SE' %in% GWAS_dat_cols & 'BETA' %in% GWAS_dat_cols){
messasge('SE measurement is missing, computing from BETA and P value')
GWAS_dat$SE <- TwoSampleMR::get_eff(GWAS_dat$BETA,GWAS_dat$P)
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'BETA')

} else if ('SE' %in% GWAS_dat_cols & 'BETA' %in% GWAS_dat_cols){
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'BETA')
}else if ('SE' %in% GWAS_dat_cols & 'OR' %in% GWAS_dat_cols){
GWAS_dat <- rename(GWAS_dat,'beta.outcome' = 'OR')
}

GWAS_out <- GWAS_dat %>% 
    dplyr::rename( 'effect_allele.outcome' = 'A1',
                   'other_allele.outcome' = 'A2',
                   'id.outcome' = 'outcome',
                   'eaf.outcome' = 'FRQ',
                    'se.outcome' = 'SE'
                    ) %>% 
    mutate(outcome = 'GWAS',
           VARIANT_POSITION = paste0(CHR,'_',BP),
           VARIANT_ID = paste(CHR,BP,effect_allele.outcome,other_allele.outcome,sep ='_')) %>% 
    mutate(RSID = SNP,SNP = VARIANT_ID)
GWAS_out
}
filter_variants_adjust_pips <- function(fm_data,variant_list){

# filter to variants in the GWAS data and adjust 
# pips to account for removed variants 
QTL_GWAS_variants <- fm_data %>% 
    filter(variant %in%variant_list) %>% 
    group_by(cs_id) %>% 
    filter(sum(pip) > 0.5) %>% 
    mutate(adjusted_pip = (pip/sum(pip) *0.95 ) )  %>% 
    ungroup()

QTL_GWAS_variants  
}


calc_I2 <- function(Q, Est) {
  Q <- Q[[1]]
  Est <- length(unique(Est))
  I2 <- if (Q > 1e-3) (Q - Est + 1) / Q else 0
  return(if (I2 < 0) 0 else I2)
}

create_MR_input <- function(harmonised_data){
MR_input <- harmonised_data %>%  
    mutate(
          bhat_x = beta.exposure / se.exposure,
          sbhat_x = 1
        )  %>% 
    group_by(molecular_trait_id, cs_id) %>%
    mutate(cpip = sum(adjusted_pip),pip = adjusted_pip)  %>% 
    dplyr::rename('bhat_y' = 'beta.outcome','sbhat_y' = 'se.outcome')  %>% 
    mutate(
        beta_yx = bhat_y / bhat_x,
        se_yx = sqrt((sbhat_y^2 / bhat_x^2) + ((bhat_y^2 * sbhat_x^2) / bhat_x^4)),
        composite_bhat = sum((beta_yx * pip) / cpip),
        composite_sbhat = sum((beta_yx^2 + se_yx^2) * pip / cpip)
      ) %>%
    ungroup()
    
MR_input   
}

format_fm_data_MR <- function(fm_data,sample_size){
cleaned_fm_data <- fm_data %>%  
    mutate(
           se.exposure = posterior_sd/sqrt(sample_size), 
           eaf.exposure = NA,
           SNP = variant,
           variant = variant,
           exposure = 'QTL',
           id.exposure = molecular_trait_id
          ) %>% 
    dplyr::rename(
        'beta.exposure' = 'posterior_mean',
        'effect_allele.exposure' = 'ref',
        'other_allele.exposure' = 'alt')  %>% 
     mutate(
           se.exposure = posterior_sd/sqrt(sample_size), 
           eaf.exposure = NA,
           SNP = variant,
           variant = variant,
           exposure = 'QTL',
           id.exposure = molecular_trait_id
          ) 
cleaned_fm_data     
}

extract_gene_set <- function(MR_output){
require(gprofiler2)

gene_set <- MR_output %>% 
        mutate(gene_id = str_remove(molecular_trait_id,'.*_')) %>% 
        mutate(gene_id = str_remove(gene_id,'\\..*')) %>% 
        select(gene_id)
gene.symbols <- gconvert(gene_set$gene_id,organism="hsapiens",target="ENTREZGENE",filter_na = F) %>% 
            select(input,name) %>% 
            distinct()
gene.symbols
}

####### PARSE COMMAND LINE ARGUMENTS #######
message('Begin analaysis')
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--MungedSumstats"), type="character", default=NULL,
                        help="Path to GWAS summary statistics that have been munged", metavar = "type"),
  optparse::make_option(c("--SusieFinemapping"), type="character", default=NULL,
                        help="Susie finemapping parquet", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="output prefix for files", metavar = "type"),
    optparse::make_option(c("--QTLSampleSize"), type="character", default=NULL,
                        help="", metavar = "type"),
    optparse::make_option(c("--QTLGroup"), type="character", default=NULL,
                        help="", metavar = "type")
  )

message('Parsinging command line arguments')
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

group <- opt$QTLGroup
output_prefix <- opt$OutputPrefix

enrichr_output <- paste0(output_prefix,'_MR.tsv')
MR_output_file <- paste0(output_prefix,'_MR_enrich.tsv')
message(paste0('Writing MR results to',MR_output_file ))


GWAS_path <- opt$MungedSumstats
fm_path <- opt$SusieFinemapping
message(paste0('Using ',GWAS_path , ' for GWAS'))
message(paste0('Using ',fm_path , ' for QTL'))
message(paste0('QTL group: ',group ))





########## LOAD DATA #########
message('Loading summary stats')
GWAS_dat <- load_gwas_data(GWAS_path)

message('Loading QTL finemapping')
fm_data <- load_finemapping_data(fm_path) %>% 
    mutate(variant = str_remove_all(variant,'chr'))

message('Cleaning finemapping data and adjust pips')
# load finemapping data and compute standard error
cleaned_fm_data <- fm_data %>% 
    format_fm_data_MR(sample_size = as.integer(opt$QTLSampleSize))  %>% 
    filter_variants_adjust_pips(GWAS_dat$VARIANT_ID)

message('Merging data')
harmonised_dat <- harmonise_data(
    exposure_dat = cleaned_fm_data %>%  select(-SNP) %>% dplyr::rename('SNP' = 'variant'),
    outcome_dat = GWAS_dat,
    action = 1
)


message('Creating MR input')
MR_input <- harmonised_dat %>%  
    create_MR_input()

###### RUN MR #######

message('Running MR')
MR_output <- MR_input %>%
    group_by(molecular_trait_id) %>% 
    mutate(
        composite_sbhat = sqrt(composite_sbhat - composite_bhat^2),
        wv = composite_sbhat^-2
      )  %>% 
    mutate(
        meta_eff = sum(unique(wv) * unique(composite_bhat)),
        sum_w = sum(unique(wv)),
        se_meta_eff = sqrt(sum_w^-1),
        num_CS = length(unique(cs_id))
      )  %>% 
      mutate(
        num_IV = length(SNP),
        meta_eff = meta_eff / sum_w,
        meta_pval = 2 * pnorm(abs(meta_eff) / se_meta_eff, lower.tail = FALSE),
        Q = sum(unique(wv) * (unique(composite_bhat) - unique(meta_eff))^2),
        I2 = calc_I2(Q, composite_bhat),
        Q_pval = pchisq(Q, df = length(unique(composite_bhat)) - 1, lower = F)
      ) %>% 
      ungroup() %>%
      distinct(molecular_trait_id, .keep_all = TRUE) %>%
      mutate(
        cpip = round(cpip, 3),
        meta_pval = round(meta_pval, 10),
        meta_eff = round(meta_eff, 5),
        se_meta_eff = round(se_meta_eff, 3),
        Q = round(Q, 3),
        Q_pval = round(Q_pval, 3),
        I2 = round(I2, 3)
      ) %>%
      arrange(meta_pval) %>% 
      select(molecular_trait_id, num_CS, num_IV, cpip, meta_eff, se_meta_eff, meta_pval, Q, Q_pval, I2)


message('Writing MR results  to output')
MR_output %>% mutate(group = group)%>% write_tsv(MR_output_file)

####### RUN GENE SET ENRICHMENT ######
message('Extracting significant genes for GSEA')
sig_MR_genes <- MR_output %>%  
    arrange(desc(abs(meta_eff))) %>%
    mutate(padj = p.adjust(meta_pval,method = 'fdr')) %>% 
    filter(padj < .1 & Q_pval > .05)
gene_list <- sig_MR_genes %>% extract_gene_set

message('Running enrichment analysis')
enrichr_res <- enrichr(gene_list$name, dbs) %>%
    imap_dfr(~bind_rows(.) %>% mutate(library = .y)) %>% 
    arrange(Adjusted.P.value)

message('Writing enrichment results to output')
enrichr_res %>% write_tsv(enrichr_output)





