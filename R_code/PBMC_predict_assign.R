### Functions for assign the cell type predicted by GPTcelltype to the leaf node
### 10/14/2024
### Wen

rm(list=ls())
setwd('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/')
gc()
####install packages
#remotes::install_github("Winnie09/GPTCelltype")
#install.packages("openai")


Sys.setenv(OPENAI_API_KEY = '')

library(GPTCelltype)
library(openai)
library(readxl)
library(Seurat)
library(org.Hs.eg.db)

source('R_code/Utility_func.R')

pbmc_df<-readRDS('data/pbmc_df.RDS')
pbmc_dis_mtx_na<-readRDS('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/pbmc_distance_mtx_na.rda')

# Define the new cell type
#new_cell_type <- "Unknown"

# Define the constant distance (e.g., 117)
#constant_distance <- 117.3096 #max(pbmc_dis_mtx)*1.5

# Add a new row and column with the constant distance value
# The new row and column will have the same distance (117) with all other cell types
#pbmc_dis_mtx_na <- rbind(
#  cbind(pbmc_dis_mtx, rep(constant_distance, nrow(pbmc_dis_mtx))),
#  rep(constant_distance, ncol(pbmc_dis_mtx) + 1)
#)


# Set the row and column names, including the new cell type
#rownames(pbmc_dis_mtx_na) <- c(rownames(pbmc_dis_mtx), 'Unknown')
#colnames(pbmc_dis_mtx_na) <- c(colnames(pbmc_dis_mtx), 'Unknown')
#pbmc_dis_mtx_na_temp<-pbmc_dis_mtx_na
#colnames(pbmc_dis_mtx_na_temp)[1] <- "Blood Cells"
#rownames(pbmc_dis_mtx_na_temp)[1] <- "Blood Cells"

#saveRDS(pbmc_dis_mtx_na_temp,'~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/pbmc_distance_mtx_na.rda')



#_________Step2: Test the cell type predicted by GPTcelltype________#
#pbmc_dis_mtx_na<-readRDS('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/pbmc_distance_mtx_na.rda')
data <- read_excel('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/data_ref_paper/GPT4_paper_supp.xlsx', sheet = 4,skip=3)
data<-as.data.frame(data)


###Deal with PBMC cell types
PBMC_B<-unlist(strsplit(data[2,]$marker,', '))
PBMC_CD4<-unlist(strsplit(data[5,]$marker,', '))
PBMC_CD8<-unlist(strsplit(data[11,]$marker,', '))
PBMC_mono<-unlist(strsplit(data[19,]$marker,', '))
PBMC_plate<-unlist(strsplit(data[27,]$marker,', '))


###________Test some PBMC Marker: downsampling
set.seed(123)
random_8<-sample(10,8)
random_6<-sample(10,6)
random_4<-sample(10,4)
random_2<-sample(10,2)

B_80<-c(PBMC_B[random_8])
B_60<-c(PBMC_B[random_6])
B_40<-c(PBMC_B[random_4])
B_20<-c(PBMC_B[random_2])

CD4_80<-c(PBMC_CD4[random_8])
CD4_60<-c(PBMC_CD4[random_6])
CD4_40<-c(PBMC_CD4[random_4])
CD4_20<-c(PBMC_CD4[random_2])

CD8_80<-c(PBMC_CD8[random_8])
CD8_60<-c(PBMC_CD8[random_6])
CD8_40<-c(PBMC_CD8[random_4])
CD8_20<-c(PBMC_CD8[random_2])

mono_80<-c(PBMC_mono[random_8])
mono_60<-c(PBMC_mono[random_6])
mono_40<-c(PBMC_mono[random_4])
mono_20<-c(PBMC_mono[random_2])

plate_80<-c(PBMC_plate[random_8])
plate_60<-c(PBMC_plate[random_6])
plate_40<-c(PBMC_plate[random_4])
plate_20<-c(PBMC_plate[random_2])

downsample_PBMC<-list(PBMC_B,B_80,B_60,B_40,B_20,
                      PBMC_CD4,CD4_80,CD4_60,CD4_40,CD4_20,
                      PBMC_CD8,CD4_80,CD4_60,CD4_40,CD4_20,
                      PBMC_mono,mono_80,mono_60,mono_40,mono_20,
                      PBMC_plate,plate_80,plate_60,plate_40,plate_20
)

downs_res_test<-CT_GPTpredict(N=2,downsample_PBMC[c(1,6)],tissueName='PBMC',model='o3')

downs_res<-CT_GPTpredict(N=20,downsample_PBMC,tissueName='PBMC',model='o3')
downs_test_res<-apply(downs_res,1,unique)
downs_table_res<-apply(downs_res,1,table)

assign_downs_res<-Assign_Multiple_Cell_Types(downs_table_res,pbmc_df,model = 'o3') # using tree prompt
downs_node_Freq<-Process_AssignRes(assign_downs_res,downs_table_res) 

downs_Dis_res<-Get_avg_dis(assigned_result_list=assign_downs_res,frequency_table_list=downs_table_res,
                           distance_mtx=pbmc_dis_mtx_na,D_Unknown =117.3096 )
downs_Dis<-unlist(downs_Dis_res)
names(downs_Dis)<-c('PBMC_B',paste0('B_',c(80,60,40,20)),
                    'PBMC_CD4',paste0('CD4_',c(80,60,40,20)),
                    'PBMC_CD8',paste0('CD8_',c(80,60,40,20)),
                    'PBMC_mono',paste0('mono_',c(80,60,40,20)),
                    'PBMC_plate',paste0('plate_',c(80,60,40,20)))

downs_Dis
saveRDS(downs_Dis,'result/pbmc_res_data/GPT-o3/Downs_score_o3.RDS')
saveRDS(downs_res,'result/pbmc_res_data/GPT-o3/Downs_pred_25res_o3.RDS')
saveRDS(assign_downs_res,'result/pbmc_res_data/GPT-o3/Downs_Assign_25res_o3.RDS')

###________Test some PBMC Marker: Mixed with other markers from package ____________###
gene<-select(org.Hs.eg.db,columns = 'SYMBOL',keys=keys(org.Hs.eg.db))

set.seed(123)
random_2gene<-sample(nrow(gene),2)
random_4gene<-sample(nrow(gene),4)
random_6gene<-sample(nrow(gene),6)
random_8gene<-sample(nrow(gene),8)

B_80<-c(PBMC_B[random_8],gene[random_2gene,]$SYMBOL)
B_60<-c(PBMC_B[random_6],gene[random_4gene,]$SYMBOL)
B_40<-c(PBMC_B[random_4],gene[random_6gene,]$SYMBOL)
B_20<-c(PBMC_B[random_2],gene[random_8gene,]$SYMBOL)

CD4_80<-c(PBMC_CD4[random_8],gene[random_2gene,]$SYMBOL)
CD4_60<-c(PBMC_CD4[random_6],gene[random_4gene,]$SYMBOL)
CD4_40<-c(PBMC_CD4[random_4],gene[random_6gene,]$SYMBOL)
CD4_20<-c(PBMC_CD4[random_2],gene[random_8gene,]$SYMBOL)

CD8_80<-c(PBMC_CD8[random_8],gene[random_2gene,]$SYMBOL)
CD8_60<-c(PBMC_CD8[random_6],gene[random_4gene,]$SYMBOL)
CD8_40<-c(PBMC_CD8[random_4],gene[random_6gene,]$SYMBOL)
CD8_20<-c(PBMC_CD8[random_2],gene[random_8gene,]$SYMBOL)

mono_80<-c(PBMC_mono[random_8],gene[random_2gene,]$SYMBOL)
mono_60<-c(PBMC_mono[random_6],gene[random_4gene,]$SYMBOL)
mono_40<-c(PBMC_mono[random_4],gene[random_6gene,]$SYMBOL)
mono_20<-c(PBMC_mono[random_2],gene[random_8gene,]$SYMBOL)

plate_80<-c(PBMC_plate[random_8],gene[random_2gene,]$SYMBOL)
plate_60<-c(PBMC_plate[random_6],gene[random_4gene,]$SYMBOL)
plate_40<-c(PBMC_plate[random_4],gene[random_6gene,]$SYMBOL)
plate_20<-c(PBMC_plate[random_2],gene[random_8gene,]$SYMBOL)

Addran_PBMC<-list(PBMC_B,B_80,B_60,B_40,B_20,
                  PBMC_CD4,CD4_80,CD4_60,CD4_40,CD4_20,
                  PBMC_CD8,CD4_80,CD4_60,CD4_40,CD4_20,
                  PBMC_mono,mono_80,mono_60,mono_40,mono_20,
                  PBMC_plate,plate_80,plate_60,plate_40,plate_20
)

Ran_res<-CT_GPTpredict(N=20,Addran_PBMC,tissueName='PBMC',model='o3')
Ran_test_res<-apply(Ran_res,1,unique)
Ran_table_res<-apply(Ran_res,1,table)

assign_Ran_res<-Assign_Multiple_Cell_Types(Ran_table_res,pbmc_df,model='o3') # using tree prompt
Ran_node_Freq<-Process_AssignRes(assign_Ran_res,Ran_table_res) 

Ran_Dis_res<-Get_avg_dis(assigned_result_list=assign_Ran_res,frequency_table_list=Ran_table_res,
                         distance_mtx=pbmc_dis_mtx_na,D_Unknown =117.3096)
Ran_Dis<-unlist(Ran_Dis_res)
names(Ran_Dis)<-c('PBMC_B',paste0('addRan_B_',c(80,60,40,20)),
                  'PBMC_CD4',paste0('addRan_CD4_',c(80,60,40,20)),
                  'PBMC_CD8',paste0('addRan_CD8_',c(80,60,40,20)),
                  'PBMC_mono',paste0('addRan_mono_',c(80,60,40,20)),
                  'PBMC_plate',paste0('addRan_plate_',c(80,60,40,20)))
Ran_Dis
saveRDS(Ran_Dis,'result/pbmc_res_data/GPT-o3/Ran_score_o3.RDS')
saveRDS(Ran_res,'result/pbmc_res_data/GPT-o3/Ran_pred_25res_o3.RDS')
saveRDS(assign_Ran_res,'result/pbmc_res_data/GPT-o3/Ran_Assign_25res_o3.RDS')




###________Test some PBMC Marker: Mixed with other markers from PBMC tissue ____________###
PBMC_GDT<-unlist(strsplit(data[29,]$marker,', '))
PBMC_NK<-unlist(strsplit(data[23,]$marker,', '))
PBMC_dendt<-unlist(strsplit(data[16,]$marker,', '))
PBMC_plasma<-unlist(strsplit(data[4,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

B_80<-c(PBMC_B[random_8],PBMC_NK[random_2gene])
B_60<-c(PBMC_B[random_6],PBMC_NK[random_4gene])
B_40<-c(PBMC_B[random_4],PBMC_NK[random_6gene])
B_20<-c(PBMC_B[random_2],PBMC_NK[random_8gene])

CD4_80<-c(PBMC_CD4[random_8],PBMC_NK[random_2gene])
CD4_60<-c(PBMC_CD4[random_6],PBMC_NK[random_4gene])
CD4_40<-c(PBMC_CD4[random_4],PBMC_NK[random_6gene])
CD4_20<-c(PBMC_CD4[random_2],PBMC_NK[random_8gene])

CD8_80<-c(PBMC_CD8[random_8],PBMC_NK[random_2gene])
CD8_60<-c(PBMC_CD8[random_6],PBMC_NK[random_4gene])
CD8_40<-c(PBMC_CD8[random_4],PBMC_NK[random_6gene])
CD8_20<-c(PBMC_CD8[random_2],PBMC_NK[random_8gene])

mono_80<-c(PBMC_mono[random_8],PBMC_NK[random_2gene])
mono_60<-c(PBMC_mono[random_6],PBMC_NK[random_4gene])
mono_40<-c(PBMC_mono[random_4],PBMC_NK[random_6gene])
mono_20<-c(PBMC_mono[random_2],PBMC_NK[random_8gene])

plate_80<-c(PBMC_plate[random_8],PBMC_NK[random_2gene])
plate_60<-c(PBMC_plate[random_6],PBMC_NK[random_4gene])
plate_40<-c(PBMC_plate[random_4],PBMC_NK[random_6gene])
plate_20<-c(PBMC_plate[random_2],PBMC_NK[random_8gene])

WithinTissue4_PBMC<-list(PBMC_B,B_80,B_60,B_40,B_20,
                         PBMC_CD4,CD4_80,CD4_60,CD4_40,CD4_20,
                         PBMC_CD8,CD4_80,CD4_60,CD4_40,CD4_20,
                         PBMC_mono,mono_80,mono_60,mono_40,mono_20,
                         PBMC_plate,plate_80,plate_60,plate_40,plate_20
)

WithinTissue4_res<-CT_GPTpredict(N=20,WithinTissue4_PBMC,tissueName='PBMC',model = 'o3')
WithinTissue4_test_res<-apply(WithinTissue4_res,1,unique)
WithinTissue4_table_res<-apply(WithinTissue4_res,1,table)

assign_WithinTissue4_res<-Assign_Multiple_Cell_Types(WithinTissue4_table_res,pbmc_df,model = 'o3') # using tree prompt
WithinTissue4_node_Freq<-Process_AssignRes(assign_WithinTissue4_res,WithinTissue4_table_res) 


WithinTissue4_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithinTissue4_res,frequency_table_list=WithinTissue4_table_res,
                                   distance_mtx=pbmc_dis_mtx_na,D_Unknown = 117.3096)
WithinTissue4_Dis<-unlist(WithinTissue4_Dis_res)
names(WithinTissue4_Dis)<-c('PBMC_B',paste0('addNK_B_',c(80,60,40,20)),
                            'PBMC_CD4',paste0('addNK_CD4_',c(80,60,40,20)),
                            'PBMC_CD8',paste0('addNK_CD8_',c(80,60,40,20)),
                            'PBMC_mono',paste0('addNK_mono_',c(80,60,40,20)),
                            'PBMC_plate',paste0('addNK_plate_',c(80,60,40,20)))

WithinTissue4_Dis

saveRDS(WithinTissue4_Dis,'result/pbmc_res_data/GPT-o3/Within_NK_score_o3.RDS')
saveRDS(WithinTissue4_res,'result/pbmc_res_data/GPT-o3/Within_NK_pred_res_o3.RDS')
saveRDS(assign_WithinTissue4_res,'result/pbmc_res_data/GPT-o3/Within_NK_Assign_25res_o3.RDS')


###________Test some PBMC Marker: Mixed with other markers from Other tissue ____________###
Lung_fib<-unlist(strsplit(data[133,]$marker,', '))
Kidney_endo<-unlist(strsplit(data[169,]$marker,', '))
Adipose_fat<-unlist(strsplit(data[266,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

B_80<-c(PBMC_B[random_8],Adipose_fat[random_2gene])
B_60<-c(PBMC_B[random_6],Adipose_fat[random_4gene])
B_40<-c(PBMC_B[random_4],Adipose_fat[random_6gene])
B_20<-c(PBMC_B[random_2],Adipose_fat[random_8gene])

CD4_80<-c(PBMC_CD4[random_8],Adipose_fat[random_2gene])
CD4_60<-c(PBMC_CD4[random_6],Adipose_fat[random_4gene])
CD4_40<-c(PBMC_CD4[random_4],Adipose_fat[random_6gene])
CD4_20<-c(PBMC_CD4[random_2],Adipose_fat[random_8gene])

CD8_80<-c(PBMC_CD8[random_8],Adipose_fat[random_2gene])
CD8_60<-c(PBMC_CD8[random_6],Adipose_fat[random_4gene])
CD8_40<-c(PBMC_CD8[random_4],Adipose_fat[random_6gene])
CD8_20<-c(PBMC_CD8[random_2],Adipose_fat[random_8gene])

mono_80<-c(PBMC_mono[random_8],Adipose_fat[random_2gene])
mono_60<-c(PBMC_mono[random_6],Adipose_fat[random_4gene])
mono_40<-c(PBMC_mono[random_4],Adipose_fat[random_6gene])
mono_20<-c(PBMC_mono[random_2],Adipose_fat[random_8gene])

plate_80<-c(PBMC_plate[random_8],Adipose_fat[random_2gene])
plate_60<-c(PBMC_plate[random_6],Adipose_fat[random_4gene])
plate_40<-c(PBMC_plate[random_4],Adipose_fat[random_6gene])
plate_20<-c(PBMC_plate[random_2],Adipose_fat[random_8gene])

WithoutTissue3_PBMC<-list(PBMC_B,B_80,B_60,B_40,B_20,
                          PBMC_CD4,CD4_80,CD4_60,CD4_40,CD4_20,
                          PBMC_CD8,CD4_80,CD4_60,CD4_40,CD4_20,
                          PBMC_mono,mono_80,mono_60,mono_40,mono_20,
                          PBMC_plate,plate_80,plate_60,plate_40,plate_20
)

WithoutTissue3_res<-CT_GPTpredict(N=20,WithoutTissue3_PBMC,model = 'o3')
WithoutTissue3_test_res<-apply(WithoutTissue3_res,1,unique)
WithoutTissue3_table_res<-apply(WithoutTissue3_res,1,table)

assign_WithoutTissue3_res<-Assign_Multiple_Cell_Types(WithoutTissue3_table_res,pbmc_df,model = 'o3') # using tree prompt
WithoutTissue3_node_Freq<-Process_AssignRes(assign_WithoutTissue3_res,WithoutTissue3_table_res) 

WithoutTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithoutTissue3_res,frequency_table_list=WithoutTissue3_table_res,
                                    distance_mtx=pbmc_dis_mtx_na,D_Unknown = 117.3096)
WithoutTissue3_Dis<-unlist(WithoutTissue3_Dis_res)
names(WithoutTissue3_Dis)<-c('PBMC_B',paste0('addfat_B_',c(80,60,40,20)),
                             'PBMC_CD4',paste0('addfat_CD4_',c(80,60,40,20)),
                             'PBMC_CD8',paste0('addfat_CD8_',c(80,60,40,20)),
                             'PBMC_mono',paste0('addfat_mono_',c(80,60,40,20)),
                             'PBMC_plate',paste0('addfat_plate_',c(80,60,40,20)))

WithoutTissue3_Dis

saveRDS(WithoutTissue3_Dis,'result/pbmc_res_data/GPT-o3/outside_fat_25res_o3.RDS')
saveRDS(WithoutTissue3_res,'result/pbmc_res_data/GPT-o3/outside_fat_pred_res_o3.RDS')
saveRDS(assign_WithoutTissue3_res,'result/pbmc_res_data/GPT-o3/outside_fat_assign_res_o3.RDS')

