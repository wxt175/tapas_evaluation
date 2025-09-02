### Functions for assign the cell type predicted by GPTcelltype to the leaf node for Lung
### 01/10/2024
### Wen

rm(list=ls())
setwd('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/')
gc()



Sys.setenv(OPENAI_API_KEY = 'sk-proj-xIWYBV3Yy7iXcQXr1SaHT3BlbkFJkGdD00exLCkCY2FkvacE')

library(GPTCelltype)
library(openai)
library(readxl)
library(Seurat)
library(org.Hs.eg.db)

source('R_code/Utility_func.R')

lung_df<-readRDS('data/lung_df.RDS')
lung_dis_mtx_na<-readRDS('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/lung_distance_mtx_na.RDS')

# Define the new cell type
#new_cell_type <- "Unknown"

# Define the constant distance (e.g., 117)
#constant_distance <- max(lung_dis_mtx)*1.5
#constant_distance #90.4405
# Add a new row and column with the constant distance value
# The new row and column will have the same distance (117) with all other cell types
#lung_dis_mtx_na <- rbind(
#  cbind(lung_dis_mtx, rep(constant_distance, nrow(lung_dis_mtx))),
  #  rep(constant_distance, ncol(lung_dis_mtx) + 1)
#)


# Set the row and column names, including the new cell type
#rownames(lung_dis_mtx_na) <- c(rownames(lung_dis_mtx), 'Unknown')
#colnames(lung_dis_mtx_na) <- c(colnames(lung_dis_mtx), 'Unknown')
#saveRDS(lung_dis_mtx_na,'~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/lung_distance_mtx_na.rda')

#_______Step1: Function to assign the cell type_______#
# Example usage
#Assign_ct_Tree("Lower Respiratory Tract Epithelial Cell", lung_df)



#_________Step2: Test the cell type predicted by GPTcelltype________#
data <- read_excel('data/data_ref_paper/GPT4_paper_supp.xlsx', sheet = 4,skip=3)
data<-as.data.frame(data)


###Deal with lung cell types
lung_dend<-unlist(strsplit(data[132,]$marker,', '))
lung_epi<-unlist(strsplit(data[138,]$marker,', '))
lung_fib<-unlist(strsplit(data[133,]$marker,', '))
lung_endo<-unlist(strsplit(data[154,]$marker,', '))
lung_smus<-unlist(strsplit(data[152,]$marker,', '))



marker_lung<-list(lung_dend, lung_epi, lung_fib, lung_endo, lung_smus)

###Using GPTcelltype to predict 20 times
res<-CT_GPTpredict(N=20,marker_lung,tissueName='lung')
test_res<-apply(res,1,unique)
table_res<-apply(res,1,table)

###Test the results
assign_res<-Assign_Multiple_Cell_Types(table_res,lung_df) # using tree prompt
node_Freq<-Process_AssignRes(assign_res,table_res) 

###
#node_Freq<-Process_AssignRes(assign_res,table_res) 
#dis_mtx<-Get_matrix(assign_res,lung_dis_mtx_na,idx=1:5)
Dis_res<-Get_avg_dis(assigned_result_list=assign_res,frequency_table_list=table_res,distance_mtx=lung_dis_mtx_na,D_Unknown = 90.4405)
Dis_res
saveRDS(Dis_res,'result/lung_res_data/Original_20res.RDS')
saveRDS(assign_res,'result/lung_res_data/Original_assign_res.RDS')
saveRDS(node_Freq,'result/lung_res_data/Original_node_Freq.RDS')


###________Test some lung Marker: downsampling
set.seed(123)
random_8<-sample(10,8)
random_6<-sample(10,6)
random_4<-sample(10,4)
random_2<-sample(10,2)

dend_80<-c(lung_dend[random_8])
dend_60<-c(lung_dend[random_6])
dend_40<-c(lung_dend[random_4])
dend_20<-c(lung_dend[random_2])

epi_80<-c(lung_epi[random_8])
epi_60<-c(lung_epi[random_6])
epi_40<-c(lung_epi[random_4])
epi_20<-c(lung_epi[random_2])

fib_80<-c(lung_fib[random_8])
fib_60<-c(lung_fib[random_6])
fib_40<-c(lung_fib[random_4])
fib_20<-c(lung_fib[random_2])

endo_80<-c(lung_endo[random_8])
endo_60<-c(lung_endo[random_6])
endo_40<-c(lung_endo[random_4])
endo_20<-c(lung_endo[random_2])

smus_80<-c(lung_smus[random_8])
smus_60<-c(lung_smus[random_6])
smus_40<-c(lung_smus[random_4])
smus_20<-c(lung_smus[random_2])

downsample_lung<-list(lung_dend,dend_80,dend_60,dend_40,dend_20,
                        lung_epi,epi_80,epi_60,epi_40,epi_20,
                        lung_fib,fib_80,fib_60,fib_40,fib_20,
                        lung_endo,endo_80,endo_60,endo_40,endo_20,
                        lung_smus,smus_80,smus_60,smus_40,smus_20
)

downs_res<-CT_GPTpredict(N=20,downsample_lung,tissueName='lung')
downs_table_res<-apply(downs_res,1,table)

assign_downs_res<-Assign_Multiple_Cell_Types(downs_table_res,lung_df) # using tree prompt
downs_node_Freq<-Process_AssignRes(assign_downs_res,downs_table_res) 

downs_Dis_res<-Get_avg_dis(assigned_result_list=assign_downs_res,frequency_table_list=downs_table_res,
                           distance_mtx=lung_dis_mtx_na,D_Unknown = 90.4405)
downs_Dis<-unlist(downs_Dis_res)
names(downs_Dis)<-c('lung_dend',paste0('dend_',c(80,60,40,20)),
                    'lung_epi',paste0('epi_',c(80,60,40,20)),
                    'lung_fib',paste0('fib_',c(80,60,40,20)),
                    'lung_endo',paste0('endo_',c(80,60,40,20)),
                    'lung_smus',paste0('smus_',c(80,60,40,20)))
saveRDS(downs_Dis,'result/lung_res_data/Downsampling_25res.RDS')
saveRDS(assign_downs_res,'result/lung_res_data/Downsampling_assign_res.RDS')
saveRDS(downs_node_Freq,'result/lung_res_data/Downsampling_node_Freq.RDS')
saveRDS(downs_res,'result/lung_res_data/Downsampling_pred_res.RDS')


###________Test some lung Marker: Mixed with other markers from package ____________###
gene<-select(org.Hs.eg.db,columns = 'SYMBOL',keys=keys(org.Hs.eg.db))

set.seed(123456)
random_2gene<-sample(nrow(gene),2)
random_4gene<-sample(nrow(gene),4)
random_6gene<-sample(nrow(gene),6)
random_8gene<-sample(nrow(gene),8)

dend_80<-c(lung_dend[random_8],gene[random_2gene,]$SYMBOL)
dend_60<-c(lung_dend[random_6],gene[random_4gene,]$SYMBOL)
dend_40<-c(lung_dend[random_4],gene[random_6gene,]$SYMBOL)
dend_20<-c(lung_dend[random_2],gene[random_8gene,]$SYMBOL)

epi_80<-c(lung_epi[random_8],gene[random_2gene,]$SYMBOL)
epi_60<-c(lung_epi[random_6],gene[random_4gene,]$SYMBOL)
epi_40<-c(lung_epi[random_4],gene[random_6gene,]$SYMBOL)
epi_20<-c(lung_epi[random_2],gene[random_8gene,]$SYMBOL)

fib_80<-c(lung_fib[random_8],gene[random_2gene,]$SYMBOL)
fib_60<-c(lung_fib[random_6],gene[random_4gene,]$SYMBOL)
fib_40<-c(lung_fib[random_4],gene[random_6gene,]$SYMBOL)
fib_20<-c(lung_fib[random_2],gene[random_8gene,]$SYMBOL)

endo_80<-c(lung_endo[random_8],gene[random_2gene,]$SYMBOL)
endo_60<-c(lung_endo[random_6],gene[random_4gene,]$SYMBOL)
endo_40<-c(lung_endo[random_4],gene[random_6gene,]$SYMBOL)
endo_20<-c(lung_endo[random_2],gene[random_8gene,]$SYMBOL)

smus_80<-c(lung_smus[random_8],gene[random_2gene,]$SYMBOL)
smus_60<-c(lung_smus[random_6],gene[random_4gene,]$SYMBOL)
smus_40<-c(lung_smus[random_4],gene[random_6gene,]$SYMBOL)
smus_20<-c(lung_smus[random_2],gene[random_8gene,]$SYMBOL)

Addran_lung<-list(lung_dend,dend_80,dend_60,dend_40,dend_20,
                    lung_epi,epi_80,epi_60,epi_40,epi_20,
                    lung_fib,fib_80,fib_60,fib_40,fib_20,
                    lung_endo,endo_80,endo_60,endo_40,endo_20,
                    lung_smus,smus_80,smus_60,smus_40,smus_20
)

Ran_res<-CT_GPTpredict(N=20,Addran_lung,tissueName='lung')
Ran_test_res<-apply(Ran_res,1,unique)
Ran_table_res<-apply(Ran_res,1,table)

assign_Ran_res<-Assign_Multiple_Cell_Types(Ran_table_res,lung_df) # using tree prompt
Ran_node_Freq<-Process_AssignRes(assign_Ran_res,Ran_table_res) 

Ran_Dis_res<-Get_avg_dis(assigned_result_list=assign_Ran_res,frequency_table_list=Ran_table_res,
                         distance_mtx=lung_dis_mtx_na,D_Unknown = 90.4405)
Ran_Dis<-unlist(Ran_Dis_res)
names(Ran_Dis)<-c('lung_dend',paste0('addRan_dend_',c(80,60,40,20)),
                  'lung_epi',paste0('addRan_epi_',c(80,60,40,20)),
                  'lung_fib',paste0('addRan_fib_',c(80,60,40,20)),
                  'lung_endo',paste0('addRan_endo_',c(80,60,40,20)),
                  'lung_smus',paste0('addRan_smus_',c(80,60,40,20)))
Ran_Dis
saveRDS(Ran_Dis,'result/lung_res_data/Random_25res.RDS')
saveRDS(Ran_node_Freq,'result/lung_res_data/Random_node_Freq.RDS')
saveRDS(assign_Ran_res,'result/lung_res_data/Random_assign_res.RDS')
saveRDS(Ran_res,'result/lung_res_data/Random_pred_res.RDS')



###________Test some lung Marker: Mixed with other markers from lung tissue ____________###
lung_myofib<-unlist(strsplit(data[144,]$marker,', '))
lung_mast<-unlist(strsplit(data[139,]$marker,', '))
lung_meso<-unlist(strsplit(data[157,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

dend_80<-c(lung_dend[random_8],lung_meso[random_2gene])
dend_60<-c(lung_dend[random_6],lung_meso[random_4gene])
dend_40<-c(lung_dend[random_4],lung_meso[random_6gene])
dend_20<-c(lung_dend[random_2],lung_meso[random_8gene])

epi_80<-c(lung_epi[random_8],lung_meso[random_2gene])
epi_60<-c(lung_epi[random_6],lung_meso[random_4gene])
epi_40<-c(lung_epi[random_4],lung_meso[random_6gene])
epi_20<-c(lung_epi[random_2],lung_meso[random_8gene])

fib_80<-c(lung_fib[random_8],lung_meso[random_2gene])
fib_60<-c(lung_fib[random_6],lung_meso[random_4gene])
fib_40<-c(lung_fib[random_4],lung_meso[random_6gene])
fib_20<-c(lung_fib[random_2],lung_meso[random_8gene])

endo_80<-c(lung_endo[random_8],lung_meso[random_2gene])
endo_60<-c(lung_endo[random_6],lung_meso[random_4gene])
endo_40<-c(lung_endo[random_4],lung_meso[random_6gene])
endo_20<-c(lung_endo[random_2],lung_meso[random_8gene])

smus_80<-c(lung_smus[random_8],lung_meso[random_2gene])
smus_60<-c(lung_smus[random_6],lung_meso[random_4gene])
smus_40<-c(lung_smus[random_4],lung_meso[random_6gene])
smus_20<-c(lung_smus[random_2],lung_meso[random_8gene])

WithinTissue3_lung<-list(lung_dend,dend_80,dend_60,dend_40,dend_20,
                           lung_epi,epi_80,epi_60,epi_40,epi_20,
                           lung_fib,fib_80,fib_60,fib_40,fib_20,
                           lung_endo,endo_80,endo_60,endo_40,endo_20,
                           lung_smus,smus_80,smus_60,smus_40,smus_20
)

WithinTissue3_res<-CT_GPTpredict(N=20,WithinTissue3_lung,tissueName='lung')
WithinTissue3_test_res<-apply(WithinTissue3_res,1,unique)
WithinTissue3_table_res<-apply(WithinTissue3_res,1,table)

assign_WithinTissue3_res<-Assign_Multiple_Cell_Types(WithinTissue3_table_res,lung_df) # using tree prompt
WithinTissue3_node_Freq<-Process_AssignRes(assign_WithinTissue3_res,WithinTissue3_table_res) 

WithinTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithinTissue3_res,frequency_table_list=WithinTissue3_table_res,
                                   distance_mtx=lung_dis_mtx_na,D_Unknown = 90.4405)
WithinTissue3_Dis<-unlist(WithinTissue3_Dis_res)
names(WithinTissue3_Dis)<-c('lung_dend',paste0('add_mast_dend_',c(80,60,40,20)),
                            'lung_epi',paste0('add_mast_epi_',c(80,60,40,20)),
                            'lung_fib',paste0('add_mast_fib_',c(80,60,40,20)),
                            'lung_endo',paste0('add_mast_endo_',c(80,60,40,20)),
                            'lung_smus',paste0('add_mast_smus_',c(80,60,40,20)))

WithinTissue3_Dis

saveRDS(WithinTissue3_Dis,'result/lung_res_data/Mix_lungmast_25res.RDS')
saveRDS(WithinTissue3_node_Freq,'result/lung_res_data/Mix_lungmast_node_Freq.RDS')
saveRDS(WithinTissue3_res,'result/lung_res_data/Mix_lungmast_pred_res.RDS')

###________Test some lung Marker: Mixed with other markers from Other tissue ____________###
Adipo_fat<-unlist(strsplit(data[266,]$marker,', '))
Brain_astro<-unlist(strsplit(data[31,]$marker,', '))
Muscle_glial<-unlist(strsplit(data[530,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

dend_80<-c(lung_dend[random_8],Muscle_glial[random_2gene])
dend_60<-c(lung_dend[random_6],Muscle_glial[random_4gene])
dend_40<-c(lung_dend[random_4],Muscle_glial[random_6gene])
dend_20<-c(lung_dend[random_2],Muscle_glial[random_8gene])

epi_80<-c(lung_epi[random_8],Muscle_glial[random_2gene])
epi_60<-c(lung_epi[random_6],Muscle_glial[random_4gene])
epi_40<-c(lung_epi[random_4],Muscle_glial[random_6gene])
epi_20<-c(lung_epi[random_2],Muscle_glial[random_8gene])

fib_80<-c(lung_fib[random_8],Muscle_glial[random_2gene])
fib_60<-c(lung_fib[random_6],Muscle_glial[random_4gene])
fib_40<-c(lung_fib[random_4],Muscle_glial[random_6gene])
fib_20<-c(lung_fib[random_2],Muscle_glial[random_8gene])

endo_80<-c(lung_endo[random_8],Muscle_glial[random_2gene])
endo_60<-c(lung_endo[random_6],Muscle_glial[random_4gene])
endo_40<-c(lung_endo[random_4],Muscle_glial[random_6gene])
endo_20<-c(lung_endo[random_2],Muscle_glial[random_8gene])

smus_80<-c(lung_smus[random_8],Muscle_glial[random_2gene])
smus_60<-c(lung_smus[random_6],Muscle_glial[random_4gene])
smus_40<-c(lung_smus[random_4],Muscle_glial[random_6gene])
smus_20<-c(lung_smus[random_2],Muscle_glial[random_8gene])

WithoutTissue3_lung<-list(lung_dend,dend_80,dend_60,dend_40,dend_20,
                            lung_epi,epi_80,epi_60,epi_40,epi_20,
                            lung_fib,fib_80,fib_60,fib_40,fib_20,
                            lung_endo,endo_80,endo_60,endo_40,endo_20,
                            lung_smus,smus_80,smus_60,smus_40,smus_20
)

WithoutTissue3_res<-CT_GPTpredict(N=20,WithoutTissue3_lung)
WithoutTissue3_test_res<-apply(WithoutTissue3_res,1,unique)
WithoutTissue3_table_res<-apply(WithoutTissue3_res,1,table)

assign_WithoutTissue3_res<-Assign_Multiple_Cell_Types(WithoutTissue3_table_res,lung_df) # using tree prompt
WithoutTissue3_node_Freq<-Process_AssignRes(assign_WithoutTissue3_res,WithoutTissue3_table_res) 

WithoutTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithoutTissue3_res,frequency_table_list=WithoutTissue3_table_res,
                                    distance_mtx=lung_dis_mtx_na,D_Unknown = 90.4405)
WithoutTissue3_Dis<-unlist(WithoutTissue3_Dis_res)
names(WithoutTissue3_Dis)<-c('lung_dend',paste0('addMuscle_glial_dend_',c(80,60,40,20)),
                             'lung_epi',paste0('addMuscle_glial_epi_',c(80,60,40,20)),
                             'lung_fib',paste0('addMuscle_glial_fib_',c(80,60,40,20)),
                             'lung_endo',paste0('addMuscle_glial_endo_',c(80,60,40,20)),
                             'lung_smus',paste0('addMuscle_glial_smus_',c(80,60,40,20)))

WithoutTissue3_Dis

saveRDS(WithoutTissue3_Dis,'result/lung_res_data/Mix_Muscle_glial_25res.RDS')
saveRDS(WithoutTissue3_node_Freq,'result/lung_res_data/Mix_Muscle_glial_node_Freq.RDS')
saveRDS(WithoutTissue3_res,'result/lung_res_data/Mix_Muscle_glial_pred_res.RDS')



