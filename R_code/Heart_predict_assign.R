### Functions for assign the cell type predicted by GPTcelltype to the leaf node for Heart
### 12/23/2024
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

heart_df<-readRDS('data/heart_df.RDS')
heart_dis_mtx_na<-readRDS('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/Heart_distance_mtx_na.RDS')

# Define the new cell type
#new_cell_type <- "Unknown"

# Define the constant distance (e.g., 117)
#constant_distance <- max(heart_dis_mtx)*1.5
#constant_distance<-79.81907
# Add a new row and column with the constant distance value
# The new row and column will have the same distance (117) with all other cell types
# heart_dis_mtx_na <- rbind(
#   cbind(heart_dis_mtx, rep(constant_distance, nrow(heart_dis_mtx))),
#   rep(constant_distance, ncol(heart_dis_mtx) + 1)
# )


# Set the row and column names, including the new cell type
# rownames(heart_dis_mtx_na) <- c(rownames(heart_dis_mtx), 'Unknown')
# colnames(heart_dis_mtx_na) <- c(colnames(heart_dis_mtx), 'Unknown')
# saveRDS(heart_dis_mtx_na,'~/Desktop/Tapas/data/heart_distance_mtx.rda')
# saveRDS(heart_dis_mtx_na,'~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/heart_distance_mtx_na.rda')
#_______Step1: Function to assign the cell type_______#
# Example usage
# Assign_ct_Tree("Cardiomyocyte", heart_df)



#_________Step2: Test the cell type predicted by GPTcelltype________#
data <- read_excel('data/data_ref_paper/GPT4_paper_supp.xlsx', sheet = 4,skip=3)
data<-as.data.frame(data)


###Deal with heart cell types
heart_fat<-unlist(strsplit(data[477,]$marker,', '))
heart_b<-unlist(strsplit(data[479,]$marker,', '))
heart_muscle<-unlist(strsplit(data[484,]$marker,', '))
heart_endo<-unlist(strsplit(data[487,]$marker,', '))
heart_fib<-unlist(strsplit(data[488,]$marker,', '))



marker_heart<-list(heart_fat, heart_b, heart_muscle, heart_endo, heart_fib)

###Using GPTcelltype to predict 5 times
res<-CT_GPTpredict(N=5,marker_heart,tissueName='Heart')
test_res<-apply(res,1,unique)
table_res<-apply(res,1,table)

###Test the results
assign_res<-Assign_Multiple_Cell_Types(table_res,heart_df) # using tree prompt
node_Freq<-Process_AssignRes(assign_res,table_res) 

###
#node_Freq<-Process_AssignRes(assign_res,table_res) 
#dis_mtx<-Get_matrix(assign_res,heart_dis_mtx_na,idx=1:5)
Dis_res<-Get_avg_dis(assigned_result_list=assign_res,frequency_table_list=table_res,distance_mtx=heart_dis_mtx_na,D_Unknown = 79.81907)
Dis_res
saveRDS(Dis_res,'result/heart_res_data/Original_25res.RDS')
saveRDS(assign_res,'result/heart_res_data/Original_assign_res.RDS')
saveRDS(node_Freq,'result/heart_res_data/Original_node_Freq.RDS')


###________Test some heart Marker: downsampling
set.seed(123)
random_8<-sample(10,8)
random_6<-sample(10,6)
random_4<-sample(10,4)
random_2<-sample(10,2)

fat_80<-c(heart_fat[random_8])
fat_60<-c(heart_fat[random_6])
fat_40<-c(heart_fat[random_4])
fat_20<-c(heart_fat[random_2])

b_80<-c(heart_b[random_8])
b_60<-c(heart_b[random_6])
b_40<-c(heart_b[random_4])
b_20<-c(heart_b[random_2])

muscle_80<-c(heart_muscle[random_8])
muscle_60<-c(heart_muscle[random_6])
muscle_40<-c(heart_muscle[random_4])
muscle_20<-c(heart_muscle[random_2])

endo_80<-c(heart_endo[random_8])
endo_60<-c(heart_endo[random_6])
endo_40<-c(heart_endo[random_4])
endo_20<-c(heart_endo[random_2])

fib_80<-c(heart_fib[random_8])
fib_60<-c(heart_fib[random_6])
fib_40<-c(heart_fib[random_4])
fib_20<-c(heart_fib[random_2])

downsample_heart<-list(heart_fat,fat_80,fat_60,fat_40,fat_20,
                      heart_b,b_80,b_60,b_40,b_20,
                      heart_muscle,muscle_80,muscle_60,muscle_40,muscle_20,
                      heart_endo,endo_80,endo_60,endo_40,endo_20,
                      heart_fib,fib_80,fib_60,fib_40,fib_20
)

downs_res<-CT_GPTpredict(N=20,downsample_heart,tissueName='Heart')
downs_table_res<-apply(downs_res,1,table)

assign_downs_res<-Assign_Multiple_Cell_Types(downs_table_res,heart_df) # using tree prompt
downs_node_Freq<-Process_AssignRes(assign_downs_res,downs_table_res) 

downs_Dis_res<-Get_avg_dis(assigned_result_list=assign_downs_res,frequency_table_list=downs_table_res,
                           distance_mtx=heart_dis_mtx_na,D_Unknown = 79.81907)
downs_Dis<-unlist(downs_Dis_res)
names(downs_Dis)<-c('heart_fat',paste0('fat_',c(80,60,40,20)),
                    'heart_b',paste0('b_',c(80,60,40,20)),
                    'heart_muscle',paste0('muscle_',c(80,60,40,20)),
                    'heart_endo',paste0('endo_',c(80,60,40,20)),
                    'heart_fib',paste0('fib_',c(80,60,40,20)))
saveRDS(downs_Dis,'result/heart_res_data/Downsampling_25res.RDS')
saveRDS(assign_downs_res,'result/heart_res_data/Downsampling_assign_res.RDS')
saveRDS(downs_node_Freq,'result/heart_res_data/Downsampling_node_Freq.RDS')
saveRDS(downs_res,'result/heart_res_data/Downsampling_pred_res.RDS')


###________Test some heart Marker: Mixed with other markers from package ____________###
gene<-select(org.Hs.eg.db,columns = 'SYMBOL',keys=keys(org.Hs.eg.db))

set.seed(1234)
random_2gene<-sample(nrow(gene),2)
random_4gene<-sample(nrow(gene),4)
random_6gene<-sample(nrow(gene),6)
random_8gene<-sample(nrow(gene),8)

fat_80<-c(heart_fat[random_8],gene[random_2gene,]$SYMBOL)
fat_60<-c(heart_fat[random_6],gene[random_4gene,]$SYMBOL)
fat_40<-c(heart_fat[random_4],gene[random_6gene,]$SYMBOL)
fat_20<-c(heart_fat[random_2],gene[random_8gene,]$SYMBOL)

b_80<-c(heart_b[random_8],gene[random_2gene,]$SYMBOL)
b_60<-c(heart_b[random_6],gene[random_4gene,]$SYMBOL)
b_40<-c(heart_b[random_4],gene[random_6gene,]$SYMBOL)
b_20<-c(heart_b[random_2],gene[random_8gene,]$SYMBOL)

muscle_80<-c(heart_muscle[random_8],gene[random_2gene,]$SYMBOL)
muscle_60<-c(heart_muscle[random_6],gene[random_4gene,]$SYMBOL)
muscle_40<-c(heart_muscle[random_4],gene[random_6gene,]$SYMBOL)
muscle_20<-c(heart_muscle[random_2],gene[random_8gene,]$SYMBOL)

endo_80<-c(heart_endo[random_8],gene[random_2gene,]$SYMBOL)
endo_60<-c(heart_endo[random_6],gene[random_4gene,]$SYMBOL)
endo_40<-c(heart_endo[random_4],gene[random_6gene,]$SYMBOL)
endo_20<-c(heart_endo[random_2],gene[random_8gene,]$SYMBOL)

fib_80<-c(heart_fib[random_8],gene[random_2gene,]$SYMBOL)
fib_60<-c(heart_fib[random_6],gene[random_4gene,]$SYMBOL)
fib_40<-c(heart_fib[random_4],gene[random_6gene,]$SYMBOL)
fib_20<-c(heart_fib[random_2],gene[random_8gene,]$SYMBOL)

Addran_heart<-list(heart_fat,fat_80,fat_60,fat_40,fat_20,
                  heart_b,b_80,b_60,b_40,b_20,
                  heart_muscle,b_80,b_60,b_40,b_20,
                  heart_endo,endo_80,endo_60,endo_40,endo_20,
                  heart_fib,fib_80,fib_60,fib_40,fib_20
)

Ran_res<-CT_GPTpredict(N=20,Addran_heart,tissueName='Heart')
Ran_test_res<-apply(Ran_res,1,unique)
Ran_table_res<-apply(Ran_res,1,table)

assign_Ran_res<-Assign_Multiple_Cell_Types(Ran_table_res,heart_df) # using tree prompt
Ran_node_Freq<-Process_AssignRes(assign_Ran_res,Ran_table_res) 

Ran_Dis_res<-Get_avg_dis(assigned_result_list=assign_Ran_res,frequency_table_list=Ran_table_res,
                         distance_mtx=heart_dis_mtx_na,D_Unknown = 79.81907)
Ran_Dis<-unlist(Ran_Dis_res)
names(Ran_Dis)<-c('heart_fat',paste0('addRan_fat_',c(80,60,40,20)),
                  'heart_b',paste0('addRan_b_',c(80,60,40,20)),
                  'heart_muscle',paste0('addRan_muscle_',c(80,60,40,20)),
                  'heart_endo',paste0('addRan_endo_',c(80,60,40,20)),
                  'heart_fib',paste0('addRan_fib_',c(80,60,40,20)))
Ran_Dis
saveRDS(Ran_Dis,'result/heart_res_data/Random_25res.RDS')
saveRDS(Ran_node_Freq,'result/heart_res_data/Random_node_Freq.RDS')
saveRDS(assign_Ran_res,'result/heart_res_data/Random_assign_res.RDS')
saveRDS(Ran_res,'result/heart_res_data/Random_pred_res.RDS')



###________Test some heart Marker: Mixed with other markers from heart tissue ____________###
heart_macroph<-unlist(strsplit(data[291,]$marker,', '))
heart_neuron<-unlist(strsplit(data[295,]$marker,', '))
heart_pericyte<-unlist(strsplit(data[489,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

fat_80<-c(heart_fat[random_8],heart_pericyte[random_2gene])
fat_60<-c(heart_fat[random_6],heart_pericyte[random_4gene])
fat_40<-c(heart_fat[random_4],heart_pericyte[random_6gene])
fat_20<-c(heart_fat[random_2],heart_pericyte[random_8gene])

b_80<-c(heart_b[random_8],heart_pericyte[random_2gene])
b_60<-c(heart_b[random_6],heart_pericyte[random_4gene])
b_40<-c(heart_b[random_4],heart_pericyte[random_6gene])
b_20<-c(heart_b[random_2],heart_pericyte[random_8gene])

muscle_80<-c(heart_muscle[random_8],heart_pericyte[random_2gene])
muscle_60<-c(heart_muscle[random_6],heart_pericyte[random_4gene])
muscle_40<-c(heart_muscle[random_4],heart_pericyte[random_6gene])
muscle_20<-c(heart_muscle[random_2],heart_pericyte[random_8gene])

endo_80<-c(heart_endo[random_8],heart_pericyte[random_2gene])
endo_60<-c(heart_endo[random_6],heart_pericyte[random_4gene])
endo_40<-c(heart_endo[random_4],heart_pericyte[random_6gene])
endo_20<-c(heart_endo[random_2],heart_pericyte[random_8gene])

fib_80<-c(heart_fib[random_8],heart_pericyte[random_2gene])
fib_60<-c(heart_fib[random_6],heart_pericyte[random_4gene])
fib_40<-c(heart_fib[random_4],heart_pericyte[random_6gene])
fib_20<-c(heart_fib[random_2],heart_pericyte[random_8gene])

WithinTissue3_heart<-list(heart_fat,fat_80,fat_60,fat_40,fat_20,
                         heart_b,b_80,b_60,b_40,b_20,
                         heart_muscle,muscle_80,muscle_60,muscle_40,muscle_20,
                         heart_endo,endo_80,endo_60,endo_40,endo_20,
                         heart_fib,fib_80,fib_60,fib_40,fib_20
)

WithinTissue3_res<-CT_GPTpredict(N=20,WithinTissue3_heart,tissueName='Heart')
WithinTissue3_test_res<-apply(WithinTissue3_res,1,unique)
WithinTissue3_table_res<-apply(WithinTissue3_res,1,table)

assign_WithinTissue3_res<-Assign_Multiple_Cell_Types(WithinTissue3_table_res,heart_df) # using tree prompt
WithinTissue3_node_Freq<-Process_AssignRes(assign_WithinTissue3_res,WithinTissue3_table_res) 

WithinTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithinTissue3_res,frequency_table_list=WithinTissue3_table_res,
                                   distance_mtx=heart_dis_mtx_na,D_Unknown = 79.81907)
WithinTissue3_Dis<-unlist(WithinTissue3_Dis_res)
names(WithinTissue3_Dis)<-c('heart_fat',paste0('addpericyte_fat_',c(80,60,40,20)),
                            'heart_b',paste0('addpericyte_b_',c(80,60,40,20)),
                            'heart_muscle',paste0('addpericyte_muscle_',c(80,60,40,20)),
                            'heart_endo',paste0('addpericyte_endo_',c(80,60,40,20)),
                            'heart_fib',paste0('addpericyte_fib_',c(80,60,40,20)))

WithinTissue3_Dis

saveRDS(WithinTissue3_Dis,'result/heart_res_data/Mix_Heartpericyte_25res.RDS')
saveRDS(WithinTissue3_node_Freq,'result/heart_res_data/Mix_Heartpericyte_node_Freq.RDS')
saveRDS(WithinTissue3_res,'result/heart_res_data/Mix_Heartpericyte_pred_res.RDS')


###________Test some heart Marker: Mixed with other markers from Other tissue ____________###
Brain_oligo<-unlist(strsplit(data[43,]$marker,', '))
pan_epi<-unlist(strsplit(data[60,]$marker,', '))
Lung_nk<-unlist(strsplit(data[136,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

fat_80<-c(heart_fat[random_8],Lung_nk[random_2gene])
fat_60<-c(heart_fat[random_6],Lung_nk[random_4gene])
fat_40<-c(heart_fat[random_4],Lung_nk[random_6gene])
fat_20<-c(heart_fat[random_2],Lung_nk[random_8gene])

b_80<-c(heart_b[random_8],Lung_nk[random_2gene])
b_60<-c(heart_b[random_6],Lung_nk[random_4gene])
b_40<-c(heart_b[random_4],Lung_nk[random_6gene])
b_20<-c(heart_b[random_2],Lung_nk[random_8gene])

muscle_80<-c(heart_muscle[random_8],Lung_nk[random_2gene])
muscle_60<-c(heart_muscle[random_6],Lung_nk[random_4gene])
muscle_40<-c(heart_muscle[random_4],Lung_nk[random_6gene])
muscle_20<-c(heart_muscle[random_2],Lung_nk[random_8gene])

endo_80<-c(heart_endo[random_8],Lung_nk[random_2gene])
endo_60<-c(heart_endo[random_6],Lung_nk[random_4gene])
endo_40<-c(heart_endo[random_4],Lung_nk[random_6gene])
endo_20<-c(heart_endo[random_2],Lung_nk[random_8gene])

fib_80<-c(heart_fib[random_8],Lung_nk[random_2gene])
fib_60<-c(heart_fib[random_6],Lung_nk[random_4gene])
fib_40<-c(heart_fib[random_4],Lung_nk[random_6gene])
fib_20<-c(heart_fib[random_2],Lung_nk[random_8gene])

WithoutTissue3_heart<-list(heart_fat,fat_80,fat_60,fat_40,fat_20,
                          heart_b,b_80,b_60,b_40,b_20,
                          heart_muscle,muscle_80,muscle_60,muscle_40,muscle_20,
                          heart_endo,endo_80,endo_60,endo_40,endo_20,
                          heart_fib,fib_80,fib_60,fib_40,fib_20
)

WithoutTissue3_res<-CT_GPTpredict(N=20,WithoutTissue3_heart,tissueName='Heart')
WithoutTissue3_test_res<-apply(WithoutTissue3_res,1,unique)
WithoutTissue3_table_res<-apply(WithoutTissue3_res,1,table)

assign_WithoutTissue3_res<-Assign_Multiple_Cell_Types(WithoutTissue3_table_res,heart_df) # using tree prompt
WithoutTissue3_node_Freq<-Process_AssignRes(assign_WithoutTissue3_res,WithoutTissue3_table_res) 

WithoutTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithoutTissue3_res,frequency_table_list=WithoutTissue3_table_res,
                                    distance_mtx=heart_dis_mtx_na,D_Unknown = 79.81907)
WithoutTissue3_Dis<-unlist(WithoutTissue3_Dis_res)
names(WithoutTissue3_Dis)<-c('heart_fat',paste0('addLung_nk_fat_',c(80,60,40,20)),
                             'heart_b',paste0('addLung_nk_b_',c(80,60,40,20)),
                             'heart_muscle',paste0('addLung_nk_muscle_',c(80,60,40,20)),
                             'heart_endo',paste0('addLung_nk_endo_',c(80,60,40,20)),
                             'heart_fib',paste0('addLung_nk_fib_',c(80,60,40,20)))

WithoutTissue3_Dis

saveRDS(WithoutTissue3_Dis,'result/heart_res_data/Mix_Lung_nk_25res.RDS')
saveRDS(WithoutTissue3_node_Freq,'result/heart_res_data/Mix_Lung_nk_node_Freq.RDS')
saveRDS(WithoutTissue3_res,'result/heart_res_data/Mix_Lung_nk_pred_res.RDS')



