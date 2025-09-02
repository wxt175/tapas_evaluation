### Functions for assign the cell type predicted by GPTcelltype to the leaf node for kidney
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

kidney_df<-readRDS('data/kidney_df.RDS')
kidney_dis_mtx_na<-readRDS('~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/kidney_distance_mtx_na.RDS')

# Define the new cell type
# new_cell_type <- "Unknown"

# Define the constant distance (e.g., 117)
# constant_distance <- max(kidney_dis_mtx)*1.5
# constant_distance #102.9212
# Add a new row and column with the constant distance value
# The new row and column will have the same distance (117) with all other cell types
# kidney_dis_mtx_na <- rbind(
#   cbind(kidney_dis_mtx, rep(constant_distance, nrow(kidney_dis_mtx))),
#   rep(constant_distance, ncol(kidney_dis_mtx) + 1)
# )


# Set the row and column names, including the new cell type
# rownames(kidney_dis_mtx_na) <- c(rownames(kidney_dis_mtx), 'Unknown')
# colnames(kidney_dis_mtx_na) <- c(colnames(kidney_dis_mtx), 'Unknown')
# saveRDS(kidney_dis_mtx_na,'~/Desktop/RA_case/2024_project/04262024_GPTCelltype/data/kidney_distance_mtx_na.rda')

#_______Step1: Function to assign the cell type_______#
# Example usage
# Assign_ct_Tree("podocyte", kidney_df)



#_________Step2: Test the cell type predicted by GPTcelltype________#
data <- read_excel('data/data_ref_paper/GPT4_paper_supp.xlsx', sheet = 4,skip=3)
data<-as.data.frame(data)


###Deal with kidney cell types
kidney_podo<-unlist(strsplit(data[175,]$marker,','))
kidney_epi<-unlist(strsplit(data[177,]$marker,','))
kidney_thick<-unlist(strsplit(data[179,]$marker,','))
kidney_endo<-unlist(strsplit(data[169,]$marker,','))
kidney_b<-unlist(strsplit(data[777,]$marker,','))



marker_kidney<-list(kidney_podo, kidney_epi, kidney_thick, kidney_endo, kidney_b)

###Using GPTcelltype to predict 20 times
res<-CT_GPTpredict(N=20,marker_kidney,tissueName='Kidney')
test_res<-apply(res,1,unique)
table_res<-apply(res,1,table)

###Test the results
assign_res<-Assign_Multiple_Cell_Types(table_res,kidney_df) # using tree prompt
node_Freq<-Process_AssignRes(assign_res,table_res) 

###
#node_Freq<-Process_AssignRes(assign_res,table_res) 
#dis_mtx<-Get_matrix(assign_res,kidney_dis_mtx_na,idx=1:5)
Dis_res<-Get_avg_dis(assigned_result_list=assign_res,frequency_table_list=table_res,distance_mtx=kidney_dis_mtx_na,D_Unknown = 102.9212)
Dis_res
saveRDS(Dis_res,'result/kidney_res_data/Original_25res.RDS')
saveRDS(assign_res,'result/kidney_res_data/Original_assign_res.RDS')
saveRDS(node_Freq,'result/kidney_res_data/Original_node_Freq.RDS')


###________Test some kidney Marker: downsampling
set.seed(123)
random_8<-sample(10,8)
random_6<-sample(10,6)
random_4<-sample(10,4)
random_2<-sample(10,2)

podo_80<-c(kidney_podo[random_8])
podo_60<-c(kidney_podo[random_6])
podo_40<-c(kidney_podo[random_4])
podo_20<-c(kidney_podo[random_2])

epi_80<-c(kidney_epi[random_8])
epi_60<-c(kidney_epi[random_6])
epi_40<-c(kidney_epi[random_4])
epi_20<-c(kidney_epi[random_2])

thick_80<-c(kidney_thick[random_8])
thick_60<-c(kidney_thick[random_6])
thick_40<-c(kidney_thick[random_4])
thick_20<-c(kidney_thick[random_2])

endo_80<-c(kidney_endo[random_8])
endo_60<-c(kidney_endo[random_6])
endo_40<-c(kidney_endo[random_4])
endo_20<-c(kidney_endo[random_2])

b_80<-c(kidney_b[random_8])
b_60<-c(kidney_b[random_6])
b_40<-c(kidney_b[random_4])
b_20<-c(kidney_b[random_2])

downsample_kidney<-list(kidney_podo,podo_80,podo_60,podo_40,podo_20,
                       kidney_epi,epi_80,epi_60,epi_40,epi_20,
                       kidney_thick,thick_80,thick_60,thick_40,thick_20,
                       kidney_endo,endo_80,endo_60,endo_40,endo_20,
                       kidney_b,b_80,b_60,b_40,b_20
)

downs_res<-CT_GPTpredict(N=20,downsample_kidney,tissueName='kidney')
downs_table_res<-apply(downs_res,1,table)

assign_downs_res<-Assign_Multiple_Cell_Types(downs_table_res,kidney_df) # using tree prompt
downs_node_Freq<-Process_AssignRes(assign_downs_res,downs_table_res) 

downs_Dis_res<-Get_avg_dis(assigned_result_list=assign_downs_res,frequency_table_list=downs_table_res,
                           distance_mtx=kidney_dis_mtx_na,D_Unknown = 102.9212)
downs_Dis<-unlist(downs_Dis_res)
names(downs_Dis)<-c('kidney_podo',paste0('podo_',c(80,60,40,20)),
                    'kidney_epi',paste0('epi_',c(80,60,40,20)),
                    'kidney_thick',paste0('thick_',c(80,60,40,20)),
                    'kidney_endo',paste0('endo_',c(80,60,40,20)),
                    'kidney_b',paste0('b_',c(80,60,40,20)))
saveRDS(downs_Dis,'result/kidney_res_data/Downsampling_25res.RDS')
saveRDS(assign_downs_res,'result/kidney_res_data/Downsampling_assign_res.RDS')
saveRDS(downs_node_Freq,'result/kidney_res_data/Downsampling_node_Freq.RDS')
saveRDS(downs_res,'result/kidney_res_data/Downsampling_pred_res.RDS')


###________Test some kidney Marker: Mixed with other markers from package ____________###
gene<-select(org.Hs.eg.db,columns = 'SYMBOL',keys=keys(org.Hs.eg.db))

set.seed(123456)
random_2gene<-sample(nrow(gene),2)
random_4gene<-sample(nrow(gene),4)
random_6gene<-sample(nrow(gene),6)
random_8gene<-sample(nrow(gene),8)

podo_80<-c(kidney_podo[random_8],gene[random_2gene,]$SYMBOL)
podo_60<-c(kidney_podo[random_6],gene[random_4gene,]$SYMBOL)
podo_40<-c(kidney_podo[random_4],gene[random_6gene,]$SYMBOL)
podo_20<-c(kidney_podo[random_2],gene[random_8gene,]$SYMBOL)

epi_80<-c(kidney_epi[random_8],gene[random_2gene,]$SYMBOL)
epi_60<-c(kidney_epi[random_6],gene[random_4gene,]$SYMBOL)
epi_40<-c(kidney_epi[random_4],gene[random_6gene,]$SYMBOL)
epi_20<-c(kidney_epi[random_2],gene[random_8gene,]$SYMBOL)

thick_80<-c(kidney_thick[random_8],gene[random_2gene,]$SYMBOL)
thick_60<-c(kidney_thick[random_6],gene[random_4gene,]$SYMBOL)
thick_40<-c(kidney_thick[random_4],gene[random_6gene,]$SYMBOL)
thick_20<-c(kidney_thick[random_2],gene[random_8gene,]$SYMBOL)

endo_80<-c(kidney_endo[random_8],gene[random_2gene,]$SYMBOL)
endo_60<-c(kidney_endo[random_6],gene[random_4gene,]$SYMBOL)
endo_40<-c(kidney_endo[random_4],gene[random_6gene,]$SYMBOL)
endo_20<-c(kidney_endo[random_2],gene[random_8gene,]$SYMBOL)

b_80<-c(kidney_b[random_8],gene[random_2gene,]$SYMBOL)
b_60<-c(kidney_b[random_6],gene[random_4gene,]$SYMBOL)
b_40<-c(kidney_b[random_4],gene[random_6gene,]$SYMBOL)
b_20<-c(kidney_b[random_2],gene[random_8gene,]$SYMBOL)

Addran_kidney<-list(kidney_podo,podo_80,podo_60,podo_40,podo_20,
                   kidney_epi,epi_80,epi_60,epi_40,epi_20,
                   kidney_thick,thick_80,thick_60,thick_40,thick_20,
                   kidney_endo,endo_80,endo_60,endo_40,endo_20,
                   kidney_b,b_80,b_60,b_40,b_20
)

Ran_res<-CT_GPTpredict(N=20,Addran_kidney,tissueName='kidney')
Ran_test_res<-apply(Ran_res,1,unique)
Ran_table_res<-apply(Ran_res,1,table)

assign_Ran_res<-Assign_Multiple_Cell_Types(Ran_table_res,kidney_df) # using tree prompt
Ran_node_Freq<-Process_AssignRes(assign_Ran_res,Ran_table_res) 

Ran_Dis_res<-Get_avg_dis(assigned_result_list=assign_Ran_res,frequency_table_list=Ran_table_res,
                         distance_mtx=kidney_dis_mtx_na,D_Unknown = 102.9212)
Ran_Dis<-unlist(Ran_Dis_res)
names(Ran_Dis)<-c('kidney_podo',paste0('addRan_podo_',c(80,60,40,20)),
                  'kidney_epi',paste0('addRan_epi_',c(80,60,40,20)),
                  'kidney_thick',paste0('addRan_thick_',c(80,60,40,20)),
                  'kidney_endo',paste0('addRan_endo_',c(80,60,40,20)),
                  'kidney_b',paste0('addRan_b_',c(80,60,40,20)))
Ran_Dis
saveRDS(Ran_Dis,'result/kidney_res_data/Random_25res.RDS')
saveRDS(Ran_node_Freq,'result/kidney_res_data/Random_node_Freq.RDS')
saveRDS(assign_Ran_res,'result/kidney_res_data/Random_assign_res.RDS')
saveRDS(Ran_res,'result/kidney_res_data/Random_pred_res.RDS')



###________Test some kidney Marker: Mixed with other markers from kidney tissue ____________###
kidney_distal<-unlist(strsplit(data[168,]$marker,','))
kidney_nk<-unlist(strsplit(data[781,]$marker,','))
kidney_inter<-unlist(strsplit(data[172,]$marker,','))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

podo_80<-c(kidney_podo[random_8],kidney_inter[random_2gene])
podo_60<-c(kidney_podo[random_6],kidney_inter[random_4gene])
podo_40<-c(kidney_podo[random_4],kidney_inter[random_6gene])
podo_20<-c(kidney_podo[random_2],kidney_inter[random_8gene])

epi_80<-c(kidney_epi[random_8],kidney_inter[random_2gene])
epi_60<-c(kidney_epi[random_6],kidney_inter[random_4gene])
epi_40<-c(kidney_epi[random_4],kidney_inter[random_6gene])
epi_20<-c(kidney_epi[random_2],kidney_inter[random_8gene])

thick_80<-c(kidney_thick[random_8],kidney_inter[random_2gene])
thick_60<-c(kidney_thick[random_6],kidney_inter[random_4gene])
thick_40<-c(kidney_thick[random_4],kidney_inter[random_6gene])
thick_20<-c(kidney_thick[random_2],kidney_inter[random_8gene])

endo_80<-c(kidney_endo[random_8],kidney_inter[random_2gene])
endo_60<-c(kidney_endo[random_6],kidney_inter[random_4gene])
endo_40<-c(kidney_endo[random_4],kidney_inter[random_6gene])
endo_20<-c(kidney_endo[random_2],kidney_inter[random_8gene])

b_80<-c(kidney_b[random_8],kidney_inter[random_2gene])
b_60<-c(kidney_b[random_6],kidney_inter[random_4gene])
b_40<-c(kidney_b[random_4],kidney_inter[random_6gene])
b_20<-c(kidney_b[random_2],kidney_inter[random_8gene])

WithinTissue3_kidney<-list(kidney_podo,podo_80,podo_60,podo_40,podo_20,
                          kidney_epi,epi_80,epi_60,epi_40,epi_20,
                          kidney_thick,thick_80,thick_60,thick_40,thick_20,
                          kidney_endo,endo_80,endo_60,endo_40,endo_20,
                          kidney_b,b_80,b_60,b_40,b_20
)

WithinTissue3_res<-CT_GPTpredict(N=20,WithinTissue3_kidney,tissueName='kidney')
WithinTissue3_test_res<-apply(WithinTissue3_res,1,unique)
WithinTissue3_table_res<-apply(WithinTissue3_res,1,table)

assign_WithinTissue3_res<-Assign_Multiple_Cell_Types(WithinTissue3_table_res,kidney_df) # using tree prompt
WithinTissue3_node_Freq<-Process_AssignRes(assign_WithinTissue3_res,WithinTissue3_table_res) 

WithinTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithinTissue3_res,frequency_table_list=WithinTissue3_table_res,
                                   distance_mtx=kidney_dis_mtx_na,D_Unknown = 102.9212)
WithinTissue3_Dis<-unlist(WithinTissue3_Dis_res)
names(WithinTissue3_Dis)<-c('kidney_podo',paste0('add_inter_podo_',c(80,60,40,20)),
                            'kidney_epi',paste0('add_inter_epi_',c(80,60,40,20)),
                            'kidney_thick',paste0('add_inter_thick_',c(80,60,40,20)),
                            'kidney_endo',paste0('add_inter_endo_',c(80,60,40,20)),
                            'kidney_b',paste0('add_inter_b_',c(80,60,40,20)))

WithinTissue3_Dis

saveRDS(WithinTissue3_Dis,'result/kidney_res_data/Mix_kidneyinter_25res.RDS')
saveRDS(WithinTissue3_node_Freq,'result/kidney_res_data/Mix_kidneyinter_node_Freq.RDS')
saveRDS(WithinTissue3_res,'result/kidney_res_data/Mix_kidneyinter_pred_res.RDS')


###________Test some kidney Marker: Mixed with other markers from Other tissue ____________###
Adipo_fat<-unlist(strsplit(data[266,]$marker,', '))
Brain_astro<-unlist(strsplit(data[31,]$marker,', '))
Muscle_fib<-unlist(strsplit(data[523,]$marker,', '))

set.seed(12345)
random_2gene<-sample(10,2)
random_4gene<-sample(10,4)
random_6gene<-sample(10,6)
random_8gene<-sample(10,8)

podo_80<-c(kidney_podo[random_8],Adipo_fat[random_2gene])
podo_60<-c(kidney_podo[random_6],Adipo_fat[random_4gene])
podo_40<-c(kidney_podo[random_4],Adipo_fat[random_6gene])
podo_20<-c(kidney_podo[random_2],Adipo_fat[random_8gene])

epi_80<-c(kidney_epi[random_8],Adipo_fat[random_2gene])
epi_60<-c(kidney_epi[random_6],Adipo_fat[random_4gene])
epi_40<-c(kidney_epi[random_4],Adipo_fat[random_6gene])
epi_20<-c(kidney_epi[random_2],Adipo_fat[random_8gene])

thick_80<-c(kidney_thick[random_8],Adipo_fat[random_2gene])
thick_60<-c(kidney_thick[random_6],Adipo_fat[random_4gene])
thick_40<-c(kidney_thick[random_4],Adipo_fat[random_6gene])
thick_20<-c(kidney_thick[random_2],Adipo_fat[random_8gene])

endo_80<-c(kidney_endo[random_8],Adipo_fat[random_2gene])
endo_60<-c(kidney_endo[random_6],Adipo_fat[random_4gene])
endo_40<-c(kidney_endo[random_4],Adipo_fat[random_6gene])
endo_20<-c(kidney_endo[random_2],Adipo_fat[random_8gene])

b_80<-c(kidney_b[random_8],Adipo_fat[random_2gene])
b_60<-c(kidney_b[random_6],Adipo_fat[random_4gene])
b_40<-c(kidney_b[random_4],Adipo_fat[random_6gene])
b_20<-c(kidney_b[random_2],Adipo_fat[random_8gene])

WithoutTissue3_kidney<-list(kidney_podo,podo_80,podo_60,podo_40,podo_20,
                           kidney_epi,epi_80,epi_60,epi_40,epi_20,
                           kidney_thick,thick_80,thick_60,thick_40,thick_20,
                           kidney_endo,endo_80,endo_60,endo_40,endo_20,
                           kidney_b,b_80,b_60,b_40,b_20
)

WithoutTissue3_res<-CT_GPTpredict(N=20,WithoutTissue3_kidney)
WithoutTissue3_test_res<-apply(WithoutTissue3_res,1,unique)
WithoutTissue3_table_res<-apply(WithoutTissue3_res,1,table)

assign_WithoutTissue3_res<-Assign_Multiple_Cell_Types(WithoutTissue3_table_res,kidney_df) # using tree prompt
WithoutTissue3_node_Freq<-Process_AssignRes(assign_WithoutTissue3_res,WithoutTissue3_table_res) 

WithoutTissue3_Dis_res<-Get_avg_dis(assigned_result_list=assign_WithoutTissue3_res,frequency_table_list=WithoutTissue3_table_res,
                                    distance_mtx=kidney_dis_mtx_na,D_Unknown = 102.9212)
WithoutTissue3_Dis<-unlist(WithoutTissue3_Dis_res)
names(WithoutTissue3_Dis)<-c('kidney_podo',paste0('addMuscle_fib_podo_',c(80,60,40,20)),
                             'kidney_epi',paste0('addMuscle_fib_epi_',c(80,60,40,20)),
                             'kidney_thick',paste0('addMuscle_fib_thick_',c(80,60,40,20)),
                             'kidney_endo',paste0('addMuscle_fib_endo_',c(80,60,40,20)),
                             'kidney_b',paste0('addMuscle_fib_b_',c(80,60,40,20)))

WithoutTissue3_Dis

saveRDS(WithoutTissue3_Dis,'result/kidney_res_data/Mix_Muscle_fib_25res.RDS')
saveRDS(WithoutTissue3_node_Freq,'result/kidney_res_data/Mix_Muscle_fib_node_Freq.RDS')
saveRDS(WithoutTissue3_res,'result/kidney_res_data/Mix_Muscle_fib_pred_res.RDS')



