1.admixtools2安装
install.packages("devtools") # if "devtools" is not installed already
devtools::install_github("uqrmaie1/admixtools")
library(admixtools)

2.f2-statistics计算
prefix = 'D:/ProfileAdress/93SBbedFile_FidRenamed/93SB_MasterCore_ChrandFidRenamed'
my_f2_dir = 'D:/ProfileAdress/93SB_f2_dir/'
extract_f2(prefix, my_f2_dir, fst = FALSE, afprod = FALSE, verbose = TRUE)
#从 extract_f2() 生成的目录中读取预计算的 f2 块数据
f2_data <- f2_from_precomp(my_f2_dir) 

pop_names <- dimnames(f2_data)[[1]]
cat("可用的群体名称:\n")
print(pop_names)

3.f3-statistics计算
For example:
pop1 = 'G4SubA'
pop2 = c('G2')
pop3 = c('M2subgr1', 'M2subgr2')
qp3pop(f2_data, pop1, pop2, pop3)
