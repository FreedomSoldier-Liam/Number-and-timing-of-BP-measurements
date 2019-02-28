
require(kableExtra)

files <- list.files(path='KableWon',all.files=TRUE, 
                    full.names=TRUE, pattern='.R')

for(f in files) source(f)

real_data=rbind(
  readRDS('Datasets/JHS_anly.RDS'),
  readRDS('Datasets/CARDIA_anly.RDS'))

meta_data=read.csv('Datasets/meta_data.csv',stringsAsFactors=F)

meta = create_meta_object(real_data = real_data,
                          meta_data = meta_data,
                          exposure  = 'study')

tab1=make_table_one(meta=meta,as.kable=TRUE,
                    include.pval = FALSE)

print(tab1)