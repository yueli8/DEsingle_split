
rds<-readRDS('NKT.rds')
#稀疏矩阵转小再组合----
ca1 <- subset(rds, subset = celltype == "ca1")
ca1_counts<-as.matrix(ca1@assays$RNA@counts)

ca2 <- subset(rds, subset = celltype == "ca2")
ca2_counts <- as.matrix(ca2@assays$RNA@counts)

ca3 <- subset(rds, subset = celltype == "ca3")
ca3_counts <- as.matrix(ca3@assays$RNA@counts)

ad1 <- subset(rds, subset = celltype == "ad1")
ad1_counts <- as.matrix(ad1@assays$RNA@counts)

ad2 <- subset(rds, subset = celltype == "ad2")
ad2_counts <- as.matrix(ad2@assays$RNA@counts)

ad3 <- subset(rds, subset = celltype == "ad3")
ad3_counts <- as.matrix(ad3@assays$RNA@counts)

nt1 <- subset(rds, subset = celltype == "nt1")
nt1_counts <- as.matrix(nt1@assays$RNA@counts)

nt2 <- subset(rds, subset = celltype == "nt2")
nt2_counts <- as.matrix(nt2@assays$RNA@counts)

nt3 <- subset(rds, subset = celltype == "nt3")
nt3_counts <- as.matrix(nt3@assays$RNA@counts)

counts <- cbind(
  ca1_counts,
  ca2_counts,
  ca3_counts,
  ad1_counts,
  ad2_counts,
  ad3_counts,
  nt1_counts,
  nt2_counts,
  nt3_counts
)
#counts<-as.matrix(rds@assays$RNA@counts)
NKT_ca_nt <- subset(rds, 
                    subset = celltype == "ca1" | celltype == "ca2" | celltype == "ca3" | celltype == "nt1" | celltype == "nt2" | celltype == "nt3")

NKT_ca_ad <- subset(rds, 
                    subset = celltype == "ca1" | celltype == "ca2" | celltype == "ca3" | celltype == "ad1" | celltype == "ad2" | celltype == "ad3")

NKT_ad_nt <- subset(rds, 
                    subset = celltype == "ad1" | celltype == "ad2" | celltype == "ad3" | celltype == "nt1" | celltype == "nt2" | celltype == "nt3")
counts_ca_nt <- cbind(ca1_counts,ca2_counts,ca3_counts,nt1_counts,nt2_counts,nt3_counts)
counts_ca_ad <- cbind(ca1_counts,ca2_counts, ca3_counts, ad1_counts, ad2_counts, ad3_counts)
counts_ad_nt <- cbind(ad1_counts, ad2_counts, ad3_counts, nt1_counts, nt2_counts, nt3_counts)

saveRDS(NKT_ca_nt, file = "NKT_ca_nt.rds")
saveRDS(NKT_ca_ad, file = "NKT_ca_ad.rds")
saveRDS(NKT_ad_nt, file = "NKT_ad_nt.rds")
saveRDS(counts_ca_nt, file = "counts_ca_nt.rds")
saveRDS(counts_ca_ad, file = "counts_ca_ad.rds")
saveRDS(counts_ad_nt, file = "counts_ad_nt.rds")