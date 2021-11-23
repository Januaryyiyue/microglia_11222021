# https://www.biostars.org/p/456294/

library(ggplot2)

# import data
h.GO.data <- read.xlsx("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/SCC/GO_human_patch_only.xlsx", sheetName = 'GO_human_patch_only')

h.GO.data.low <- read.xlsx("/Users/januaryyiyue/Desktop/SchoolWorkLife/year3.5/camh/project/SCC/GO_human_patch_low_only.xlsx", sheetName = 'GO_human_patch_low_only')

# add AUC_low column to h.GO.data
h.GO.data.low <- h.GO.data.low[2:4] %>% as.data.frame()
h.GO.data <- inner_join(h.GO.data, h.GO.data.low, by = 'MainTitle')

# filter for adjusted p value <= 0.05
h.GO.data <- filter(h.GO.data, adj.P.Val <= 0.05)

# get 10 biggest AUC and 10 smallest AUC ----
# order based on AUC
h.GO.data <- arrange(h.GO.data, desc(AUC.x))

# rowname to column to facilitate selection of top and bottom 15 GO terms
h.GO.data <- rownames_to_column(h.GO.data)

# keep top and bottom 15 for high
h.GO.data2 <- h.GO.data[1:15,]
h.GO.data3 <- h.GO.data[657:671,]
h.GO.hi <- bind_rows(h.GO.data2, h.GO.data3)
h.GO.hi$AUC <- h.GO.hi$AUC.x
h.GO.hi <- dplyr::select(h.GO.hi, MainTitle, AUC, adj.P.Val, ID, aspect)
h.GO.hi$group <- 'high'

# keep top and bottom 15 for low
h.GO.lo <- bind_rows(h.GO.data2, h.GO.data3)
h.GO.lo$AUC <- h.GO.lo$AUC.y
h.GO.lo <- dplyr::select(h.GO.lo, MainTitle, AUC, adj.P.Val, ID, aspect)
h.GO.lo$group <- 'low'

# combine hi and low
h.GO.final <- bind_rows(h.GO.hi, h.GO.lo)

# plot ----
ggplot(data = h.GO.final, aes(x = group, y = MainTitle, 
                        color = -log10(adj.P.Val), size = AUC)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("Human Patch-Seq GO Dotplot")
