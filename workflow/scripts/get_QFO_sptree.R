library(ape)

t <- read.nhx("resources/species_tree_QFO2020.nx")
t <- TreeTools::Preorder(t@phylo)
t$tip.label <- gsub("_.*", "", t$tip.label)
t <- TreeTools::Subtree(t, getMRCA(t, tip = c("HUMAN", "GIAIC")))
t$tip.label <- gsub("CANFA", "CANLF", t$tip.label)
t <- TreeTools::AddTip(t, getMRCA(t, tip = c("DANRE", "ORYLA")), label = "LEPOC")
t <- TreeTools::AddTip(t, getMRCA(t, tip = c("DROME", "BOMMO")), label = "TRICA")
t <- TreeTools::AddTip(t, where = "TRIVA", label = "TRIV3")
t <- TreeTools::AddTip(t, where = "SCHMA", label = "HELRO")
t <- TreeTools::AddTip(t, where = "SORBI", label = "MAIZE")
# Some species are not found in Euka dataset and this causes a bug in A-Pro
# mammals <- c("PONAB", "NOMLE", "MACMU",
#              "CALJA", "OTOGA", "CAVPO",
#              "SPETR", "RABIT", "PIG", 
#              "HORSE", "CANLF", "AILME", "MYOLU", "LOXAF", "SARHA", "ORNAN", "CHICK")
# t <- TreeTools::DropTip(t, mammals)
write.tree(t, "data/sptrees/QFO.nwk")

# p <- ggtree(t) + 
#   geom_tiplab() + 
#   xlim(c(0,30))
# 
# ggsave("test/qfo.pdf", width = 8, height = 12)
