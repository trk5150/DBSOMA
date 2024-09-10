install.packages("VennDiagram")
library(VennDiagram)

list1 <- readLines("C:/Users/tik105/OneDrive - Harvard University/Atollin Paper/Pulldown mass spec/For R studio venn diagram/176A pulldown 1.txt")
list2 <- readLines("C:/Users/tik105/OneDrive - Harvard University/Atollin Paper/Pulldown mass spec/For R studio venn diagram/176a pulldown 2.txt")
list3 <- readLines("C:/Users/tik105/OneDrive - Harvard University/Atollin Paper/Pulldown mass spec/For R studio venn diagram/IG.txt")

list1 <- toupper(list1)
list2 <- toupper(list2)
list3 <- toupper(list3)
library(RColorBrewer)
myCol <- brewer.pal(3, "Set1")


venn.plot <- venn.diagram(
  x = list(
    List1 = list1,
    List2 = list2,
    List3 = list3
  ),
  category.names = c("", "", ""),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = myCol,
  alpha = 0.50,
  cex = 2,
  fontfamily = "serif",
  fontface = "bold",
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.pos = 0,
  cat.dist = 0.05,
  cat.default.pos = "outer"
)

# To plot it in R
grid.draw(venn.plot)
