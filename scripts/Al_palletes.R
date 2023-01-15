Al_palettes <- list(
  'main'  = c("#305099","#779F28","#88275E","#CB921E","#4D1893","#226D4D","#914324","#CABA1A"),
  'blues' =c("#8CA1D0","#516CAB","#305099","#17367E","#09225E"),
  'greens' =c("#BBE46B","#95BE43","#779F28","#54770D","#2E4500"),
  'pinks' =c("#D279AB","#AB4780","#88275E","#5D0A39","#33001D"),
  'dark-yellows' =c("#FFC552","#EEB239","#CB921E","#9C6F14","#5B3D01"),
  'purples' =c("#7E5BAC","#63379D","#4D1893","#3B0B7A","#2E0663"),
  'turquoises' =c("#4DAD84","#368D67","#226D4D","#134E35","#09301F"),
  'browns' =c("#F9AB8C","#C4704F","#914324","#61220A","#330F00"),
  'yellows' =c("#FFF058","#E9D937","#CABA1A","#A3960E","#786D00"),
  'main2'=c("#8CA1D0","#BBE46B","#D279AB","#FFC552","#4DAD84","#7E5BAC","#F9AB8C","#FFF058"),
  'main3'=c("#516CAB","#95BE43","#AB4780","#EEB239","#368D67","#C4704F","#63379D","#E9D937"),
  "all"=c("#8CA1D0","#BBE46B","#D279AB","#FFC552","#7E5BAC","#4DAD84","#F9AB8C","#FFF058","#516CAB","#95BE43","#AB4780","#EEB239","#63379D","#368D67","#C4704F","#E9D937","#305099","#779F28","#88275E","#CB921E","#4D1893","#226D4D","#914324","#CABA1A","#17367E","#54770D","#5D0A39","#9C6F14","#3B0B7A","#134E35","#61220A","#A3960E","#09225E","#2E4500","#33001D","#5B3D01","#2E0663","#09301F","#330F00","#786D00")
  
  
  
  
  )

Al_palettes2 <- list(
  'main'  = c("#B4CC00","#038CDD","#FF8C00","#BE01DF","#1A9E31","#44248D","#CEB522","#CE2E22"),
  'pairs' = c("#B4CC00","#ECFB7F","#038CDD","#78BCE4","#FF8C00","#FFC57F","#BE01DF","#D575E6","1A9E31","#80D58E","#44248D","#957FC6","#CEB522","#FFF099","#CE2E22","#FFA099"),
  'olives' =c("#ECFB7F","#E1F92C","#B4CC00","#7F9000","#353C00"),
  'blue_ocean' =c("#78BCE4","#289DE3","#038CDD","#025687","#012F4B"),
  'oranges' =c("#FFC57F","#FF9B20","#FF8C00","#D47500","#764100"),
  'violet' =c("#D575E6","#C721E4","#BE01DF","#77018B","#42004D"),
  'greens' =c("#80D58E","#42B356","#1A9E31","#07761A","#00480D"),
  'purples' =c("#957FC6","#6146A0","#44248D","#2C116A","#170441"),
  'yellows' =c("#FFF099","#EAD456","#CEB522","#9A8509","#5E5100"),
  'reds' =c("#FFA099","#EA6156","#CE2E22","#9A1309","#5E0700"),
  'all' =c("#ECFB7F","#78BCE4","#FFC57F","#D575E6","#80D58E", "#957FC6","#FFF099","#FFA099",
    "#E1F92C","#289DE3","#FF9B20","#C721E4","#42B356", "#6146A0","#EAD456","#EA6156",
    "#B4CC00","#038CDD","#FF8C00","#BE01DF","#1A9E31", "#44248D","#CEB522","#CE2E22",
    "#7F9000","#025687","#D47500","#77018B","#07761A", "#2C116A","#9A8509","#9A1309",
    "#353C00","#012F4B","#764100","#42004D","#00480D", "#170441","#5E5100","#5E0700"),
  'main1'=c("#ECFB7F","#78BCE4","#FFC57F","#D575E6","#80D58E", "#957FC6","#FFF099","#FFA099"),
  'main2'=c("#E1F92C","#289DE3","#FF9B20","#C721E4","#42B356", "#6146A0","#EAD456","#EA6156"),
  'main4'=c("#7F9000","#025687","#D47500","#77018B","#07761A", "#2C116A","#9A8509","#9A1309"),
  'main5'=c("#353C00","#012F4B","#764100","#42004D","#00480D", "#170441","#5E5100","#5E0700")
    )


# 
# A=data.frame(A=c(Al_palettes2$olives,Al_palettes2$blue_ocean,Al_palettes2$oranges,Al_palettes2$violet,
#                  Al_palettes2$greens,Al_palettes2$purples,Al_palettes2$yellows,Al_palettes2$reds), B=rep(1:5, times=8), C=rep(1:8, each=5))
# A$B=factor(A$B); A$C=factor(A$C); A$A=factor(A$A, levels = A$A)
# 
# A %>% ggplot() + geom_tile(aes(x=C, y=B, fill=A))+
#     scale_fill_manual(values = c(Al_palettes2$olives,Al_palettes2$blue_ocean,Al_palettes2$oranges,Al_palettes2$violet,
#                                  Al_palettes2$greens,Al_palettes2$purples,Al_palettes2$yellows,Al_palettes2$reds))+
#   theme_minimal()+ theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank())



