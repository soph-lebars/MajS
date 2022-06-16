
library(ggplot2)
library(ggsignif)
library(dplyr)

png("comp_prob_majs_enz.png")
newdat= data.frame(
  Group = c("HIF1A= - ","HIF1A= 0 ","HIF1A= + "),
  Name = rep(c("SLC2A1","HK2","HK3","ALDOA","GAPDH","PGK1","ENO1","ENO3","LDHA","PDHB","PDHA1","HK1","PDHA2","PFKL","ENO2"),each=3),
  Value =c(0.779151716590602,1.10183766591446,1.7562773458332
           ,0.759471951491222,0.892220894389505,1.16144915140515
           ,0.958588372837362,1.01022132296783,1.11493815810566
           ,0.81280271283367,0.883242125261548,1.02610036442429
           ,0.796911471416812,0.856965065733679,0.978759817333647
           ,1.05226817228769,1.22397165385555,1.57220364769865
           ,0.71153889935849,0.851325296529999,1.1348262212328
           ,0.960245460436472,0.984309496341361,1.03311379038685
           ,0.836795916050019,0.863949617720173,0.919020065835003
           ,0.833932103259344,0.779482168311164,0.669052203592375
           ,0.968690829257304,0.88590323209899,0.71800162777115
           ,0.954022738261349,0.922709928123812,0.859204388181014
           ,0.978870586445989,0.984748620682376,0.996669867493006
           ,1.05781389894616,1.05317508022382,1.04376708757274
           ,0.847149073753584,0.84179364990855,0.830932309728363),
  Value2=c("(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(-,100)","(-,100)",
           "(0,100)","(-,100)","(-,100)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(-,100)","(-,100)",
           "(0,100)","(+,25)","(+,50)",
           "(0,100)","(+,25)","(+,50)"
  )
  
  
)





head(newdat)



newdat$Name <- factor(newdat$Name,levels = c("SLC2A1","HK2","HK3","ALDOA","GAPDH","PGK1","ENO1","ENO3","LDHA","PDHB","PDHA1","HK1","PDHA2","PFKL","ENO2"))
newdat$Name <- factor(newdat$Name,levels = unique(newdat$Name))

a <- ifelse(unique(newdat$Name)  %in% c("HK1","PDHA2","PFKL","ENO2") , "orange", "purple")

G=factor(newdat$Group, levels = c("HIF1A= - ","HIF1A= 0 ","HIF1A= + "))



Value3=as.numeric(factor(newdat$Value2, levels = c("(-,100)","(0,100)","(+,25)","(+,50)")))*0.5


Value3[Value3==0.5]=0.75 #Improve proportion on y axis MajS 

Value3[Value3==1.5]=1.125 

Value3[Value3==2]=1.25 

Nam=as.numeric(newdat$Name)



p <- ggplot(newdat, aes(x=Name,y=Value,fill= G) )+ 
  geom_bar(width=0.80,position="dodge",  colour="black",  stat="identity")+ theme(legend.title=element_blank(), axis.text.y = element_text(size=12), axis.text.x = element_text(size=12,angle = 45, hjust = 1.2, colour = a), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line.y = element_line(colour = "black")) + labs(x="Enzyme") + geom_point(aes(x=Name,y=Value3, color=G),size=3,position =position_dodge(width = 0.9), stat = "identity")+
  geom_line(aes(x=Name,y=Value3),group=Nam, position =position_dodge(width = 0.9), stat = "identity", linetype = "dashed")+
  scale_y_continuous(
    
    name = "FC value predicted by probregnet",

    # Add a second axis and specify its features
    sec.axis = sec_axis(~.,breaks=c(0.75,1,1.125,1.25),labels=c("(-,100)","(0,100)","(+,25)","(+,50)"),name="MajS prediction")
  ) + scale_fill_manual(values=c("lightpink","lightblue","lightgreen")) + scale_colour_manual(values=c("red","blue","darkgreen"))+ theme(legend.position = c(0.9, 0.85),text = element_text(size=15, face= "bold"), legend.text = element_text(size=10),legend.background = element_rect(fill="white",size=0.5, linetype="solid",  colour ="black"))
  
  
  
print(p)

dev.off()

