library(ggmuller)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(scales)
library(gridExtra)
setwd("/home/guillemet/Documents/crispr/work/Martin/scripts")


tailletitre=15

myplots <- list()
counter=1

for (r in c('W1','W2','W3','W4','W5','W6','W7','W8')){
  print(r)
  edges_file=paste0("../steps/muller/csv/edges_",r,".csv")
  edges=read.csv(edges_file,row.names = NULL, colClasses = 'character')
  edges=edges[,c(2,3)]
  edges[edges$Parent=='DGCC','Parent']='NA'
  edges=rbind(edges, c('NA', 'DGCC'))


  pop_file=paste0("../steps/muller/csv/first_apop_",r,".csv")
  pop=read.csv(pop_file,row.names = NULL, colClasses = 'character')
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop=pop[,c(2,3,4)]
  pop=rbind(pop,c('NA', 1e-10, 0))
  pop=rbind(pop,c('NA', 1e-10, 1))
  pop=rbind(pop,c('NA', 1e-10, 2))
  pop=rbind(pop,c('NA', 1e-10, 3))
  pop=rbind(pop,c('NA', 1e-10, 4))
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop[pop$Population!=0,]


  ###colors
  getPalette = colorRampPalette(brewer.pal(20, "Spectral"))
  colors=c(getPalette(20)[seq(3,18)],'#CCCCCC','white')
  coeff0=1
  alpha0=200
  for (i in seq(1,length(colors))){
    parentcol=colors[i]
    colors[i]=rgb(col2rgb(parentcol)['red',1]*coeff0, col2rgb(parentcol)['green',1]*coeff0, col2rgb(parentcol)['blue',1]*coeff0, max=255, alpha = alpha0)
  }

  #orderbim=c('971', '16236', '31065', '7037', '21039', '34608', '29998', '23084', '25461', '24343', '31149', '3233', '30386', '27013', '1209', '31725','DGCC','NA')
  orderbim=c('31725','1209','27013','30386','3233','25461','24343','29998','31149','23084','21039','34608','7037','16236','31065','971','DGCC','NA')
  names(colors)=orderbim




  for (i in unique(orderbim)) {
    diff=setdiff(unique(pop$Identity), orderbim)
    parentcol=colors[i]
    coeff=0.6
    newcol= rgb(col2rgb(parentcol)['red',1]*coeff, col2rgb(parentcol)['green',1]*coeff, col2rgb(parentcol)['blue',1]*coeff, max=255, alpha = 255)
    reps=0
    nams=c()
    for (gen in diff){
      if (grepl("-",gen)){
        if (str_split(gen, "-")[[1]][1]==i){
          reps=reps+1
          nams=c(nams,gen)
        }
      }
    }
    win=replicate(reps, newcol)
    names(win)=nams

    for (nam in names(win)){
      if (str_count(nam, "-")>1){
        win[nam] = rgb(col2rgb(parentcol)['red',1]*coeff*coeff, col2rgb(parentcol)['green',1]*coeff*coeff, col2rgb(parentcol)['blue',1]*coeff*coeff, max=255, alpha = 255)
      }
    }

    colors=c(colors,win)
  }

  ndiff=setdiff(unique(pop$Identity), names(colors))
  getPaletteother = colorRampPalette(brewer.pal(20, "Greys"))
  other=replicate(length(ndiff),"#666666")
  names(other)=ndiff
  colors=c(colors, other)
  colors=colors[sort(names(colors))]

  myorder=c('971','31065','16236','7037','34608','21039','23084','31149','29998','24343','25461','3233','30386','27013','1209','31725','DGCC','NA')
  myorder=c(myorder, setdiff(pop$Identity, myorder))

  pop=pop[order(factor(pop$Identity, levels=myorder)),]


  Muller_df <- get_Muller_df(edges, pop)


  Muller_df$lines=rep(0,length(rownames(Muller_df)))

  Muller_df$lines=rep(0,length(rownames(Muller_df)))
  for (i in seq(2, length(rownames(Muller_df))-2)){
    if(Muller_df[i,'Identity'] == Muller_df[i+1,'Identity']){
      Muller_df[i,'lines']=0
    }
    else if(Muller_df[i,'Identity'] == Muller_df[i-1,'Identity']){
      Muller_df[i,'lines']=2*Muller_df[i,'Population']
    }
    else {
      Muller_df[i,'lines']=Muller_df[i,'Population']
    }
  }
  Muller_df_pop <- add_empty_pop(Muller_df)


  for (i in seq(1, length(rownames(Muller_df_pop))-1)){
    if(is.na(Muller_df_pop[i,'Identity'])){
      Muller_df_pop[i,'lines']=Muller_df_pop[i,'Population']
    }
  }
  id_list <- sort(unique(Muller_df_pop$Identity)) # list of legend entries, omitting NA
  ccount=length(unique(pop$Identity))

  contour=c()
  for (i in pop$Identity){
    nam=names(contour)
    if (str_count(i,'-')==0){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')==1){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')>1){
      contour=c(contour, rgb(255, 255 ,255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
  }

  Muller_df_pop[nrow(Muller_df_pop), "lines"]=0
  #Muller_df_pop=Muller_df_pop[(Muller_df_pop$Population>1e-10),]
  Muller_df_pop$Population=Muller_df_pop$lines

  Muller_df_pop$x=-10
  Muller_df_pop$y=-10


  recombi=c()
  #find recombination
  for (gen in unique(Muller_df_pop$Identity)){
    cpt=0
    for (sp in orderbim){
      if (grepl(sp, gen)){
        cpt=cpt+1
      }
    }
    if (cpt>1){
      recombi=c(recombi,gen)
    }
  }

  Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=0
  Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=0

  x=c()
  y=c()
  i=18
  for (gen in recombi){
    subdf= Muller_df_pop[which(Muller_df_pop$Identity==gen),]
    maxpop=max(subdf$Population)
    maxgen=round(subdf[which(subdf$Population==maxpop),'Generation'][[1]])
    tdf=Muller_df_pop[Muller_df_pop$Generation==maxgen,]
    nsubdf=tdf[seq(1, which(tdf$Identity==gen)[1]),]
    nsubdf1=tdf[seq(1, which(tdf$Identity==gen)[2]),]
    Muller_df_pop[i, "x"]=maxgen
    Muller_df_pop[i, "y"]=0.5 * (sum(nsubdf$lines) + sum(nsubdf1$lines))
    i=i+1
  }

  Muller_df_pop[Muller_df_pop[,"x"]==4, "x"]=3.97


  # subdf_plot=Muller_df_pop[Muller_df_pop$Generation==0,]
  # sum0=sum(subdf_plot$Population)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']+0.5*(10-sum0)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']+0.5*(10-sum0)
  # Muller_df_pop[,'y']=Muller_df_pop[,'y']+0.5*(10-sum0)

  r1=r
  if (grepl('W', r)){r=paste0('B', str_split(r,'W')[[1]][2])}
  if (grepl('R', r)){r=paste0('C', str_split(r,'R')[[1]][2])}

  if (grepl('1', r) || grepl('3', r) || grepl('5', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('2', r) || grepl('4', r) || grepl('6', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('7', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80)) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('8', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    #guides(linetype = FALSE, color = FALSE) +
    ylab("") +
    xlab("") +
    theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
    geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
    #geom_point(aes(x,y), fill="white", color='white') +
    theme(legend.position="none", axis.text=element_text(size=25),axis.title=element_text(size=24), axis.title.x=element_blank(), axis.text.x=element_text(size=7), axis.title.y=element_blank(), axis.text.y=element_text(size=7), plot.title = element_text(size = tailletitre, hjust=0.06, vjust=-5, color=rgb(1.0, 0.4980392156862745, 0.054901960784313725)), plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"),) +
    scale_x_continuous(expand=c(0,0), limits=c(0,4.1)) +
    scale_y_continuous(breaks=c(0,2,4,6,8), labels=math_format(expr = 10^.x, format = force), limits=c(0,9), expand=c(0,0))

  myplots[[counter]]=a
  counter=counter+1

  # outfile=paste0("../steps/muller/plots/ttt_muller_",r)
  # pdf(file = outfile, width = 6.5, height = 5)
  # print(a)
  # dev.off()

  if (!grepl('C',r)){
    for (t in c(0,1,2,3,4)){
      variant=read.csv(paste0("../steps/muller/variants/",r1,"_variant.csv"))
      p<-ggplot(variant[variant$Time==t,], aes(x=Genotype, y=f_variant, fill=Genotype, colour=Genotype)) +
        geom_bar(stat="identity")+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=colors) + labs(title=paste0(r,'_t',t)) + ylim(0,1)+scale_color_manual(values=contour)
      #pdf(paste0("../steps/muller/variants/",r,'_',t,".pdf"),width   = 9, height = 5)
      #print(p)
      #dev.off()
    }
  }
}

pdf(paste0("../steps/muller/first_B.pdf"),width = 6, height = 9)
grid.arrange(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]],  ncol=2, nrow = 4)
dev.off()



library(ggmuller)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(scales)
library(gridExtra)
setwd("Documents/crispr/work/Martin/scripts/")

myplots <- list()
counter=1

for (r in c('R1','R2','R3','R4','R5','R6','R7','R8')){
  print(r)
  edges_file=paste0("../steps/muller/csv/edges_",r,".csv")
  edges=read.csv(edges_file,row.names = NULL, colClasses = 'character')
  edges=edges[,c(2,3)]
  edges[edges$Parent=='DGCC','Parent']='NA'
  edges=rbind(edges, c('NA', 'DGCC'))


  pop_file=paste0("../steps/muller/csv/first_apop_",r,".csv")
  pop=read.csv(pop_file,row.names = NULL, colClasses = 'character')
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop=pop[,c(2,3,4)]
  pop=rbind(pop,c('NA', 1e-10, 0))
  pop=rbind(pop,c('NA', 1e-10, 1))
  pop=rbind(pop,c('NA', 1e-10, 2))
  pop=rbind(pop,c('NA', 1e-10, 3))
  pop=rbind(pop,c('NA', 1e-10, 4))
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop[pop$Population!=0,]


  ###colors
  getPalette = colorRampPalette(brewer.pal(20, "Spectral"))
  colors=c(getPalette(20)[seq(3,18)],'#CCCCCC','white')
  coeff0=1
  alpha0=200
  for (i in seq(1,length(colors))){
    parentcol=colors[i]
    colors[i]=rgb(col2rgb(parentcol)['red',1]*coeff0, col2rgb(parentcol)['green',1]*coeff0, col2rgb(parentcol)['blue',1]*coeff0, max=255, alpha = alpha0)
  }

  #orderbim=c('971', '16236', '31065', '7037', '21039', '34608', '29998', '23084', '25461', '24343', '31149', '3233', '30386', '27013', '1209', '31725','DGCC','NA')
  orderbim=c('31725','1209','27013','30386','3233','25461','24343','29998','31149','23084','21039','34608','7037','16236','31065','971','DGCC','NA')
  names(colors)=orderbim




  for (i in unique(orderbim)) {
    diff=setdiff(unique(pop$Identity), orderbim)
    parentcol=colors[i]
    coeff=0.6
    newcol= rgb(col2rgb(parentcol)['red',1]*coeff, col2rgb(parentcol)['green',1]*coeff, col2rgb(parentcol)['blue',1]*coeff, max=255, alpha = 255)
    reps=0
    nams=c()
    for (gen in diff){
      if (grepl("-",gen)){
        if (str_split(gen, "-")[[1]][1]==i){
          reps=reps+1
          nams=c(nams,gen)
        }
      }
    }
    win=replicate(reps, newcol)
    names(win)=nams

    for (nam in names(win)){
      if (str_count(nam, "-")>1){
        win[nam] = rgb(col2rgb(parentcol)['red',1]*coeff*coeff, col2rgb(parentcol)['green',1]*coeff*coeff, col2rgb(parentcol)['blue',1]*coeff*coeff, max=255, alpha = 255)
      }
    }

    colors=c(colors,win)
  }

  ndiff=setdiff(unique(pop$Identity), names(colors))
  getPaletteother = colorRampPalette(brewer.pal(20, "Greys"))
  other=replicate(length(ndiff),"#666666")
  names(other)=ndiff
  colors=c(colors, other)
  colors=colors[sort(names(colors))]

  myorder=c('971','31065','16236','7037','34608','21039','23084','31149','29998','24343','25461','3233','30386','27013','1209','31725','DGCC','NA')
  myorder=c(myorder, setdiff(pop$Identity, myorder))

  pop=pop[order(factor(pop$Identity, levels=myorder)),]


  Muller_df <- get_Muller_df(edges, pop)


  Muller_df$lines=rep(0,length(rownames(Muller_df)))

  Muller_df$lines=rep(0,length(rownames(Muller_df)))
  for (i in seq(2, length(rownames(Muller_df))-2)){
    if(Muller_df[i,'Identity'] == Muller_df[i+1,'Identity']){
      Muller_df[i,'lines']=0
    }
    else if(Muller_df[i,'Identity'] == Muller_df[i-1,'Identity']){
      Muller_df[i,'lines']=2*Muller_df[i,'Population']
    }
    else {
      Muller_df[i,'lines']=Muller_df[i,'Population']
    }
  }
  Muller_df_pop <- add_empty_pop(Muller_df)


  for (i in seq(1, length(rownames(Muller_df_pop))-1)){
    if(is.na(Muller_df_pop[i,'Identity'])){
      Muller_df_pop[i,'lines']=Muller_df_pop[i,'Population']
    }
  }
  id_list <- sort(unique(Muller_df_pop$Identity)) # list of legend entries, omitting NA
  ccount=length(unique(pop$Identity))

  contour=c()
  for (i in pop$Identity){
    nam=names(contour)
    if (str_count(i,'-')==0){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')==1){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')>1){
      contour=c(contour, rgb(255, 255 ,255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
  }

  Muller_df_pop[nrow(Muller_df_pop), "lines"]=0
  #Muller_df_pop=Muller_df_pop[(Muller_df_pop$Population>1e-10),]
  Muller_df_pop$Population=Muller_df_pop$lines

  Muller_df_pop$x=-10
  Muller_df_pop$y=-10


  recombi=c()
  #find recombination
  for (gen in unique(Muller_df_pop$Identity)){
    cpt=0
    for (sp in orderbim){
      if (grepl(sp, gen)){
        cpt=cpt+1
      }
    }
    if (cpt>1){
      recombi=c(recombi,gen)
    }
  }

  Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=0
  Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=0

  x=c()
  y=c()
  i=18
  for (gen in recombi){
    subdf= Muller_df_pop[which(Muller_df_pop$Identity==gen),]
    maxpop=max(subdf$Population)
    maxgen=round(subdf[which(subdf$Population==maxpop),'Generation'][[1]])
    tdf=Muller_df_pop[Muller_df_pop$Generation==maxgen,]
    nsubdf=tdf[seq(1, which(tdf$Identity==gen)[1]),]
    nsubdf1=tdf[seq(1, which(tdf$Identity==gen)[2]),]
    Muller_df_pop[i, "x"]=maxgen
    Muller_df_pop[i, "y"]=0.5 * (sum(nsubdf$lines) + sum(nsubdf1$lines))
    i=i+1
  }

  Muller_df_pop[Muller_df_pop[,"x"]==4, "x"]=3.97


  # subdf_plot=Muller_df_pop[Muller_df_pop$Generation==0,]
  # sum0=sum(subdf_plot$Population)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']+0.5*(10-sum0)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']+0.5*(10-sum0)
  # Muller_df_pop[,'y']=Muller_df_pop[,'y']+0.5*(10-sum0)

  r1=r
  if (grepl('W', r)){r=paste0('B', str_split(r,'W')[[1]][2])}
  if (grepl('R', r)){r=paste0('C', str_split(r,'R')[[1]][2])}

  if (grepl('1', r) || grepl('3', r) || grepl('5', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('2', r) || grepl('4', r) || grepl('6', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('7', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80)) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('8', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    #guides(linetype = FALSE, color = FALSE) +
    ylab("") +
    xlab("") +
    theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
    geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
    #geom_point(aes(x,y), fill="white", color='white') +
    theme(legend.position="none", axis.text=element_text(size=25),axis.title=element_text(size=24), axis.title.x=element_blank(), axis.text.x=element_text(size=7), axis.title.y=element_blank(), axis.text.y=element_text(size=7), plot.title = element_text(size = tailletitre, hjust=0.06, vjust=-5, color=rgb(0.8392156862745098, 0.15294117647058825, 0.1568627450980392)), plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"),) +
    scale_x_continuous(expand=c(0,0), limits=c(0,4.1)) +
    scale_y_continuous(breaks=c(0,2,4,6,8), labels=math_format(expr = 10^.x, format = force), limits=c(0,9), expand=c(0,0))

  myplots[[counter]]=a
  counter=counter+1

  # outfile=paste0("../steps/muller/plots/ttt_muller_",r)
  # pdf(file = outfile, width = 6.5, height = 5)
  # print(a)
  # dev.off()

  if (!grepl('C',r)){
    for (t in c(0,1,2,3,4)){
      variant=read.csv(paste0("../steps/muller/variants/",r1,"_variant.csv"))
      p<-ggplot(variant[variant$Time==t,], aes(x=Genotype, y=f_variant, fill=Genotype, colour=Genotype)) +
        geom_bar(stat="identity")+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=colors) + labs(title=paste0(r,'_t',t)) + ylim(0,1)+scale_color_manual(values=contour)
      #pdf(paste0("../steps/muller/variants/",r,'_',t,".pdf"),width   = 9, height = 5)
      #print(p)
      #dev.off()
    }
  }
}

#pdf(paste0("../steps/muller/first_C.pdf"),width = 6, height = 9)
grid.arrange(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]], myplots[[8]],  ncol=2, nrow = 4)
#dev.off()





library(ggmuller)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(scales)
library(gridExtra)
setwd("Documents/crispr/work/Martin/scripts/")

myplots <- list()
counter=1

for (r in c('C1','C2','C3','C4','C5','C6','C7')){
  print(r)
  edges_file=paste0("../steps/muller/csv/edges_",r,".csv")
  edges=read.csv(edges_file,row.names = NULL, colClasses = 'character')
  edges=edges[,c(2,3)]
  edges[edges$Parent=='DGCC','Parent']='NA'
  edges=rbind(edges, c('NA', 'DGCC'))


  pop_file=paste0("../steps/muller/csv/apop_",r,".csv")
  pop=read.csv(pop_file,row.names = NULL, colClasses = 'character')
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop=pop[,c(2,3,4)]
  pop=rbind(pop,c('NA', 1e-10, 0))
  pop=rbind(pop,c('NA', 1e-10, 1))
  pop=rbind(pop,c('NA', 1e-10, 2))
  pop=rbind(pop,c('NA', 1e-10, 3))
  pop=rbind(pop,c('NA', 1e-10, 4))
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop[pop$Population!=0,]


  ###colors
  getPalette = colorRampPalette(brewer.pal(20, "Spectral"))
  colors=c(getPalette(20)[seq(3,18)],'#CCCCCC','white')
  coeff0=1
  alpha0=200
  for (i in seq(1,length(colors))){
    parentcol=colors[i]
    colors[i]=rgb(col2rgb(parentcol)['red',1]*coeff0, col2rgb(parentcol)['green',1]*coeff0, col2rgb(parentcol)['blue',1]*coeff0, max=255, alpha = alpha0)
  }

  #orderbim=c('971', '16236', '31065', '7037', '21039', '34608', '29998', '23084', '25461', '24343', '31149', '3233', '30386', '27013', '1209', '31725','DGCC','NA')
  orderbim=c('31725','1209','27013','30386','3233','25461','24343','29998','31149','23084','21039','34608','7037','16236','31065','971','DGCC','NA')
  names(colors)=orderbim




  for (i in unique(orderbim)) {
    diff=setdiff(unique(pop$Identity), orderbim)
    parentcol=colors[i]
    coeff=0.6
    newcol= rgb(col2rgb(parentcol)['red',1]*coeff, col2rgb(parentcol)['green',1]*coeff, col2rgb(parentcol)['blue',1]*coeff, max=255, alpha = 255)
    reps=0
    nams=c()
    for (gen in diff){
      if (grepl("-",gen)){
        if (str_split(gen, "-")[[1]][1]==i){
          reps=reps+1
          nams=c(nams,gen)
        }
      }
    }
    win=replicate(reps, newcol)
    names(win)=nams

    for (nam in names(win)){
      if (str_count(nam, "-")>1){
        win[nam] = rgb(col2rgb(parentcol)['red',1]*coeff*coeff, col2rgb(parentcol)['green',1]*coeff*coeff, col2rgb(parentcol)['blue',1]*coeff*coeff, max=255, alpha = 255)
      }
    }

    colors=c(colors,win)
  }

  ndiff=setdiff(unique(pop$Identity), names(colors))
  getPaletteother = colorRampPalette(brewer.pal(20, "Greys"))
  other=replicate(length(ndiff),"#666666")
  names(other)=ndiff
  colors=c(colors, other)
  colors=colors[sort(names(colors))]

  myorder=c('971','31065','16236','7037','34608','21039','23084','31149','29998','24343','25461','3233','30386','27013','1209','31725','DGCC','NA')
  myorder=c(myorder, setdiff(pop$Identity, myorder))

  pop=pop[order(factor(pop$Identity, levels=myorder)),]


  Muller_df <- get_Muller_df(edges, pop)


  Muller_df$lines=rep(0,length(rownames(Muller_df)))

  Muller_df$lines=rep(0,length(rownames(Muller_df)))
  for (i in seq(2, length(rownames(Muller_df))-2)){
    if(Muller_df[i,'Identity'] == Muller_df[i+1,'Identity']){
      Muller_df[i,'lines']=0
    }
    else if(Muller_df[i,'Identity'] == Muller_df[i-1,'Identity']){
      Muller_df[i,'lines']=2*Muller_df[i,'Population']
    }
    else {
      Muller_df[i,'lines']=Muller_df[i,'Population']
    }
  }
  Muller_df_pop <- add_empty_pop(Muller_df)


  for (i in seq(1, length(rownames(Muller_df_pop))-1)){
    if(is.na(Muller_df_pop[i,'Identity'])){
      Muller_df_pop[i,'lines']=Muller_df_pop[i,'Population']
    }
  }
  id_list <- sort(unique(Muller_df_pop$Identity)) # list of legend entries, omitting NA
  ccount=length(unique(pop$Identity))

  contour=c()
  for (i in pop$Identity){
    nam=names(contour)
    if (str_count(i,'-')==0){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')==1){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')>1){
      contour=c(contour, rgb(255, 255 ,255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
  }

  Muller_df_pop[nrow(Muller_df_pop), "lines"]=0
  #Muller_df_pop=Muller_df_pop[(Muller_df_pop$Population>1e-10),]
  Muller_df_pop$Population=Muller_df_pop$lines

  Muller_df_pop$x=-10
  Muller_df_pop$y=-10


  recombi=c()
  #find recombination
  for (gen in unique(Muller_df_pop$Identity)){
    cpt=0
    for (sp in orderbim){
      if (grepl(sp, gen)){
        cpt=cpt+1
      }
    }
    if (cpt>1){
      recombi=c(recombi,gen)
    }
  }

  Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=0
  Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=0

  x=c()
  y=c()
  i=18
  for (gen in recombi){
    subdf= Muller_df_pop[which(Muller_df_pop$Identity==gen),]
    maxpop=max(subdf$Population)
    maxgen=round(subdf[which(subdf$Population==maxpop),'Generation'][[1]])
    tdf=Muller_df_pop[Muller_df_pop$Generation==maxgen,]
    nsubdf=tdf[seq(1, which(tdf$Identity==gen)[1]),]
    nsubdf1=tdf[seq(1, which(tdf$Identity==gen)[2]),]
    Muller_df_pop[i, "x"]=maxgen
    Muller_df_pop[i, "y"]=0.5 * (sum(nsubdf$lines) + sum(nsubdf1$lines))
    i=i+1
  }

  Muller_df_pop[Muller_df_pop[,"x"]==4, "x"]=3.97


  # subdf_plot=Muller_df_pop[Muller_df_pop$Generation==0,]
  # sum0=sum(subdf_plot$Population)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']+0.5*(10-sum0)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']+0.5*(10-sum0)
  # Muller_df_pop[,'y']=Muller_df_pop[,'y']+0.5*(10-sum0)

  r1=r
  if (grepl('C', r)){r=paste0('A', str_split(r,'C')[[1]][2])}
  if (grepl('R', r)){r=paste0('P', str_split(r,'R')[[1]][2])}

  if (grepl('1', r) || grepl('3', r) || grepl('5', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('2', r) || grepl('4', r) || grepl('6', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('7', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80)) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  if (grepl('8', r)){
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
      #guides(linetype = FALSE, color = FALSE) +
      ylab("") +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) +
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }

  a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) +
    #guides(linetype = FALSE, color = FALSE) +
    ylab("Bacteria Population (/ml)") +
    xlab("Time (days)") +
    theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
    geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+
    #geom_point(aes(x,y), fill="white", color='white') +
    theme(legend.position="none", axis.text=element_text(size=25),axis.title=element_text(size=24), axis.title.x=element_blank(), axis.text.x=element_text(size=7), axis.title.y=element_blank(), axis.text.y=element_text(size=7), plot.title = element_text(size = tailletitre, hjust=0.06, vjust=-5, color=rgb(0.12156862745098039, 0.4666666666666667, 0.7058823529411765)), plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"),) +
    scale_x_continuous(expand=c(0,0), limits=c(0,4.1)) +
    scale_y_continuous(breaks=c(0,2,4,6,8), labels=math_format(expr = 10^.x, format = force), limits=c(0,9), expand=c(0,0))

  myplots[[counter]]=a
  counter=counter+1
}

p1 <- ggplot() + theme_void()
#pdf(paste0("../steps/muller/first_A.pdf"),width = 6, height = 9)
grid.arrange(p1,myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[5]], myplots[[6]], myplots[[7]],  ncol=2, nrow = 4)
#dev.off()






############## legend

library(ggmuller)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(scales)
library(gridExtra)
setwd("Documents/crispr/work/Martin/scripts/")

myplots <- list()
counter=1

for (r in c('C1','C2','C3','C4','C5','C6','C7')){
  print(r)
  edges_file=paste0("../steps/muller/csv/edges_",r,".csv")
  edges=read.csv(edges_file,row.names = NULL, colClasses = 'character')
  edges=edges[,c(2,3)]
  edges[edges$Parent=='DGCC','Parent']='NA'
  edges=rbind(edges, c('NA', 'DGCC'))
  
  
  pop_file=paste0("../steps/muller/csv/apop_",r,".csv")
  pop=read.csv(pop_file,row.names = NULL, colClasses = 'character')
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop=pop[,c(2,3,4)]
  pop=rbind(pop,c('NA', 1e-10, 0))
  pop=rbind(pop,c('NA', 1e-10, 1))
  pop=rbind(pop,c('NA', 1e-10, 2))
  pop=rbind(pop,c('NA', 1e-10, 3))
  pop=rbind(pop,c('NA', 1e-10, 4))
  pop$Population=as.numeric(pop$Population)
  pop$Generation=as.numeric(pop$Generation)
  pop[pop$Population!=0,]
  
  
  ###colors
  getPalette = colorRampPalette(brewer.pal(20, "Spectral"))
  colors=c(getPalette(20)[seq(3,18)],'#CCCCCC','white')
  coeff0=1
  alpha0=200
  for (i in seq(1,length(colors))){
    parentcol=colors[i]
    colors[i]=rgb(col2rgb(parentcol)['red',1]*coeff0, col2rgb(parentcol)['green',1]*coeff0, col2rgb(parentcol)['blue',1]*coeff0, max=255, alpha = alpha0)
  }
  
  #orderbim=c('971', '16236', '31065', '7037', '21039', '34608', '29998', '23084', '25461', '24343', '31149', '3233', '30386', '27013', '1209', '31725','DGCC','NA')
  orderbim=c('31725', '1209', '27013', '30386', '3233', '25461', '24343', '29998', '31149', '23084', '21039', '34608', '7037', '16236', '31065', '971','DGCC','NA')
  names(colors)=orderbim
  
  
  
  
  for (i in unique(orderbim)) {
    diff=setdiff(unique(pop$Identity), orderbim)
    parentcol=colors[i]
    coeff=0.6
    newcol= rgb(col2rgb(parentcol)['red',1]*coeff, col2rgb(parentcol)['green',1]*coeff, col2rgb(parentcol)['blue',1]*coeff, max=255, alpha = 255)
    reps=0
    nams=c()
    for (gen in diff){
      if (grepl("-",gen)){
        if (str_split(gen, "-")[[1]][1]==i){
          reps=reps+1
          nams=c(nams,gen)
        }
      }
    }
    win=replicate(reps, newcol)
    names(win)=nams
    
    for (nam in names(win)){
      if (str_count(nam, "-")>1){
        win[nam] = rgb(col2rgb(parentcol)['red',1]*coeff*coeff, col2rgb(parentcol)['green',1]*coeff*coeff, col2rgb(parentcol)['blue',1]*coeff*coeff, max=255, alpha = 255)
      }
    }
    
    colors=c(colors,win)
  }
  
  ndiff=setdiff(unique(pop$Identity), names(colors))
  getPaletteother = colorRampPalette(brewer.pal(20, "Greys"))
  other=replicate(length(ndiff),"#666666")
  names(other)=ndiff
  colors=c(colors, other)
  colors=colors[sort(names(colors))]
  
  myorder=c('971','31065','16236','7037','34608','21039','23084','31149','29998','24343','25461','3233','30386','27013','1209','31725','DGCC','NA')
  myorder=c(myorder, setdiff(pop$Identity, myorder))
  
  pop=pop[order(factor(pop$Identity, levels=myorder)),]
  
  
  Muller_df <- get_Muller_df(edges, pop)
  
  
  Muller_df$lines=rep(0,length(rownames(Muller_df)))
  
  Muller_df$lines=rep(0,length(rownames(Muller_df)))
  for (i in seq(2, length(rownames(Muller_df))-2)){
    if(Muller_df[i,'Identity'] == Muller_df[i+1,'Identity']){
      Muller_df[i,'lines']=0
    }
    else if(Muller_df[i,'Identity'] == Muller_df[i-1,'Identity']){
      Muller_df[i,'lines']=2*Muller_df[i,'Population']
    }
    else {
      Muller_df[i,'lines']=Muller_df[i,'Population']
    }
  }
  Muller_df_pop <- add_empty_pop(Muller_df)
  
  
  for (i in seq(1, length(rownames(Muller_df_pop))-1)){
    if(is.na(Muller_df_pop[i,'Identity'])){
      Muller_df_pop[i,'lines']=Muller_df_pop[i,'Population']
    }
  }
  id_list <- sort(unique(Muller_df_pop$Identity)) # list of legend entries, omitting NA
  levels(id_list)=myorder
  ccount=length(unique(pop$Identity))
  
  contour=c()
  for (i in pop$Identity){
    nam=names(contour)
    if (str_count(i,'-')==0){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')==1){
      contour=c(contour, rgb(255, 255, 255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
    if (str_count(i,'-')>1){
      contour=c(contour, rgb(255, 255 ,255, max=255, alpha = 255))
      names(contour)=c(nam, i)
    }
  } 
  
  Muller_df_pop[nrow(Muller_df_pop), "lines"]=0
  #Muller_df_pop=Muller_df_pop[(Muller_df_pop$Population>1e-10),]
  Muller_df_pop$Population=Muller_df_pop$lines
  
  Muller_df_pop$x=-10
  Muller_df_pop$y=-10
  
  
  recombi=c()
  #find recombination
  for (gen in unique(Muller_df_pop$Identity)){
    cpt=0
    for (sp in orderbim){
      if (grepl(sp, gen)){
        cpt=cpt+1
      }
    }
    if (cpt>1){
      recombi=c(recombi,gen)
    }
  }
  
  Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=0
  Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=0
  
  x=c()
  y=c()
  i=18
  for (gen in recombi){
    subdf= Muller_df_pop[which(Muller_df_pop$Identity==gen),]
    maxpop=max(subdf$Population)
    maxgen=round(subdf[which(subdf$Population==maxpop),'Generation'][[1]])
    tdf=Muller_df_pop[Muller_df_pop$Generation==maxgen,]
    nsubdf=tdf[seq(1, which(tdf$Identity==gen)[1]),]
    nsubdf1=tdf[seq(1, which(tdf$Identity==gen)[2]),]
    Muller_df_pop[i, "x"]=maxgen
    Muller_df_pop[i, "y"]=0.5 * (sum(nsubdf$lines) + sum(nsubdf1$lines))
    i=i+1
  }
  
  Muller_df_pop[Muller_df_pop[,"x"]==4, "x"]=3.97
  
  
  # subdf_plot=Muller_df_pop[Muller_df_pop$Generation==0,]
  # sum0=sum(subdf_plot$Population)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']=Muller_df_pop[is.na(Muller_df_pop$Identity),'Population']+0.5*(10-sum0)
  # Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']=Muller_df_pop[is.na(Muller_df_pop$Identity),'lines']+0.5*(10-sum0)
  # Muller_df_pop[,'y']=Muller_df_pop[,'y']+0.5*(10-sum0)
  colors=colors[myorder]
  r1=r
  if (grepl('W', r)){r=paste0('M', str_split(r,'W')[[1]][2])}
  if (grepl('R', r)){r=paste0('P', str_split(r,'R')[[1]][2])}
  levels(Muller_df_pop$Identity)=myorder
  if (grepl('1', r) || grepl('3', r) || grepl('5', r)){ 
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      #guides(linetype = FALSE, color = FALSE) + 
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) + 
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+ 
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }
  
  if (grepl('2', r) || grepl('4', r) || grepl('6', r)){ 
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      #guides(linetype = FALSE, color = FALSE) + 
      ylab("") +
      xlab("") +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) + 
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+ 
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank()) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }
  
  if (grepl('7', r)){ 
    a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
      #guides(linetype = FALSE, color = FALSE) + 
      ylab(paste0("Population in ",r,"  (/ml)")) +
      xlab(paste0("Time (days)")) +
      theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
      xlim(0,4) + 
      geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+ 
      #geom_point(aes(x,y), fill="white", color='white') +
      theme(legend.position="none", axis.text=element_text(size=30),axis.title=element_text(size=30), plot.title = element_text(size=80)) +
      scale_y_continuous(breaks=c(0,2,4,6,8,10), labels=math_format(expr = 10^.x, format = force), limits=c(0,10), expand=c(0,0))
  }
  
  Muller_df_pop$Identity = factor(Muller_df_pop$Identity, levels= myorder)
  a=ggplot(Muller_df_pop, aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
    geom_area(data=Muller_df_pop, mapping=aes_string(x = "Generation", y = "Population", group = "Group_id", fill = "Identity", colour='Identity')) + 
    #guides(linetype = FALSE, color = FALSE) + 
    ylab("Bacteria Population (/ml)") +
    xlab("Time (days)") +
    theme_classic() +scale_fill_manual(values=colors) +ggtitle(r)+scale_color_manual(values=colors) +
    xlim(0,4) +
    geom_area(aes(x=Generation, y=lines), fill=rgb(255, 255 ,255, max=255, alpha = 0), color="white", size=0.1)+ 
    #geom_point(aes(x,y), fill="white", color='white') +
    theme(plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),legend.position=c(1.08,1.08),axis.text=element_text(size=25),axis.title=element_text(size=12), axis.text.x=element_text(size=12), axis.text.y=element_text(size=10), plot.title = element_text(size = 15, hjust=0.06, vjust=-5)) +
    scale_y_continuous(breaks=c(0,2,4,6,8), labels=math_format(expr = 10^.x, format = force), limits=c(0,9), expand=c(0,0)) + labs(fill='Strain')
  
  myplots[[counter]]=a
  counter=counter+1
  
  outfile=paste0("../steps/muller/plots/alegend_muller_",r)
  pdf(file = outfile, width = 6.5, height = 5)
  print(a)
  dev.off()
  
  if (!grepl('C',r)){
    for (t in c(0,1,2,3,4)){
      variant=read.csv(paste0("../steps/muller/variants/",r1,"_variant.csv"))
      p<-ggplot(variant[variant$Time==t,], aes(x=Genotype, y=f_variant, fill=Genotype, colour=Genotype)) +
        geom_bar(stat="identity")+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=colors) + labs(title=paste0(r,'_t',t)) + ylim(0,1)+scale_color_manual(values=contour)
      #pdf(paste0("../steps/muller/variants/",r,'_',t,".pdf"),width   = 9, height = 5)
      #print(p)
      #dev.off()
    }
  }
}

#, plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"),

