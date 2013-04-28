#! /usr/bin/env Rscript

library(ggplot2)
library(reshape)
library(grid)
library(gridExtra)

one_col_width=3.15
two_col_width=7.25

# colors
cmi_main_blue="#0071b2"
cmi_grey="#929d9e"
cmi_light_blue="#00c4d9"
cmi_pea_green="#b5bf00"

cmi_rich_green="#73933d"
cmi_rich_purple="#8e7fac"
cmi_rich_red="#d75920"
cmi_rich_blue="#4c87a1"
cmi_rich_aqua="#66c7c3"
cmi_rich_orange="#eebf42"

cmi_vibrant_yellow="#ffd457"
cmi_vibrant_orange="#f58025"
cmi_vibrant_green="#78a22f"
cmi_vibrant_garnet="#e6006f"
cmi_vibrant_purple="#9A4d9e"
cmi_vibrant_blue="#19398a"

cmi_year_colors=c(cmi_vibrant_blue,
                  cmi_rich_blue,
                  cmi_vibrant_purple,
                  cmi_vibrant_garnet,
                  cmi_rich_red,
                  cmi_vibrant_orange,
                  cmi_vibrant_yellow,
                  cmi_vibrant_green)
# TF-conditionalDF scatter plot
if (T) {
  file = paste("tfidf_processed.txt",sep="");
  x = read.table(file, as.is=T, header=T, sep="\t")
  x = x[(nrow(x)-15):nrow(x),]
  x$term <- gsub("_", " ", x$term)
  
  text.data = x
  text.data$tf = x$y
  text.data$df = x$x
  pdf(file=paste("tfidf_top.pdf",sep=""),onefile=TRUE,width=8,height=4,
      family="Times",title="tfidf",colormodel="rgb",paper="special")
  bp <- ggplot(NULL,aes(y=tf,x=df,label=term,color=Source))+
    geom_text(data=text.data) +
    geom_point(data=x) +
    theme_bw() + 
    opts(legend.position="right") +
    ylab("conditional tf") +
    xlab("df") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=10)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=10, angle=90)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=8, angle=0,vjust=0.5)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=8, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.85,0.35)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  print(bp)
  dev.off()
}

# Coauthorship graph components after removing nodes/edges
if (T) {
  all1 = read.table("connectedAfterRemovingAuthors.txt", sep="\t", as.is=T, header=T)
  all2 = read.table("connectedAfterRemovingPapers.txt", sep="\t", as.is=T, header=T)
  all1$labels = "Authors Removed"
  all2$labels = "Papers Removed"
  all = rbind(all1, all2)
  
  pdf(file="connected_after_removing.pdf",onefile=TRUE,width=9,height=6,
      family="Times",title="methods",colormodel="rgb",paper="special")
  bp <- ggplot(all,aes(y=count,x=removed,colour=labels)) + ylab("") + xlab("Quantity Removed") +
    theme_bw() + 
    scale_colour_discrete(name = "")+
    xlim(c(0,1000))+
    geom_line(size=1) + geom_point(size=3) +
    facet_grid(type~.,scale="free")+
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.84,0.865)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  print(bp)
  dev.off()
}

# Growth of methods
if (T) {
  citations = read.table("methods_growth.txt", sep="\t", as.is=T, header=T)
  citations = citations[which(citations$year>=2005),]
  pdf(file="methods_growth_hist.pdf",onefile=TRUE,width=9,height=6,
      family="Times",title="citations",colormodel="rgb",paper="special")
  bp <- ggplot(citations,aes(y=count,x=reorder(method, count, sum),fill=factor(year))) +
    theme_bw() +
    geom_bar(alpha=1.0)+labs(fill='Year')+coord_flip()+
    scale_y_continuous() +
    ylab("Count") +
    xlab("Method") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=0,vjust=0.5)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.85,0.32)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  print(bp)
  dev.off()
}

# Growth of tags
if (T) {
  x = read.table("clinical_bytag_growth.txt", as.is=T, header=T, sep="\t", allowEscapes=T, quote="\"")
  pdf(file="clinical_bytag_hist.pdf",onefile=TRUE,width=two_col_width,height=5,
      family="Times",title="Clinical Tags",colormodel="rgb",paper="special")
  bp <- ggplot(x,aes(y=count,x=reorder(Tag, count, sum),fill=factor(Year))) +
    theme_bw()+
    geom_bar()+labs(fill='Year')+
    scale_fill_manual(name="Year", values=cmi_year_colors)+
    #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
    scale_y_continuous() +
    ylab("Count") +
    xlab("Tag") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=30,vjust=1,hjust=1)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.1,0.65)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  print(bp)
  dev.off()
}

# Growth of 'connectome' usage
if (T) {
  citations = read.table("connectome_growth.txt", sep="\t", as.is=T, header=T)
  citations = citations[which(citations$year>=2005),]
  pdf(file="connectome_growth.pdf",onefile=TRUE,width=9,height=4,
      family="Times",title="methods",colormodel="rgb",paper="special")
  bp <- ggplot(citations,aes(y=count,x=year,colour=metric)) +labs(colour='Metric')+
    #geom_bar(fill="#660000",alpha=.5)+
    theme_bw()+
    geom_line(size=1) + geom_point(size=3) +
    #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
    scale_y_continuous() +
    ylab("Value") +
    xlab("Year") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.15,0.8)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  
  print(bp)
  dev.off()
}

# Journal counts by year
if (T) {
  citations = read.table("journal_counts_byear.txt", sep="\t", as.is=T, header=T)
  pdf(file="journal_dist.pdf",onefile=TRUE,width=two_col_width,height=5,
      family="Times",title="growth_rate_journal",colormodel="rgb",paper="special")
  bp <- ggplot(citations,aes(y=count,x=reorder(journal, count, sum),fill=factor(year))) +
    theme_bw() +
    scale_fill_manual(name="Year", values=cmi_year_colors) +
    geom_bar(alpha=1.0)+
    labs(fill='Year')+
    scale_x_discrete(expand=c(0.01,0.01))+
    scale_y_continuous(expand=c(0.01,0.01))+
    ylab("Publication Count") +
    xlab("Journal") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12, hjust=0.78, vjust=0)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90,vjust=.3)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=30,vjust=1,hjust=1)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
    opts(legend.background = theme_rect(fill = 'white', size = 0, colour='white', linetype='dashed')) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(legend.position=c(0.1,0.65)) +
    opts(panel.margin = unit(c(0.0, 0.0,0.0,0.0), "lines"))+
    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
  
  #coord_flip()+
  print(bp)
  dev.off()
  
}

# Overall corpus growth
if (T) {
  cutoff = 2006
  startYear = 2000
  x = read.table("overall_growth.txt", as.is=T, header=T, sep="\t")
  x = x[which(x$year>1994),]
  yearFreqs = x
  yearFreqs=melt(yearFreqs,id="year")
  yearFreqs$year=yearFreqs$year
  yearFreqs$lvalue=log(yearFreqs$value)
  yearFreqs$lvalue[is.infinite(yearFreqs$lvalue)]=0
  
  yearFreqs$labels=yearFreqs$variable
  levels(yearFreqs$labels)=c("fMRI","Resting State")
  
  yearFreqs_rs=subset(yearFreqs,variable=="rsTotal",drop=TRUE)
  yearFreqs_rs$vol=cumsum(yearFreqs_rs$value)
  yearFreqs_pre2005_rs=subset(yearFreqs_rs,(year>startYear)&(year<=cutoff),drop=TRUE)
  yearFreqs_post2004_rs=subset(yearFreqs_rs,year>=cutoff,drop=TRUE)
  
  rs_model_all=lm(log(vol)~year,data=yearFreqs_rs)
  rs_model_pre2005=lm(log(vol)~year,data=yearFreqs_pre2005_rs)
  rs_model_post2004=lm(log(vol)~year,data=yearFreqs_post2004_rs)
  
  print(rs_model_pre2005)
  print(rs_model_post2004)
  
  yearFreqs$model_fit[yearFreqs$variable=="rsTotal"]=exp(fitted(rs_model_all))
  yearFreqs$model_fit_pre2005[(yearFreqs$variable=="rsTotal")&(yearFreqs$year<=cutoff)&(yearFreqs$year>startYear)]=exp(fitted(rs_model_pre2005))
  yearFreqs$model_fit_post2004[(yearFreqs$variable=="rsTotal")&(yearFreqs$year>=cutoff)]=exp(fitted(rs_model_post2004))
  yearFreqs$vol[yearFreqs$variable=="rsTotal"]=yearFreqs_rs$vol
  
  yearFreqs_neuro=subset(yearFreqs,variable=="neuroTotal",drop=TRUE)
  yearFreqs_neuro$vol=cumsum(yearFreqs_neuro$value)
  yearFreqs_pre2005_neuro=subset(yearFreqs_neuro,(year>startYear)&(year<=cutoff),drop=TRUE)
  yearFreqs_post2004_neuro=subset(yearFreqs_neuro,year>=cutoff,drop=TRUE)
  neuro_model_all=lm(log(vol)~year,data=yearFreqs_neuro)
  neuro_model_pre2005=lm(log(vol)~year,data=yearFreqs_pre2005_neuro)
  neuro_model_post2004=lm(log(vol)~year,data=yearFreqs_post2004_neuro)
  
  print(neuro_model_pre2005$residuals)
  print(neuro_model_pre2005)
  print(neuro_model_post2004)
  print(neuro_model_all)
  
  yearFreqs$model_fit[yearFreqs$variable=="neuroTotal"]=exp(fitted(neuro_model_all))
  yearFreqs$model_fit_post2004[(yearFreqs$variable=="neuroTotal")&(yearFreqs$year>=cutoff)]=exp(fitted(neuro_model_post2004))
  yearFreqs$model_fit_pre2005[(yearFreqs$variable=="neuroTotal")&(yearFreqs$year>startYear)&(yearFreqs$year<=cutoff)]=exp(fitted(neuro_model_pre2005))
  
  
  yearFreqs$vol[yearFreqs$variable=="neuroTotal"]=yearFreqs_neuro$vol
  yearFreqs$pre_2005_vol[yearFreqs$year<=cutoff]=yearFreqs$vol[yearFreqs$year<=cutoff]
  yearFreqs$post_2004_vol[yearFreqs$year>=cutoff]=yearFreqs$vol[yearFreqs$year>=cutoff]
  
  p=ggplot(yearFreqs)+
    theme_bw()+
    geom_area(aes(x=year,y=pre_2005_vol),fill=cmi_grey,alpha=.7)+
    geom_line(aes(x=year,y=model_fit_pre2005),color=cmi_grey)+
    geom_area(aes(x=year,y=post_2004_vol),fill=cmi_main_blue,alpha=.7)+
    geom_line(aes(x=year,y=model_fit_post2004),color=cmi_main_blue)+
    facet_grid(labels~.,scale="free")+
    ylab("Paper Volume") +
    xlab("Year") +
    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
    opts(strip.text.x = theme_text(family = "Times", face = "plain", size=12)) +
    opts(strip.text.y = theme_text(family = "Times", face = "plain", size=12, angle=270)) +
    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10)) +
    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=90)) +
    opts(axis.ticks.length = unit(.15, "lines")) +
    opts(axis.ticks.margin=unit(.15,"lines")) +
    opts(plot.margin = unit(c(0.75, 0.25,0.5,0.75), "lines"))+
    opts(strip.background=element_blank())
  
  pdf(file="overall_growth.pdf",onefile=TRUE,width=one_col_width,height=4,
      family="Times",title="Overall Growth",colormodel="rgb",paper="special")
  print(p)
  dev.off()
}

# Comparison of 'connectome' usage between DTI and RS libraries
if (T) {
        all = read.table("dti_growth.txt", sep="\t", as.is=T, header=T)
	connect = all[which(all$Tag=="Connectome"),]
	all = all[which(all$Tag=="All"),]
        pdf(file="dti_vs_rs_all_growth.pdf",onefile=TRUE,width=9,height=6,
                family="Times",title="methods",colormodel="rgb",paper="special")
        bp <- ggplot(all,aes(y=count,x=Year,colour=Library)) + ylab("Count") +
            #geom_bar(fill="#660000",alpha=.5)+
            geom_line(size=1) + geom_point(size=3) +
            #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
            scale_y_continuous() +
            xlab("Year") +
	    
            opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
            opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
            opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
            opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
            opts(axis.ticks.length = unit(.15, "lines")) +
            opts(axis.ticks.margin=unit(.15,"lines")) +
            opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
	    print(bp)
        dev.off()
        pdf(file="dti_vs_rs_connectome_growth.pdf",onefile=TRUE,width=9,height=6,
                family="Times",title="methods",colormodel="rgb",paper="special")
        bp <- ggplot(connect,aes(y=count,x=Year,colour=Library)) + ylab("Count") +
            #geom_bar(fill="#660000",alpha=.5)+
            geom_line(size=1) + geom_point(size=3) +
            #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
            scale_y_continuous() +
            xlab("Year") +
	    
            opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
            opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
            opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
            opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
            opts(axis.ticks.length = unit(.15, "lines")) +
            opts(axis.ticks.margin=unit(.15,"lines")) +
            opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
	    print(bp)
        dev.off()
}

# Impact factor by author
if (T) {
	citations = read.table("impactByAuthor.txt", as.is=T, header=T, sep="\t", allowEscapes=T)
	pdf(file="impactByAuthor.pdf",onefile=TRUE,width=9,height=6,
   		family="Times",title="citations",colormodel="rgb",paper="special")
	bp <- ggplot(citations,aes(y=impactFactor,reorder(author, impactFactor))) +
	    geom_bar(fill="#660000",alpha=.5)+
	    #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
	    scale_y_continuous() +
	    ylab("Impact Factor") +
	    xlab("Author") +
	    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
	    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
	    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
	    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
	    opts(axis.ticks.length = unit(.15, "lines")) +
	    opts(axis.ticks.margin=unit(.15,"lines")) +
	    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
	print(bp)
   	dev.off()
}

# Citations by author
if (T) {
	citations = read.table("citationsByAuthor.txt", as.is=T, header=T, sep="\t", allowEscapes=T)
	pdf(file="citationsByAuthor.pdf",onefile=TRUE,width=9,height=6,
   		family="Times",title="citations",colormodel="rgb",paper="special")
	bp <- ggplot(citations,aes(y=citationCount,reorder(author,citationCount))) +
	    geom_bar(fill="#660000",alpha=.5)+
	    #geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
	    scale_y_continuous() +
	    ylab("Citations") +
	    xlab("Author") +
	    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=12)) +
	    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=12, angle=90)) +
	    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=10, angle=45,vjust=0.5)) +
	    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=10, angle=0,hjust=1)) +
	    opts(axis.ticks.length = unit(.15, "lines")) +
	    opts(axis.ticks.margin=unit(.15,"lines")) +
	    opts(plot.margin = unit(c(0.75, 0.75,0.75,0.75), "lines"))
	print(bp)
   	dev.off()
}

if (T) {
	citations = read.table("bag_report_formatted.txt", as.is=T, header=T, sep=",",quote="")
  print(citations$title)
	citations$title <- sub("!", "\n", citations$title)

	pdf(file="citation_pageranks.pdf",onefile=TRUE,width=9,height=6,
   		family="Times",title="citations",colormodel="rgb",paper="special")
	bp <- ggplot(citations,aes(y=mean,x=reorder(title,mean)))+
	    theme_bw()+
	    geom_bar(fill="#660000",alpha=0.7)+
	    geom_errorbar(aes(ymin = mean - stddev, ymax = mean + stddev),color="#660000",width=0.5)+
	    coord_flip()+
	    scale_y_continuous() +
	    ylab("Pagerank") +
	    xlab("Title") +
	    opts(axis.title.x = theme_text(family = "Times", face = "plain", size=14, hjust=0.85)) +
	    opts(axis.title.y = theme_text(family = "Times", face = "plain", size=14, angle=90)) +
	    opts(axis.text.x  = theme_text(family = "Times", face = "plain", size=12, angle=0,vjust=0.5)) +
	    opts(axis.text.y  = theme_text(family = "Times", face = "plain", size=12, angle=0,hjust=1)) +
	    opts(axis.ticks.length = unit(.15, "lines")) +
	    opts(axis.ticks.margin=unit(.15,"lines")) +
	    opts(plot.margin = unit(c(0.75, 0.25,0.25,0.75), "lines"))
	print(bp)
   	dev.off()
}
