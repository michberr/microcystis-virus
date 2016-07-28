library(ggplot2)
library(gtable)
library(gridExtra)

setwd("/Users/duhaimem/Box Sync/manuscripts_in_prep/PSA_phage_ecogenomics_manuscript/PSAphage_git_wd/PSA_TOV_Mapping/")

#### First, read the coverage file to make a x-y plot with x-axis -> Total cumulated coverage, and y-axis -> % of genome covered (across all TOV 90 viromes)

# import data
# Note that you may have to remove the first comma on the first line of All_TOV_90_sorted_coverage.csv, depending on your version of R
All_TOV_90_sorted_coverage <- read.csv("Coverage_files/All_TOV_90_sorted_coverage.csv")
#View(All_TOV_90_sorted_coverage)

theme<-theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"),legend.position="bottom", legend.text = element_text(size = 8))

ggplot(data=All_TOV_90_sorted_coverage) + geom_point(aes(x=coverage,y=length_covered,fill=row.names(All_TOV_90_sorted_coverage)),pch=21,size=5) + theme
### -> May be possible (and a good idea) to do some cosmetic work on this


#### Second, we read the detailed coverage file (listing all hits) to check the distribution of Id% and be able to do recruitment plots
col<-c("100_SUR"="#A2CD5A","102_SUR"="#A9A9A9","109_DCM"="#8B795E","109_SUR"="#CD3700","111_MES"="#CD2990","122_DCM"="#8B7D6B","122_MES"="#282828","122_SUR"="#9E9E9E","123_SUR"="#F8F8FF","124_SUR"="#B0E2FF","137_DCM"="#F5F5DC","137_MES"="#4A4A4A","137_SUR"="#B0B0B0","138_MES"="#6B8E23","18_DCM"="#5B5B5B","18_SUR"="#66CD00","22_SUR"="#838B8B","23_DCM"="#BDBDBD","25_SUR"="#EE9A00","30_DCM"="#00FFFF","31_SUR"="#CD3278","32_DCM"="#CDCD00","32_SUR"="#EECBAD","34_DCM"="#8E388E","36_DCM"="#292421","37_OMZ"="#969696","38_DCM"="#EE82EE","38_OMZ"="#556B2F","39_OMZ"="#D6D6D6","41_DCM"="#8B8378","56_MES"="#8B814C","62_SUR"="#BA55D3","64_DCM"="#F7F7F7","64_MES"="#CDAF95","65_DCM"="#A1A1A1","65_SUR"="#FFE7BA","66_DCM"="#9AC0CD","66_SUR"="#363636","67_SUR"="#EE4000","68_DCM"="#0F0F0F","68_MES"="#9C9C9C","68_SUR"="#9AFF9A","70_MES"="#FAFAFA","72_DCM"="#4F4F4F","72_MES"="#607B8B","76_DCM"="#FFFACD","76_MES"="#C7C7C7","78_DCM"="#CD5555","82_DCM"="#FFEC8B","82_SUR"="#EE7621","84_SUR"="#EEEEE0","85_DCM"="#3B3B3B","85_MES"="#668B8B","85_SUR"="#8B6508","Other_TOV90"="gray64")

All_genomes_recruitment <- read.csv("All_genomes_recruitment.csv")
All_genomes_recruitment <- read.csv("/mnt/DATA/Donnees/ResultatsPseudomonasphages/Final_curated_files/All_genomes_recruitment.csv")


# First the distribution of Id%
theme_2<-theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),legend.position="none")
ggplot(data=All_genomes_recruitment) + geom_boxplot(aes(x=Genome,y=Id,fill=Genome)) + theme_2

# Then we generate recruitment plots for each genome
list_genome<-unique(All_genomes_recruitment$Genome)
ng<-1
for (ng in 1:length(list_genome)){
  genome_id<-as.character(list_genome[ng])
  genome_id<-sub("/","_",genome_id)
  genome_id<-sub(",","_",genome_id)
  print(paste("generating RP for",genome_id,sep=" "))
  recruit<-All_genomes_recruitment[All_genomes_recruitment$Genome==as.character(list_genome[ng]),]
  recruit<-droplevels(recruit)
  title=paste(genome_id);
  plot_file=paste(genome_id,"_TOV_90_recruit.pdf");
  
  # Preparing the vector of colors
  list_sample<-unique(recruit$Sample)
  vec_col<-c()
  for (j in 1:length(list_sample)){
    vec_col<-c(vec_col,col[as.character(list_sample[j])])
  }
  print(vec_col)
  
  # Calculating coverage for bottom plot
  print("calculating coverage");
  max_length=max(recruit$Stop)
  coverage_all=matrix(0,max_length,1)
  for (i in 1:nrow(recruit)){
    for (j in recruit[i,2]:recruit[i,3]){
      coverage_all[j,1]<-coverage_all[j,1]+1;
    }
  }
  
  coverage<-matrix(data=NA,nrow=0,ncol=3)
  window_size_half<-max_length/50;  # Change this value to modify the size of the sliding window, larger sliding window (e.g. max_length/20) will be more smoothed, smaller window (e.g. max_length/100) will be more detailed
  step<-max_length/500; # Can lower this value (e.g. max_length/100) if the coverage computation is too long
  print(paste("max_length ",max_length,", window size half ",window_size_half,", step ",step,sep=""))
  i<-1;
  while(i<max_length){
    total<-window_size_half+window_size_half;
    boundary_start<-i-window_size_half;
    if(boundary_start<0){
      total<-total+boundary_start;
      boundary_start<-0
    }
    boundary_stop<-i+window_size_half;
    if (boundary_stop>max_length){
      total<-total-(boundary_stop-max_length);
      boundary_stop<-max_length;
    }
    cover<-sum(coverage_all[boundary_start:boundary_stop,1])/total;
    new_vec<-c(i,cover,"Total");
    coverage<-rbind(coverage,new_vec);
    i<-i+step;
  }
  
  coverage<-data.frame(coverage,row.names=NULL);
  colnames(coverage)<-c("Coord","Coverage","Sample");
  coverage$Coord<-as.numeric(as.character(coverage$Coord))
  coverage$Coverage<-as.numeric(as.character(coverage$Coverage))
  print("we have the coverage on sliding windows, let's do the plots");
  # Now doing the 3 plots (recruitment, coverage below, distribution of Id% on the right)
  n_col<-2
  # Start with the actual recruitment plot
  plot= ggplot() + scale_y_continuous(limits=c(60,101),breaks=c(60,80,100)) + ylab("% Identity")+ 
    xlab("Genomic position (bp)") + scale_x_continuous(limits=c(0,max_length)) + 
    geom_rect(data=recruit,aes(xmin=Start,xmax=Stop,ymin=Id,ymax=Id,colour=factor(Sample),name=Sample),size=0.6, alpha=0.8) + 
    guides(fill=FALSE, col = guide_legend(override.aes = list(alpha = 1, size=1 ), ncol=n_col)) + ggtitle(title) + 
    theme(title=element_text(size=8,lineheight=.8), text=element_text(size=8,lineheight=.8),legend.position = 'left',legend.justification = 'center',legend.text=element_text(size=6),legend.key.size=unit(0.5,"line"),axis.text=element_text(color="black",size=8),axis.ticks=element_line(color="black")) + scale_color_manual(name="Samples",values=vec_col);
  ggsave("plot1.pdf")
  print("plot 1 done")
  
  # Then the bottom coverage plot
  ratio.display <- 50
  ratio.values <- (max(coverage$Coord)-min(coverage$Coord))/(max(coverage$Coverage)-min(coverage$Coverage))
  plot_2 <- ggplot()  + scale_x_continuous(limits=c(0,max_length)) + 
    scale_y_log10(limits=c(0.001,100000),breaks=c(0.1,1,10,100,1000,10000)) + ylab("Log10(coverage)") +  xlab("Genomic position (bp)") +
    geom_line(data=coverage,aes(x=Coord,y=Coverage,colour=factor(Sample),name=Sample),size=0.4, alpha=0.8) + 
    guides(fill=FALSE, col = guide_legend(override.aes = list(alpha = 1, size=3 ), ncol=n_col)) + 
    theme(text=element_text(size=8,lineheight=.8),  legend.position = 'left',legend.justification = 'center',axis.text=element_text(color="black",size=8),axis.ticks=element_line(color="black")) + 
    scale_color_manual(values=c("black")) # + coord_fixed(10000)  # + coord_fixed(ratio.values/ratio.display);
  ggsave("plot2.pdf")
  print("plot 2 done");
  
  # Then the read Id% distribution
  plot_3 <- ggplot() + geom_histogram(data=recruit,aes(x=Id-1),position="identity",alpha=0.6)  + 
    scale_x_continuous(limits=c(60,101)) + coord_flip()  +  ylab("Counts") + 
    guides(fill=FALSE, col = guide_legend(override.aes = list(alpha = 1, size=3 ), ncol=1)) + ggtitle("") + 
    theme(title=element_text(size=8,lineheight=.8), text=element_text(size=8,lineheight=.8),legend.position = 'right',
          legend.justification = 'center', axis.title.y = element_blank(),axis.text=element_text(color="black",size=8),axis.ticks=element_line(color="black"));
  ggsave("plot3.pdf")
  print("plot 3 done");
  
  # Then we put all of these together as nicely as we can
  
  gp <- ggplotGrob(plot)
  gp2 <- ggplotGrob(plot_2)
  gp3 <- ggplotGrob(plot_3)
  maxWidth = grid::unit.pmax(gp$widths[2:5], gp2$widths[2:5])
  gp$widths[2:5] <- as.list(maxWidth)
  gp2$widths[2:5] <- as.list(maxWidth)
  pdf(file=plot_file, width=13, height=6,useDingbats=F)
  final_plot<-grid.arrange(gp, gp3, gp2,nrow=2,ncol=2,widths=c(8,2))
  dev.off()
  print(paste("Final plot printed in",plot_file))
}