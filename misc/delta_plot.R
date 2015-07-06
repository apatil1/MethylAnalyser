
#Function for plotting delta beta plots
delta_plot <- function(tss, sheet, Dir, type) {

  #Read the sheet consisting of questions
  dat <- read.csv(sheet, header=T)

  subDir_plot <- "Delta-Plots"
  
  #Create directory for plots and files generated from the analysis
  dir.create(file.path(Dir, subDir_plot), showWarnings = FALSE)
  
  for(i in 1:length(dat$x1)) {
    #Delta plot 

    #Calculate Dela Beta
    x <- tss[,colnames(tss)==dat$x1[i]] - tss[,colnames(tss)==dat$x2[i]]
    y <- tss[,colnames(tss)==dat$y1[i]] - tss[,colnames(tss)==dat$y2[i]]

    #Calculate Delta Beta for probes in the control having Beta values > 0.5
    red <- tss[tss[,colnames(tss)==dat$x2[i]] > 0.5, ]
    point.red.x <- red[,colnames(red)==dat$x1[i]] - red[,colnames(red)==dat$x2[i]]
    point.red.y <- red[,colnames(red)==dat$y1[i]] - red[,colnames(red)==dat$y2[i]]

    #Calculate Delta Beta for probes in the control having Beta values < 0.2
    blue <- tss[tss[,colnames(tss)==dat$x2[i]] < 0.2, ]
    point.blue.x <- blue[,colnames(blue)==dat$x1[i]] - blue[,colnames(blue)==dat$x2[i]]
    point.blue.y <- blue[,colnames(blue)==dat$y1[i]] - blue[,colnames(blue)==dat$y2[i]]

    png(paste(Dir, subDir_plot,"/", type, "_",dat$x1[i],"-",dat$x2[i],"-",dat$y1[i], "-",dat$y2[i],".jpg",sep=""), width = 1000, height = 1000)

    #plot the Delta Beta values
    plot(y~x, cex=0.4,ylim=range(-1,1),xlim=range(-1,1), col=rgb(0.6,0.6,0.6, alpha=0.4), pch=19, xlab=paste(dat$x1[i],"-",dat$x2[i],sep=""), ylab=paste(dat$y1[i],"-",dat$y2[i],sep=""))
    points(point.red.x, point.red.y, col=NULL, bg=rgb(0,0,0.702, alpha=0.4, maxColorValue = 1), pch=21, cex=0.4)
    points(point.blue.x, point.blue.y, col=NULL, bg=rgb(1,0,0, alpha=0.4, maxColorValue = 1), pch=21, cex=0.4)

    #Add lines to x and y axis
    abline(v=-0.05, col="green"); abline(v=0.05, col="green"); abline(v=-0.2, col="pink"); abline(v=0.2, col="pink")
    abline(h=0.2, col="pink"); abline(h=-0.2, col="pink")

    #Add title to the plot
    title(main=paste( type, "_", dat$x1[i],"-",dat$x2[i]," vs ",dat$y1[i],"-",dat$y2[i],sep=""), cex.main=1, col="brown")

    #Add legend to the plot
    legend("topright", inset=c(0,0), legend=c("All probes", paste(dat$x2[i]," > 0.5",sep=""), paste(dat$x2[i]," < 0.2",sep="")), col=c("grey","darkblue","red"), pt.bg=c("grey","darkblue","red"), pch=21, ncol=1, cex=1) # inset might need adjusting based on the width of the legend

    #Add counts of probes in each block in the plot
    text(-0.3,0.6,length(x[(x< -0.2) & (y > 0.2)]))
    text(-0.3,0,length(x[(x< -0.2) & ((y> -0.2) & (y< 0.2))]))
    text(-0.3,-0.6,length(x[(x < -0.2)  & (y < -0.2)]))
    text(-0.1,-0.6,length(x[(x < -0.05) & (x > -0.2)  & (y < -0.2)]))
    text(-0.1,0.6,length(x[(x < -0.05) & (x > -0.2)  & (y > 0.2)]))
    text(-0.1,0,length(x[(x < -0.05) & (x > -0.2) & (y < 0.2) & (y > -0.2)]))
    text(-0,0,length(x[(x < 0.05)  & (x > -0.05) & (y < 0.2) & (y > -0.2)]))
    text(-0,0.6,length(x[(x < 0.05)  & (x > -0.05) & (y > 0.2)]))
    text(-0,-0.6,length(x[(x < 0.05)  & (x > -0.05) & (y < -0.2)]))
    text(0.1,-0.6,length(x[(x < 0.2)  & (x > 0.05) & (y < -0.2)]))
    text(0.1,0.6,length(x[(x < 0.2)  & (x > 0.05) & (y > 0.2)]))
    text(0.1,0,length(x[(x < 0.2)  & (x > 0.05) & (y < 0.2) & (y > -0.2)]))
    text(0.3,0,length(x[(x > 0.2) & (y < 0.2) & (y > -0.2)]))
    text(0.3,0.6,length(x[(x > 0.2) & (y > 0.2)]))
    text(0.3,-0.6,length(x[(x > 0.2) & (y < -0.2) ]))

    dev.off()
  }
}