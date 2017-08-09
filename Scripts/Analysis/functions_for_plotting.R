library(plot3Drgl)

mean.expression.plot <- function(mean.coords, all.coords, col.var, cl, 
                                 brain.cl, col.var.title, add.cluster=NULL, cluster.col=NULL, theta=25, phi=25,
                                 xlab="x", ylab="y", zlab="z", bty="o"){
  ### Function to plot mean coords with color based on gene expression
  # Parameters
  # mean.coords: Matrix of mean voxel positions to be plotted
  # all.coords: Determines range of plot
  # col.var: colvar, should have the same ordering as mean.coords
  # cl: colors to be used
  # brain.cl: color of brain voxels
  # col.var.title: Title to be displayed above colkey
  # optional: add.cluster: vector of voxels to be colored differently
  # optional: cluster.col: color of the added cluster
  
  ##### Function in function
  #panelfirst <- function(pmat){
    ### Function to generate panels before plotting
    # Parameters:
    # mean.coords: Matrix of mean voxel positions to be plotted
    # all.coords: Determines range of plot
    # col.var: colvar, should have the same ordering as mean.coords
    # cl: colors to be used
    #zmin <- min(all.coords[,3])
    #XY <- trans3D(x = mean.coords[1,], 
    #              y = mean.coords[2,],
    #              z = rep(zmin, ncol(mean.coords)), 
    #              pmat=pmat
    #)
    #scatter2D(XY$x, XY$y, colvar=col.var, 
    #          pch=16, cex=0.5, add=TRUE, 
    #          colkey=FALSE, col = cl, 
    #          xlim = c(min(all.coords[,1]), max(all.coords[,1])),
    #          ylim = c(min(all.coords[,2]), max(all.coords[,2]))
    #)
    
    #ymax <- max(all.coords[,2])
    #XY <- trans3D(x = mean.coords[1,], 
     #             y = rep(ymax, ncol(mean.coords)),
    #              z = mean.coords[3,],
    #              pmat=pmat
    #)
    #scatter2D(XY$x, XY$y, colvar=col.var, 
    #          pch=16, cex=0.5, add=TRUE, 
    #          colkey=FALSE, col = cl,
    #          xlim = c(min(all.coords[,1]), max(all.coords[,1])),
    #          ylim = c(min(all.coords[,3]), max(all.coords[,3]))
    #)
    
    #xmin <- min(all.coords[,1])
    #XY <- trans3D(x = rep(xmin, ncol(mean.coords)), 
    #              y = mean.coords[2,],
    #              z = mean.coords[3,],
    #              pmat=pmat
    #)
    #scatter2D(XY$x, XY$y, colvar=col.var, 
    #          pch=16, cex=0.5, add=TRUE, 
    #          colkey=FALSE, col = cl,
    #          xlim = c(min(all.coords[,2]), max(all.coords[,2])),
    #          ylim = c(min(all.coords[,3]), max(all.coords[,3]))
    #)
  #}
  
  
  coords.for.plotting <- cbind(all.coords, rep((min(col.var)-2), length(all.coords[,1])))
  m <- t(rbind(mean.coords, col.var))
  colnames(m) <- colnames(coords.for.plotting) <- c("X", "Y", "Z", "Col Var")
  coords.for.plotting <- rbind(coords.for.plotting, m)
  
  if (!is.null(add.cluster)){
    coords.for.plotting[add.cluster, 4] <- ((min(col.var)-1))
    col.seq = c(brain.cl, cluster.col, cl)
    size <- ifelse(coords.for.plotting[,4]==(min(col.var)-2), 0.7, ifelse(coords.for.plotting[,4]==(min(col.var)-1), 0.7, 1))
  }
  else{
    col.seq = c(brain.cl, cl)
    size <- ifelse(coords.for.plotting[,4]==(min(col.var)-2), 0.7, 1)
  }
  
  scatter3D(x = coords.for.plotting[,1], y = coords.for.plotting[,2], 
            z = coords.for.plotting[,3], 
            colvar = coords.for.plotting[,4], col = col.seq,
            pch = 16, cex = size, bty=bty,
            #panel.first = panelfirst,
            xlab=xlab, ylab=ylab, zlab=zlab,
            theta=theta, phi=phi,
            colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
            clab = col.var.title
  )
}

###### Plot all voxels per cell

all.voxels.expression.plot <- function(cell.coords, all.coords, col.var, cl, 
                                       brain.cl, col.var.title, t, p){
  
  all.voxels <- as.matrix(cbind(all.coords, rep((min(col.var)-0.1), nrow(all.coords))))
  colnames(all.voxels) <- c("X", "Y", "Z", "Col.Var")
  
  #for(i in 1:length(keys(cell.coords))){
  #  all.voxels[unlist(lapply(cell.coords[[keys(cell.coords)[i]]], rownames)), 4] <- col.var[keys(cell.coords)[i]]
  #}
  
  for(i in 1:length(keys(cell.coords))){
    all.voxels[rownames(cell.coords[[keys(cell.coords)[i]]]), 4] <- col.var[keys(cell.coords)[i]]
  }
  
  scatter3D(all.voxels[,1], 
            all.voxels[,2], 
            all.voxels[,3],
            colvar=all.voxels[,4],
            col = c(brain.cl, cl),
            pch = 16, cex = ifelse(all.voxels[,4]==(min(col.var)-0.1), 0.7, 1.0),
            theta = t, phi = p,
            colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75),
            clab = col.var.title, xlab=xlab, ylab=ylab, zlab=zlab        
  )
}
