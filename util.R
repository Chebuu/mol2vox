library(ChemmineOB)

argmin.1 <- function(X){
  # TODO::
  # - !!! TIES ARE RESOLVED BY SELECTING THE FIRST INDEX
  # # - Rounding should be done instead
  which(X == min(X))[1]
}

smi2xyz <- function(smi){
  convertFormat('smi','xyz', smi)
}

mol2xyz <- function(mol){
  convertFormat(from='mol', to='xyz', mol)
}

xyz2mat <- function(xyz){
  path.xyz <- tempfile('temp', fileext='.xyz')
  writeLines(xyz, path.xyz)
  df <- read.table(path.xyz, header=FALSE, quote=NULL, sep='\n', numerals='no.loss', check.names=FALSE, stringsAsFactors=FALSE)
  df <- do.call(rbind, apply(df, 1, function(r){ trimws(strsplit(r, '       ')[[1]]) }))
  df[1,c(1,2,3)] <- ''
  file.remove(path.xyz)
  return(df)
}
# xyz <- xyz2mat(smi2xyz('C1CCCCC1'))

dimRange <- function(xyz){
  if(!is.matrix(xyz)){
    xyz <- xyz2mat(xyz)
  }
  
  xyz.coords <- xyz[2:nrow(xyz), 2:ncol(xyz)]

  dims.max <- apply(xyz.coords, 2, function(dim){
    max(as.numeric(dim))
  })
  names(dims.max) <-  c('x','y','z')
  
  dims.min <- apply(xyz.coords, 2, function(dim){
    min(as.numeric(dim))
  })
  names(dims.min) <- c('x','y','z')
  
  t(data.frame(
    min=dims.min,
    max=dims.max
  ))
}

dimRanges <- function(xys.list){
  dim.min <- lapply(xyz.list, function(xyz){
    dimRange(xyz)['min',]
  })
  dim.min <- do.call(rbind,dim.min)
  dim.min <- apply(dim.min, 2, min)
  
  dim.max <-lapply(xyz.list, function(xyz){
    dimRange(xyz)['max',]
  })
  dim.max <- do.call(rbind, dim.max)
  dim.max <- apply(dim.max, 2, max) 
  
  data.frame(
    min=dim.min,
    max=dim.max
  )
}

zeroXYZ <- function(xyz){
  if(!is.matrix(xyz)) xyz <- xyz2mat(xyz)
  coords <- apply(xyz[-(1), -(1)], 2, as.numeric)
  mins <- dimRange(xyz)[1,]
  trans <- t(sapply(1:nrow(coords), function(i){ mins }))
  coords.trans <- coords - trans
  xyz[-(1), -(1)] <- coords.trans
  return(xyz)
}
# library(scatterplot3d)
# xyz <- xyz2mat(smi2xyz('C1CCCCC1'))
# xyz.coords <- xyz[-(1), -(1)]
# xyz.zeroed <- zeroXYZ(xyz)[-(1), -(1)]
# par(mfrow=c(1,2))
# scatterplot3d(xyz.coords[,1],xyz.coords[,2],xyz.coords[,3], xlim=c(-2,1),ylim=c(-1,3),zlim=c(0,1))
# scatterplot3d(xyz.zeroed[,1],xyz.zeroed[,2],xyz.zeroed[,3], xlim=c(0,3),ylim=c(-0,3),zlim=c(0,1))

recenterXYZ <- function(xyz, center){
  if(!is.matrix(xyz)) xyz <- xyz2mat(xyz)
  xyz.coords <- apply(xyz[-(1), -(1)], 2, as.numeric)
  xyz.dims <- dimRange(xyz)
  xyz.center <- xyz.dims[1,] + (xyz.dims[2,] - xyz.dims[1,]) / 2
  trans <- as.vector(center - xyz.center)
  xyz.recenter <- t(t(xyz.coords) + trans)
  xyz[-(1), -(1)] <- xyz.recenter
  return(xyz)
}
# library(scatterplot3d)
# xyz <- xyz2mat(smi2xyz('C1CCCCC1'))
# xyz.coords <- xyz[-(1), -(1)]
# xyz.recentered <- recenterXYZ(xyz, c(2,2,2))[-(1), -(1)]
# par(mfrow=c(1,2))
# scatterplot3d(xyz.coords[,1],xyz.coords[,2],xyz.coords[,3], xlim=c(-3,3),ylim=c(-3,3),zlim=c(-3,3))
# scatterplot3d(xyz.recentered[,1],xyz.recentered[,2],xyz.recentered[,3], xlim=c(-3,3),ylim=c(-3,3),zlim=c(-3,3))

scaleXYZ <- function(xyz, scale){
  if(!is.matrix(xyz)) xyz <- xyz2mat(xyz)
  xyz.coords <- apply(xyz[-(1), -(1)], 2, as.numeric)
  xyz.scaled <- xyz.coords * scale
  xyz[-(1), -(1)] <- xyz.scaled
  return(xyz)
}
# library(scatterplot3d)
# xyz <- xyz2mat(smi2xyz('C1CCCCC1'))
# xyz.coords <- xyz[-(1), -(1)]
# xyz.scaled <- scaleXYZ(xyz, 2.8)[-(1), -(1)]
# par(mfrow=c(1,2))
# scatterplot3d(xyz.coords[,1],xyz.coords[,2],xyz.coords[,3], xlim=c(-5,5),ylim=c(-5,5),zlim=c(-5,5))
# scatterplot3d(xyz.scaled[,1],xyz.scaled[,2],xyz.scaled[,3], xlim=c(-5,5),ylim=c(-5,5),zlim=c(-5,5))

nearestVox <- function(coords, dims){
  # @ param coords c(<x>, <y>, <z>)
  # @ param block c(<x_lim>, <y_lim>, <z_lim>)
  setNames(sapply(1:3, function(i){
    drange <- dims[[i]]
    didx <- argmin.1(abs(drange - coords[[i]])) # TODO:: TIES ARE NOT RESOLVED !!!
    drange[[didx]]
  }), c('x', 'y', 'z'))
}
# nearestVox(c(1,2,3), rep(seq(1,10,3),3)) # TODO:: TIES ARE NOT RESOLVED !!!

placeAtoms <- function(xyz, block){
  if(!is.matrix(xyz)) xyz <- xyz2mat(xyz)
  block.dims <- lapply(dimnames(block), as.numeric)
  atoms <- xyz[-(1), 1]
  coords <- apply(xyz[-(1), -(1)], 2, as.numeric)
  
  if(any(unlist(coords) < 0)){
    # xyz must be translated to at least 0,0,0
    coords <- apply(zeroXYZ(xyz)[-(1), -(1)], 2, as.numeric)
  }
  
  for(n in 1:nrow(coords)){
    atom.x <- coords[n,1]
    atom.y <- coords[n,2]
    atom.z <- coords[n,3]
    
    atom.coords <- c(atom.x, atom.y, atom.z)
    atom.vox <- nearestVox(atom.coords, block.dims)

    block[
      as.character(atom.vox[[1]]), 
      as.character(atom.vox[[2]]), 
      as.character(atom.vox[[3]])
    ] <- atoms[[n]]
  }
  return(block)
}
# xyz <- recenterXYZ(smi2xyz('C1CCCCC1'), c(2,2,2))
# dims <- c(5,5,5)
# dnames <- list(1:dims[1]-1, 1:dims[2]-1, 1:dims[3]-1)
# block <- array(rep('', prod(dims)), dim=dims, dimnames=dnames)
# placeAtoms(xyz, block)
