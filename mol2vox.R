source('util.R')

xyz2vox <- function(xyz, dims){
  center <- (dims-1)/2
  dims.names <- list(1:dims[1]-1, 1:dims[2]-1, 1:dims[3]-1)
  block <- array(rep('', prod(dim)), dim=dims,  dimnames=dims.names)
  xyz <- recenterXYZ(xyz, center)
  placeAtoms(xyz, block)
}
# xyz <- smi2xyz('C1CCCCC1')
# dims <- c(10,10,10)
# xyz2vox(scaleXYZ(xyz, 3), dims)
