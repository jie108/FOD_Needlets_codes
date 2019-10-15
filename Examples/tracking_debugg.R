rm(list=ls())
library(R.matlab)
library(rgl)
library(compositions)
source("dwi_fit.R")
source("dwi_track.R")


path_load = '/Users/hao/Dropbox/stats_project/FOD_codes_simulation/Real_data/S110933/fitting/space_indexx108-123y124-139z37-42/'
num_fib_cut = 4
temp = readMat(paste0(path_load,'for_tracking_cut',toString(num_fib_cut),'.mat'))


v.obj = temp
eig = -v.obj$vec
loc = v.obj$loc

tracks1 <- list()
tracks2 <- list()
all.pvox <- NULL
all.pdir <- NULL
all.pdis <- NULL
all.ppdis <- NULL
n.use.iind <- array(0, dim=length(v.obj$n.fiber2))
n.iinds <- array(0,dim=length(v.obj$n.fiber2))
lens <- array(0, dim=length(v.obj$n.fiber2))

braingrid = temp$braingrid
xgrid.sp = temp$xgrid.sp
ygrid.sp = temp$ygrid.sp
zgrid.sp = temp$zgrid.sp
map = temp$map
rmap = temp$rmap
n.fiber = temp$n.fiber
n.fiber2 = temp$n.fiber2

max.line = 100
nproj = 2
thres.ang = 0.5235988
vorient=c(1,1,1)
elim = T
elim.thres = 1



for (iind in which(v.obj$n.fiber2>0)){
  cat(iind,"\n")
  tracks1[[iind]] <- fiber.track(iind=iind, eig=v.obj$vec, loc=v.obj$loc,
                                 map=v.obj$map, rmap=v.obj$rmap,
                                 n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                                 ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp, braingrid=braingrid,
                                 max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                                 vorient=vorient)
  
  tracks2[[iind]] <- fiber.track(iind=iind, eig=-v.obj$vec, loc=v.obj$loc,
                                 map=v.obj$map, rmap=v.obj$rmap,
                                 n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                                 braingrid=braingrid,
                                 ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp,
                                 max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                                 vorient=vorient)
  
  #all.pvox <- c(all.pvox, tracks1[[iind]]$pvox, tracks2[[iind]]$pvox)
  #all.pdir <- rbind(all.pdir, tracks1[[iind]]$pdir, tracks2[[iind]]$pdir)
  #all.pdis <- c(all.pdis, tracks1[[iind]]$pdis, tracks2[[iind]]$pdis)
  #all.ppdis <- c(all.ppdis, tracks1[[iind]]$ppdis, tracks2[[iind]]$ppdis)
  n.use.iind[tracks1[[iind]]$iinds] <- n.use.iind[tracks1[[iind]]$iinds] + 1
  n.use.iind[tracks2[[iind]]$iinds] <- n.use.iind[tracks2[[iind]]$iinds] + 1
  n.use.iind[iind] <- n.use.iind[iind] - 1
  n.iinds[iind] <- length(union(tracks1[[iind]]$iinds, tracks2[[iind]]$iinds))
  lens[iind] <- get.fdis(tracks1[[iind]]$inloc) + get.fdis(tracks2[[iind]]$inloc)
  
  if (length(all.pdis)!=length(all.pvox)){
    break
  }
}
#len.ord <- order(n.iinds, decreasing=T)
len.ord <- order(lens, decreasing=T)
if (max(lens[n.iinds<=1])> elim.thres){
  cat("elim.thres is too small: it should be set at least", max(lens[n.iinds<=1]),"\n") 
}

if (elim){
  update.ind <- rep(T, length(v.obj$n.fiber2))
  #update.ind[as.logical((v.obj$n.fiber2==0)+(n.use.iind<=elim.thres))] <- F
  #update.ind[as.logical((v.obj$n.fiber2==0)+(n.iinds<=elim.thres))] <- F
  update.ind[as.logical((v.obj$n.fiber2==0)+(lens<=elim.thres))] <- F
  nv.obj <- update.v.obj(v.obj, list(vec=v.obj$vec, update.ind=update.ind))$obj
} else {
  nv.obj <- v.obj
  update.ind <- rep(T, length(v.obj$n.fiber2))
  #update.ind[as.logical((v.obj$n.fiber2==0)+(n.iinds<=elim.thres))] <- F
  update.ind[as.logical((v.obj$n.fiber2==0)+(lens<=elim.thres))] <- F
}
sorted.iinds <- (1:length(v.obj$n.fiber2))[len.ord]
sorted.update.ind <- update.ind[len.ord]


#############################
## v.track
#############################

idx_fiber_track = which(v.obj$n.fiber2>0)
iind = idx_fiber_track[1472]
iind = 19
cat(iind,"\n")
tracks1[[iind]] <- fiber.track(iind=iind, eig=v.obj$vec, loc=v.obj$loc,
                               map=v.obj$map, rmap=v.obj$rmap,
                               n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                               ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp, braingrid=braingrid,
                               max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                               vorient=vorient)

tracks2[[iind]] <- fiber.track(iind=iind, eig=-v.obj$vec, loc=v.obj$loc,
                               map=v.obj$map, rmap=v.obj$rmap,
                               n.fiber=v.obj$n.fiber, xgrid.sp=xgrid.sp,
                               ygrid.sp=ygrid.sp, zgrid.sp=zgrid.sp, braingrid=braingrid,
                               max.line=max.line, nproj=nproj, thres.ang=thres.ang,
                               vorient=vorient)

#############################
## fiber.track
#############################
braindim <- dim(braingrid)[-1]
nvox <- prod(braindim)
dimens <- c(xgrid.sp, ygrid.sp, zgrid.sp)

path.voxel <- array(dim=max.line)
path.dir <- array(dim=c(max.line, 3))
path.in <- array(dim=c(max.line, 3))
path.change <- array(dim=max.line)
path.iind <- array(dim=max.line)
pass.vox <- NULL
pass.dir <- NULL
pass.dis <- NULL
pass.pdis <- NULL # perpendicular distance

iind = 19

# initialization
path.voxel[1] <- map[iind]
path.dir[1,] <- eig[iind,]
path.in[1,] <- loc[iind,]
path.change[1] <- T
path.iind[1] <- iind

ii <- 1
while ((ii<max.line)){
  #     if (T){
  #       cat(ii,"\n")
  #       spheres3d(path.in[ii,], radius=0.002, col="red")
  #     }
  
  # fio <- fiber.in.out(inc=path.in[ii,]-loc[path.iind[ii],], direct=path.dir[ii,], dimens=dimens)
  
  inc = path.in[ii,]-loc[path.iind[ii],]
  direct=path.dir[ii,]
  if (sum(dimens==0)){
    stop("directions has zero component, not yet supported! Please modify fiber.in.out\n")
  }
  # compute the distance of the current fiber directon to each face of the current voxel
  tempdiff <- (round(cbind(dimens/2-inc,-inc-dimens/2),5)/direct)  ## Hao: add round5
  cbind(dimens/2-inc,-inc-dimens/2)
  tempdiff
  # Hao
  # tempdiff[tempdiff==Inf]=1e10
  tempdiff[tempdiff==-Inf]=Inf
  # tempdiff[is.nan(tempdiff)]=1e10
  # tempdiff[tempdiff==Inf]=1e10
  
  
  # get which axis is the current fiber direction hitting face of the current voxel first
  # 1:x  2:y  3:z
  index1 <- which.min(diag(tempdiff[,2-(direct>=0)]))  # Hao change direct>0 to direct>=0
  # which direction it is hitting 1:positive  2:negative
  index <- c(index1, (2-(direct>0))[index1])
  const <- tempdiff[index[1],index[2]]
  outc <- round(inc + const*direct,5)     ## Hao: add round5
  
  fio = list(outc=outc,index=as.vector(index))
  
  path.in[ii+1,] <- fio$outc + loc[path.iind[ii],]
  
  # for previous pass.dis and pass.pdis, using the previous "change"
  if ((!path.change[ii])&&(n.fiber[path.voxel[ii]]>0)){
    pass.pdis <- c(pass.pdis, dist.line(loc[path.iind[ii],], path.in[ii,], path.in[ii+1,]))
    pass.dis <- c(pass.dis, sqrt(sum((path.in[ii,]-path.in[ii+1,])^2)))
  }
  
  # determine which voxel it is going to
  next.vox <- get.out.vox(fio$index, path.voxel[ii], braindim=braindim, vorient=vorient)
  
  if (is.na(next.vox)){
    break
  }
  
  # determine if we should stop
  pro.res <- project.proceed(inc0=path.in[ii+1,], vox0=next.vox,
                             dir0=path.dir[ii,], loc, eig, rmap, n.fiber,
                             braindim, dimens, nproj=nproj,
                             thres.ang=thres.ang, vorient=vorient)
  change <- pro.res$first
  good <- pro.res$last
  
  if (!good){
    break
  }
  
  # update voxel
  path.voxel[ii+1] <- next.vox
  
  # update dir, iind and change
  if (n.fiber[next.vox]<=1){
    path.iind[ii+1] <- rmap[next.vox]
    path.change[ii+1] <- change
    if (change){
      path.dir[ii+1,] <- eig[path.iind[ii+1],]
    } else {
      path.dir[ii+1,] <- path.dir[ii,]
      if (n.fiber[next.vox]==1){
        pass.vox <- c(pass.vox,next.vox)
        pass.dir <- rbind(pass.dir, path.dir[ii,])
      }
    }
  } else {
    # thresholding rule -> determine stop or not, and within the thresholding rule, choose the closest
    
    if (change){
      # decide which directions
      tiind <- rmap[next.vox]
      chosen <- which.max(abs(eig[tiind+(0:(n.fiber[next.vox]-1)),]%*%path.dir[ii,]))
      path.iind[ii+1] <- tiind+chosen-1
      path.dir[ii+1,] <- eig[path.iind[ii+1],]
      path.change[ii+1] <- T
    } else {
      path.iind[ii+1] <- rmap[next.vox]
      path.change[ii+1] <- F
      path.dir[ii+1,] <- path.dir[ii,]
      pass.vox <- c(pass.vox,next.vox)
      pass.dir <- rbind(pass.dir, path.dir[ii,])
    }
  }
  
  # align directions
  path.dir[ii+1,] <- sign(sum(path.dir[ii+1,]*path.dir[ii,]))*path.dir[ii+1,]
  
  ii <- ii+1
}

if (ii<max.line){
  path.in <- path.in[1:(ii+1),]
  path.iind <- path.iind[1:ii]
  path.dir <- path.dir[1:ii,]
  path.change <- path.change[1:ii]
}


#############################
## project.proceed
#############################
vox0 = next.vox
dir0=path.dir[ii,]
first <- proceed(vox0, dir0, eig, rmap, n.fiber, thres.ang)

#############################
## proceed
#############################
good <- T
if (n.fiber[vox0]==0){
  good <- F
} else if (n.fiber[vox0]==1) {
  good <- acos(min(abs(eig[rmap[vox0],]%*%dir0),1))<thres.ang
} else {
  good <- as.logical(sum(as.vector(acos(pmin(abs(eig[rmap[vox0]+(0:(n.fiber[vox0]-1)),]%*%dir0),1)))<thres.ang))
}

#############################
## fiber.in.out
#############################
inc = path.in[ii,]-loc[path.iind[ii],]
direct=path.dir[ii,]


if (sum(dimens==0)){
  stop("directions has zero component, not yet supported! Please modify fiber.in.out\n")
}

# compute the distance of the current fiber directon to each face of the current voxel
tempdiff <- (cbind(dimens/2-inc,-inc-dimens/2)/direct) 
index1 <- which.min(diag(tempdiff[,2-(direct>=0)]))
index <- c(index1, (2-(direct>0))[index1])
const <- tempdiff[index[1],index[2]]
outc <- inc + const*direct
return(list(outc=outc, index=as.vector(index)))

#####

# compute the distance of the current fiber directon to each face of the current voxel
tempdiff <- (cbind(dimens/2-inc,-inc-dimens/2)/direct)

# Hao
tempdiff[tempdiff==Inf]=1e10
tempdiff[tempdiff==-Inf]=-1e10
tempdiff[is.nan(tempdiff)]=1e10
tempdiff[tempdiff==Inf]=1e10


# get which axis is the current fiber direction hitting face of the current voxel first
# 1:x  2:y  3:z
index1 <- which.min(diag(tempdiff[,2-(direct>=0)]))  # Hao change direct>0 to direct>=0
# which direction it is hitting 1:positive  2:negative
index <- c(index1, (2-(direct>0))[index1])
const <- tempdiff[index[1],index[2]]
outc <- inc + const*direct
return(list(outc=outc, index=as.vector(index)))


#############################
## get.out.vox
#############################
cvox = path.voxel[ii]

cvoxindex <- as.vector(arrayInd(cvox, braindim))
if (index[2]==1){
  # positive sides
  cvoxindex[index[1]] <- cvoxindex[index[1]] + vorient[index[1]]
} else {
  # negative sides
  cvoxindex[index[1]] <- cvoxindex[index[1]] - vorient[index[1]]
}
if ((cvoxindex[index[1]]<1)||(cvoxindex[index[1]]>braindim[index[1]])){
  return(NA)
} else {
  return(ArrayIndex(braindim, cvoxindex[1], cvoxindex[2], cvoxindex[3]))
}




#########

eig_min = rep(0,dim(eig)[1])
for(iind in 1:dim(eig)[1]){
  eig_min[iind] = min(abs(eig))
}
min(abs(eig),na.rm=T)