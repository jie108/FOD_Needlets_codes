rm(list=ls())
library(R.matlab)
library(xtable)
#### load library ####
library(rgl)
library(compositions)
source("dwi_fit.R")
source("dwi_track.R")


path_load = '/Users/hao/Dropbox/stats_project/FOD_codes_simulation/Real_data/S110933/fitting/space_indexx108-123y124-139z37-42/'
num_fib_cut = 4
x_subr = 1:16;
y_subr = 1:16;
z_subr = 4:4;
voxel_sub = paste0('x',toString(x_subr[1]),'-',toString(x_subr[length(x_subr)]),'y',toString(y_subr[1]),'-',toString(y_subr[length(y_subr)]),'z',toString(z_subr[1]),'-',toString(z_subr[length(z_subr)]));

temp = readMat(paste0(path_load,'for_tracking_cut',toString(num_fib_cut),'_',voxel_sub,'.mat'))
temp$n.fiber2 = c(temp$n.fiber2)



our.track <- v.track(v.obj=temp, xgrid.sp=1.3672, ygrid.sp=1.3672,
                     zgrid.sp=2.7, braingrid=array(temp$braingrid,dim=c(3,length(x_subr),length(y_subr),length(z_subr))), elim=T, nproj=2,
                     vorient=c(1,1,1), elim.thres=5.5)


############################
#### fiber realizations ####
############################
tobj <- our.track

length(tobj$sorted.iinds[tobj$sorted.update.ind])
ndis <- length(tobj$sorted.iinds[tobj$sorted.update.ind])  #length(tobj$sorted.iinds[tobj$sorted.update.ind]) # number of fibers

open3d()
for (iind in (tobj$sorted.iinds[tobj$sorted.update.ind])[1:ndis]){
  cat(iind,"\n")
  # plot
  plot.fib(tobj$tracks1[[iind]]$inloc, tobj$tracks1[[iind]]$dir)
  plot.fib(tobj$tracks2[[iind]]$inloc, tobj$tracks2[[iind]]$dir)
}
decorate3d(xlim=range(temp$braingrid[1,,,]), ylim=range(temp$braingrid[2,,,]),
           zlim=range(temp$braingrid[3,,,]), box=F)
title3d(main=paste0('S110933_cut',toString(num_fib_cut),'_numfib',toString(ndis)),pos=c(15,15,15))
