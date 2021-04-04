library(rTensor)
Objective = function (x, mu.array, Cs, Ds, Es) {
return(sum((x - mu.array[Cs, Ds, Es, drop=FALSE])^2,na.rm = T))
}

makediagNA = function(n){
    tnsr = array(1,dim = c(n,n,n))
    tnsr = array( apply( tnsr, 3, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
    tnsr = array( apply( tnsr, 2, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
    tnsr = array( apply( tnsr, 1, function(x) {x[ row(x) == col(x) ] <-NA; x} ), dim(tnsr) )
    return(tnsr)
}

tensor_unfold = function(tensor,dim=1){
    if (dim == 1) unfold = aperm(tensor,c(3,2,1))
    if (dim == 2) unfold = aperm(tensor,c(1,3,2))
    if (dim == 3) unfold = tensor
    unfold = apply(unfold,3,c)
    if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}


ReNumber = function (Cs){
    newCs <- rep(NA, length(Cs))
    uniq <- sort(unique(Cs))
    for (i in 1:length(uniq)) {
        newCs[Cs == uniq[i]] <- i
    }
    return(newCs)
}

## cut function
cut = function(tnsr){
    cond = (dim(tnsr)[1]==dim(tnsr)[2])&(dim(tnsr)[1]==dim(tnsr)[3])
    if (!(cond)){
        stop("* tensor must have same dimension for each mode")
    }
    tnsr = array( apply( tnsr, 3, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )
    tnsr = array( apply( tnsr, 2, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )
    tnsr = array( apply( tnsr, 1, function(x) {x[ row(x) == col(x) ] <-0; x} ), dim(tnsr) )

    return(tnsr)
}

UpdateMus_tensor = function (x, Cs, Ds, Es) {
    uniqCs = sort(unique(Cs))
    uniqDs = sort(unique(Ds))
    uniqEs = sort(unique(Es))
    
    mus = array(NA, c(length(uniqCs), length(uniqDs), length(uniqEs)))
        for (k in uniqCs){
            for (r in uniqDs){
                for (l in uniqEs){
                    mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l],na.rm = T)
                }
            }
        }
    return(mus)
}

UpdateClusters_tensor = function (x, mus, curCs, curDs) {
    Cs.new <- rep(NA, length(curCs))
    uniq <- 1:max(curCs)
    mus.expandcolumns <- mus[, curDs, drop = FALSE]
    for (i in 1:nrow(x)) {
        dist2.clust <- NULL
        for (k in 1:length(uniq)) {
            #see which cluster is closest to one sample
            dist2.clust <- c(dist2.clust, sum((x[i, , drop = FALSE] -
            mus.expandcolumns[k, , drop = FALSE])^2,na.rm = T))
        }
        wh <- which(dist2.clust == min(dist2.clust))
        Cs.new[i] <- wh[1]
    }
    return(Cs.new)
}


HSC = function(x,k,l,r){
    result = list()
    u1 = svd(tensor_unfold(x,1))$u[,1:k]
    u2 = svd(tensor_unfold(x,2))$u[,1:l]
    u3 = svd(tensor_unfold(x,3))$u[,1:r]
    hu1 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u2),t(u3)),ms = c(2,3))@data,1))$u[,1:k]
    hu2 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u3)),ms = c(1,3))@data,2))$u[,1:l]
    hu3 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u2)),ms = c(1,2))@data,3))$u[,1:r]
    Y1 = hu1%*%t(hu1)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu2),t(hu3)),ms = c(2,3))@data,1)
    Y2 = hu2%*%t(hu2)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu3)),ms = c(1,3))@data,2)
    Y3 = hu3%*%t(hu3)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu2)),ms = c(1,2))@data,3)
    result$Cs  = kmeans(Y1,k,nstart = 100)$cluster
    result$Ds  = kmeans(Y2,r,nstart = 100)$cluster
    result$Es  = kmeans(Y3,l,nstart = 100)$cluster
    
    return(result)
}

#Classify => tbmClustering
#Followings are changes
#1) I change the order of updates: from Cs=>Es=>Ds=>core to Cs=>core=>Es=>core=>Ds=>core
#This is because previous order of update makes the number of clustering in Es and Ds small

#2) Initialization really matters for the performance. I change nstart option in Kmeans from 25=>100 to make sure we have good initialization
#However, current one makes algorithm slow. I constructed HSC initialization but did not work well by now.
#I will figure out how to initialize group efficiently.

#3) Fix some bugs which happened when initial values (Cs.init,Ds.init,Es.init) are given.

tbmClustering = function(x,k,r,l,sym = F,diagP = T,max.iter=100,threshold = 1e-15,trace=FALSE,Cs.init=NULL,Ds.init=NULL,Es.init=NULL,nstart=100){
    n = dim(x)[1]; p = dim(x)[2]; q = dim(x)[3]
    if(is.null(Cs.init)){
        if(k==1) Cs = rep(1,n) else {Cs  = kmeans(tensor_unfold(x,1),k,nstart = nstart)$cluster}
    } else {
        Cs = ReNumber(Cs.init)
        k = length(unique(Cs))
    }
    if(is.null(Ds.init)){
        if(r==1) Ds = rep(1,p) else {Ds  = kmeans(tensor_unfold(x,2),r,nstart = nstart)$cluster}
    } else {
        Ds = ReNumber(Ds.init)
        r = length(unique(Ds))
    }
    if(is.null(Es.init)){
        if(l==1) Es = rep(1,q) else {Es  = kmeans(tensor_unfold(x,3),l,nstart = nstart)$cluster}
    } else {
        Es = ReNumber(Es.init)
        l = length(unique(Es))
    }
    if(diagP==F){
        x = x*makediagNA(n)
    }
    objs <- 1e+15
    improvement <- 1e+10
    i <- 1
    mu.array = UpdateMus_tensor(x,Cs,Ds,Es)
    while((improvement > threshold) & (i <= max.iter)){
        ### first mode
        Cs = UpdateClusters_tensor(tensor_unfold(x),tensor_unfold(mu.array),Cs,(rep(Ds,each=q)-1)*l+rep(Es,times=p))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        Cs <- ReNumber(Cs)
        k = length(unique(Cs))
        mu.array = UpdateMus_tensor(x,Cs,Ds,Es)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        ### second mode
        Ds = UpdateClusters_tensor(tensor_unfold(x,2),tensor_unfold(mu.array,2),Ds,(rep(Es,each=n)-1)*k+rep(Cs,times=q))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        Ds <- ReNumber(Ds)
        r = length(unique(Ds))
        mu.array = UpdateMus_tensor(x,Cs,Ds,Es)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        ### third mode
        Es = UpdateClusters_tensor(tensor_unfold(x,3),tensor_unfold(mu.array,3),Es,(rep(Ds,each=n)-1)*k+rep(Cs,times=p))
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        Es <- ReNumber(Es)
        l = length(unique(Es))
        mu.array = UpdateMus_tensor(x,Cs,Ds,Es)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
        ### end of update
        improvement <- abs(objs[length(objs)] - objs[length(objs) -
        6])/abs(objs[length(objs) - 6])
        i <- i + 1
        if(trace) cat("step",i,",improvement=",improvement,".\n")
        if (is.na(improvement)) break
    }
    if (sym ==T){
        symobj = NULL
        
        temp = UpdateMus_tensor(x,Cs,Cs,Cs)
        symobj = c(symobj,Objective(x, temp, Cs, Cs, Cs))
        temp = UpdateMus_tensor(x,Ds,Ds,Ds)
        symobj = c(symobj,Objective(x, temp, Ds, Ds, Ds))
        temp = UpdateMus_tensor(x,Es,Es,Es)
        symobj = c(symobj,Objective(x, temp, Es, Es, Es))
        Cs = Ds = Es = list(Cs,Ds,Es)[[which.min(symobj)]]
        mu.array = UpdateMus_tensor(x,Cs,Ds,Es)
        objs <- c(objs, Objective(x, mu.array, Cs, Ds, Es))
    }
    if (i > max.iter) {
        warning("The algorithm has not converged by the specified maximum number of iteration.\n")
    }
    if(diagP==F){
        judgeX = cut(mu.array[Cs,Ds,Es, drop=FALSE])
    }
    return(list("judgeX"=judgeX,"Cs"=Cs,"Ds"=Ds,"Es"=Es,"objs"=objs[length(objs)], "mus"=mu.array))
}

