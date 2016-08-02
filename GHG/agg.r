#index file has NA cells for blanks due to uneven length of each industry index. If the blanks are filled with 0, then use agg2.r

agg=function(indr,indc,data){
# input data indr, indc are data.frame, data should be matrix to apply colSums, rowSums
	K=dim(indr)[1]
	L=dim(indc)[1] 	
# Preparing storage
	stack_r=matrix(nrow=K,ncol=dim(data)[2])
	stack_c=matrix(nrow=K,ncol=L)
#row integration
      for (i in 1:K){
#			print(c("i=",i),quote=F)
# choose index for industry/section i
			indr_i=as.vector(na.omit(as.matrix(indr)[i,]))
#			print(c("indr_i=",indr_i),quote=F)
# if industry i composed of multi industry, take sum over column(vertical sum) and store in stack_r[i,]
			if (length(indr_i)>1){sum_i=colSums(data[indr_i,]) 
						    } 
#if industry i composed of single industry, take the data and store in stack_r[i,]
			if (length(indr_i)==1){sum_i=data[indr_i,]
						    }
#			print(c("sum_i=",sum_i),quote=F)
			stack_r[i,]=sum_i}

#column integrateion using row-integrated data
	for (j in 1:L){
#			print(c("j=",j),quote=F)
# choose index for industry/section j
			indc_j=as.vector(na.omit(as.matrix(indc)[j,]))
#			print(c("indc_j=",indc_j),quote=F)
#We are going to use stack_r, row sum data.
#For multi industry, take sum over column(horizontal sum) of the row sum and store in stack_c[,.j]
			if (length(indc_j)>1){sum_j=rowSums(stack_r[,indc_j])
						   }
#For single industry/sector, put row sum data in stack_c[,.j]
			if (length(indc_j)==1){sum_j=stack_r[,indc_j]
						   }		
#			print(c("sum_j=",sum_j),quote=F)
			stack_c[,j]=sum_j
			}
#report the integrated rejult
	return(stack_c)
}

#[test drive]
#indr=read.csv(file="indr.csv",header=F)
#indc=read.csv(file="indc.csv",header=F)
#xdata=read.csv(file="xdata.csv",header=F)
#xdata.m=as.matrix(xdata)

#data1=agg(indr,indc,xdata.m)
#save(data1,file="data1.Rdata")


#
agg2=function(indr,indc,data){
# input data indr, indc are data.frame, data should be matrix to apply colSums, rowSums
	K=dim(indr)[1]
	L=dim(indc)[1] 	
# Preparing storage
	stack_r=matrix(nrow=K,ncol=dim(data)[2])
	stack_c=matrix(nrow=K,ncol=L)
#row integration
	for (i in 1:K){
#			print(c("i=",i),quote=F)
# choose index for industry/section i
		indr_i=indr[i,]
		indr_i=indr_i[indr_i>0]
#			print(c("indr_i=",indr_i),quote=F)
# if industry i composed of multi industry, take sum over column(vertical sum) and store in stack_r[i,]
		if (length(indr_i)>1){sum_i=colSums(data[indr_i,]) 
		} 
#if industry i composed of single industry, take the data and store in stack_r[i,]
		if (length(indr_i)==1){sum_i=data[indr_i,]
		}
#			print(c("sum_i=",sum_i),quote=F)
		stack_r[i,]=sum_i}
	
#column integrateion using row-integrated data
	for (j in 1:L){
#			print(c("j=",j),quote=F)
# choose index for industry/section j
#			indc_j=as.vector(na.omit(as.matrix(indc)[j,]))
		indc_j=indc[j,]
		indc_j=indc_j[indc_j>0]
#			print(c("indc_j=",indc_j),quote=F)
#We are going to use stack_r, row sum data.
#For multi industry, take sum over column(horizontal sum) of the row sum and store in stack_c[,.j]
		if (length(indc_j)>1){sum_j=rowSums(stack_r[,indc_j])
		}
#For single industry/sector, put row sum data in stack_c[,.j]
		if (length(indc_j)==1){sum_j=stack_r[,indc_j]
		}		
#			print(c("sum_j=",sum_j),quote=F)
		stack_c[,j]=sum_j
	}
#report the integrated rejult
	return(stack_c)
}

  
