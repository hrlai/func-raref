#---------------------------------------------------------------------------------------------
# R code for computing functional diversity in the paper 
# Chiu, C.-H. and Chao, A. (2014). Distance-based functional diversity measures 
# and their decomposition: a framework based on Hill numbers.PLoS ONE 
# 9(7):e100014.
#--------------------------------------------------------------------------------------------

Func2014=function (Dis,abun,q)  
  # input:
  # Dis: species pairwise functional distance matrix.
  # abun: one assemblage species abundance vector data or multi-assemblage
  # species by assemblages matrix data, where row is species and column is 
  # community or site(plot).
  # q: a non-negative value specifying diversity order.
  
  #output:
  # Q: Rao's quadratic entropy for each community.
  # FuncD = functional diversity of each community. 
  # Gamma = functional gamma diversity.
  # Alpha = functional alpha diversity.
  # Beta = functional beta diversity.
  # FunCqN = functional overlap index CqN(similarity index) 
  # FunUqN = functional overlap index UqN(similarity index)
{
  if(is.vector(abun)){
    Dis=as.matrix(Dis); 
    n=sum(abun);
    p=abun/n;I=which(p>0);p=p[I];
    Q=c(t(p)%*%Dis[I,I]%*%p);
    temp=p%*%t(p);
    if(q==1){FD= exp(sum(-Dis[I,I]*temp/Q*log(temp/Q)))}
    else{FD=(t(p^q)%*%Dis[I,I]%*%(p/Q)^q)^(1/(1-q));}
    #output=matrix(ncol=2,nrow=1);
    output=c(Q,FD)
    names(output)=c("Q","FuncD")
    return(output);
  }else{
    abun=as.matrix(abun);  
    N=ncol(abun);n=colSums(abun)
    FuncD=numeric(N);Q=numeric(N)
    Dis=as.matrix(Dis); 
    for(i in 1:N){            
      Q[i]=t(abun[,i]/n[i])%*%Dis%*%(abun[,i]/n[i]);
      temp=(abun[,i]/n[i])%*%t(abun[,i]/n[i]);
      I=which(temp>0);
      
      if(q==1){ FuncD[i]= exp(sum(-Dis[I]*temp[I]/Q[i]*log(temp[I]/Q[i]))); }
      else{ FuncD[i]=sum(Dis[I]*(temp[I]/Q[i])^q)^(1/(1-q)); }
    }
    
    gn=sum(abun);
    pop=abun/gn;  
    p=rowSums(pop);gI=which(p>0);p=p[gI];
    gQ=c(t(p)%*%Dis[gI,gI]%*%p);
    
    # return( list(Q=Q,FuncD=FuncD, output=gQ))               
    return( list(Q=Q,FuncD=FuncD))               
  }
}


###pairwise functional similarity indices for CqN and UqN similarity indices
# input:
# abun:species by community matrix dataframe, where row is species and
# column is community or site.
# q: order q should be any non-negative value.
# output: pairwise functional similarity matrix 
# CqN: Sorensen-type overlap index (from a local view)(similarity index) 
# UqN: Jaccard-type overlap index (from a regional view)(similarity index)

pairFunc=function(Dis, abun,q){
  N=ncol(abun);
  CqN=matrix(1,ncol=N,nrow=N);UqN=CqN;
  for(i in 1:(N-1)){
    for(j in i:N){
      o=Func2014(Dis,abun[,c(i,j)],q);
      CqN[i,j]=o[[3]][5];CqN[j,i]=CqN[i,j];
      UqN[i,j]=o[[3]][6];UqN[j,i]=UqN[i,j];
    }
  }
  return(list(CqN=CqN,UqN=UqN));
}






# Modified by Hao Ran Lai to only return Q --------------------------------
pairQ=function(Dis, abun, q=0){
  N=ncol(abun);
  Q=matrix(1,ncol=N,nrow=N,dimnames=list(colnames(abun),colnames(abun)))
  for(i in 1:N){
    for(j in i:N){
      o=Func2014(Dis,abun[,c(i,j)],q);
      Q[i,j]=o[[3]][1];Q[j,i]=Q[i,j];
    }
  }
  return(Q);
}
