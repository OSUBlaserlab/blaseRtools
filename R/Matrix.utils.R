# From https://github.com/cran/Matrix.utils/blob/master/R/Matrix.utils.R 
#' @import Matrix
#' @import grr
#' @importFrom stats aggregate as.formula terms contrasts na.omit na.pass
#' @importFrom methods is as
NULL

#' Casts or pivots a long \code{data frame} into a wide sparse matrix
#' 
#' Similar in function to \code{\link[reshape2]{dcast}}, but produces a sparse 
#' \code{\link{Matrix}} as an output. Sparse matrices are beneficial for this 
#' application because such outputs are often very wide and sparse. Conceptually
#' similar to a \code{pivot} operation.
#' 
#' Casting formulas are slightly different than those in \code{dcast} and follow
#' the conventions of \code{\link{model.matrix}}. See \code{\link{formula}} for 
#' details.  Briefly, the left hand side of the \code{~} will be used as the 
#' grouping criteria.  This can either be a single variable, or a group of 
#' variables linked using \code{:}.  The right hand side specifies what the 
#' columns will be. Unlike \code{dcast}, using the \code{+} operator will append
#' the values for each variable as additional columns.  This is useful for 
#' things such as one-hot encoding.  Using \code{:} will combine the columns as 
#' interactions.
#' 
#' @param data a data frame
#' @param formula casting \code{\link[stats]{formula}}, see details for specifics.
#' @param fun.aggregate name of aggregation function.  Defaults to 'sum'
#' @param value.var name of column that stores values to be aggregated numerics
#' @param as.factors if TRUE, treat all columns as factors, including
#' @param factor.nas if TRUE, treat factors with NAs as new levels.  Otherwise, 
#'  rows with NAs will receive zeroes in all columns for that factor
#' @param drop.unused.levels should factors have unused levels dropped? Defaults to TRUE, 
#'  in contrast to \code{\link{model.matrix}}
#' @return a sparse \code{Matrix}
#' @seealso \code{\link[reshape]{cast}}
#' @seealso \code{\link[reshape2]{dcast}}
dMcast<-function(data,formula,fun.aggregate='sum',value.var=NULL,as.factors=FALSE,factor.nas=TRUE,drop.unused.levels=TRUE)
{
  values<-1
  if(!is.null(value.var))
    values<-data[,value.var]
  alltms<-terms(formula,data=data)
  response<-rownames(attr(alltms,'factors'))[attr(alltms,'response')]
  tm<-attr(alltms,"term.labels")
  interactionsIndex<-grep(':',tm)
  interactions<-tm[interactionsIndex]
  simple<-setdiff(tm,interactions)
  i2<-strsplit(interactions,':')
  newterms<-unlist(lapply(i2,function (x) paste("paste(",paste(x,collapse=','),",","sep='_'",")")))
  newterms<-c(simple,newterms)
  newformula<-as.formula(paste('~0+',paste(newterms,collapse='+')))
  allvars<-all.vars(alltms)
  data<-data[,c(allvars),drop=FALSE]
  if(as.factors)
    data<-data.frame(lapply(data,as.factor))
  characters<-unlist(lapply(data,is.character))
  data[,characters]<-lapply(data[,characters,drop=FALSE],as.factor)
  factors<-unlist(lapply(data,is.factor))
  #Prevents errors with 1 or fewer distinct levels
  data[,factors]<-lapply(data[,factors,drop=FALSE],function (x) 
  {
    if(factor.nas)
      if(any(is.na(x)))
      {
        levels(x)<-c(levels(x),'NA')
        x[is.na(x)]<-'NA'
      }
    if(drop.unused.levels)
        if(nlevels(x)!=length(na.omit(unique(x))))
          x<-factor(as.character(x))
    y<-contrasts(x,contrasts=FALSE,sparse=TRUE)
    attr(x,'contrasts')<-y
    return(x)
  })
  #Allows NAs to pass
  attr(data,'na.action')<-na.pass
  result<-sparse.model.matrix(newformula,data,drop.unused.levels = FALSE,row.names=FALSE)
  brokenNames<-grep('paste(',colnames(result),fixed = TRUE)
  colnames(result)[brokenNames]<-lapply(colnames(result)[brokenNames],function (x) {
    x<-gsub('paste(',replacement='',x=x,fixed = TRUE) 
    x<-gsub(pattern=', ',replacement='_',x=x,fixed=TRUE) 
    x<-gsub(pattern='_sep = \"_\")',replacement='',x=x,fixed=TRUE)
    return(x)
  })

  result<-result*values
  if(isTRUE(response>0))
  {
    responses=all.vars(terms(as.formula(paste(response,'~0'))))
    result<-aggregate.Matrix(result,data[,responses,drop=FALSE],fun=fun.aggregate)
  }
  return(result)
}

#' Compute summary statistics of a Matrix
#' 
#' Similar to \code{\link[stats]{aggregate}}.  Splits the matrix into groups as 
#' specified by groupings, which can be one or more variables. Aggregation 
#' function will be applied to all columns in data, or as specified in formula. 
#' Warning: groupings will be made dense if it is sparse, though data will not.
#' 
#' \code{aggregate.Matrix} uses its own implementations of functions and should
#' be passed a string in the \code{fun} argument.
#' 
#' @param x a \code{\link{Matrix}} or matrix-like object
#' @param groupings an object coercible to a group of factors defining the
#'   groups
#' @param form \code{\link[stats]{formula}}
#' @param fun character string specifying the name of aggregation function to be
#'   applied to all columns in data.  Currently "sum", "count", and "mean"
#'   are supported.
#' @param ... arguments to be passed to or from methods.  Currently ignored
#' @return A sparse \code{Matrix}.  The rownames correspond to the values of the
#'   groupings or the interactions of groupings joined by a \code{_}.
#'   
#'   There is an attribute \code{crosswalk} that includes the groupings as a
#'   data frame.  This is necessary because it is not possible to include
#'   character or data frame groupings in a sparse Matrix.  If needed, one can 
#'   \code{cbind(attr(x,"crosswalk"),x)} to combine the groupings and the
#'   aggregates.
#'   
#' @seealso \code{\link[dplyr]{summarise}}
#' @seealso \code{\link[plyr]{summarise}}
#' @seealso \code{\link[stats]{aggregate}}
aggregate.Matrix<-function(x,groupings=NULL,form=NULL,fun='sum',...)
{
  if(!is(x,'Matrix'))
    x<-Matrix(as.matrix(x),sparse=TRUE)
  if(fun=='count')
    x<-x!=0
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-data.frame(interaction(groupings2,sep = '_'))
  colnames(groupings2)<-'A'
  if(is.null(form))
    form<-as.formula('~0+.')
  form<-as.formula(form)
  mapping<-dMcast(groupings2,form)
  colnames(mapping)<-substring(colnames(mapping),2)
  result<-t(mapping) %*% x
  if(fun=='mean')
    result@x<-result@x/(aggregate.Matrix(x,groupings2,fun='count'))@x
  attr(result,'crosswalk')<-grr::extract(groupings,match(rownames(result),groupings2$A))
  return(result)
}

aggregate2.Matrix<-function(x,groupings=NULL,form=NULL,fun=sum,...)
{
  #if(!is(x,'Matrix'))
    x<-as.matrix(x)
  groupings2<-groupings
  if(!is(groupings2,'data.frame'))
    groupings2<-as(groupings2,'data.frame')
  groupings2<-data.frame(lapply(groupings2,as.factor))
  groupings2<-interaction(groupings2,sep = '_')
  index<-grr::order2(groupings2)
  breaks<-which(!duplicated(groupings2[index]))
  results1<-fun(x[1:breaks[1],])
  results<-matrix(results1,ncol=length(results1),nrow=length(breaks))
  for (i in seq_len(length(breaks)-1)) 
  {
    results[i+1]<-fun(x[(breaks[i]+1):breaks[i+1],])
  }
  return(results)
}
  

#'Merges two Matrices or matrix-like objects
#'
#'Implementation of \code{\link{merge}} for \code{\link{Matrix}}.  By explicitly
#'calling \code{merge.Matrix} it will also work for \code{matrix}, for 
#'\code{data.frame}, and \code{vector} objects as a much faster alternative to 
#'the built-in \code{merge}.
#'
#'#' \code{all.x/all.y} correspond to the four types of database joins in the 
#'following way:
#'
#'\describe{ \item{left}{\code{all.x=TRUE}, \code{all.y=FALSE}} 
#'\item{right}{\code{all.x=FALSE}, \code{all.y=TRUE}} 
#'\item{inner}{\code{all.x=FALSE}, \code{all.y=FALSE}} 
#'\item{full}{\code{all.x=TRUE}, \code{all.y=TRUE}} }
#'
#'Note that \code{NA} values will match other \code{NA} values.
#'
#'@param x,y \code{Matrix} or matrix-like object
#'@param by.x vector indicating the names to match from \code{Matrix} x
#'@param by.y vector indicating the names to match from \code{Matrix} y
#'@param all.x logical; if \code{TRUE}, then each value in \code{x} will be 
#'  included even if it has no matching values in \code{y}
#'@param all.y logical; if \code{TRUE}, then each value in \code{y} will be 
#'  included even if it has no matching values in \code{x}
#'@param out.class the class of the output object.  Defaults to the class of x. 
#'  Note that some output classes are not possible due to R coercion
#'  capabilities, such as converting a character matrix to a Matrix.
#'@param fill.x,fill.y the value to put in merged columns where there is no match.
#'  Defaults to 0/FALSE for sparse matrices in order to preserve sparsity, NA for
#'  all other classes
#'@param ... arguments to be passed to or from methods.  Currently ignored
merge.Matrix<-function(x,y,by.x,by.y,all.x=TRUE,all.y=TRUE,out.class=class(x)[1],
                       fill.x=ifelse(is(x,'sparseMatrix'),FALSE,NA),fill.y=fill.x,...)
{
  requireNamespace('grr')
  if(is.null(dim(x)))
    return(grr::matches(by.x,by.y,all.x,all.y,indexes=FALSE))
  indices<-grr::matches(by.x,by.y,all.x,all.y,nomatch = NULL)
  x<-rbind(x,fill.x)
  x<-as(grr::extract(x,indices$x),out.class)
  
  y<-rbind(y,fill.y)
  if(!is.null(colnames(x)) & !is.null(colnames(y)))
    colnames(y)[colnames(y) %in% colnames(x)]<-paste('y',colnames(y)[colnames(y) %in% colnames(x)],sep='.')

  y<-as(grr::extract(y,indices$y),out.class)
  result<-cbind2(x,y)
  return(result)
}

#' @rdname merge.Matrix
join.Matrix<-merge.Matrix

#' Combine matrixes by row, fill in missing columns
#' 
#' rbinds a list of Matrix or matrix like objects, filling in missing columns.
#' 
#' Similar to \code{\link[plyr]{rbind.fill.matrix}}, but works for
#' \code{\link{Matrix}} as well as all other R objects.  It is completely
#' agnostic to class, and will produce an object of the class of the first input
#' (or of class \code{matrix} if the first object is one dimensional).
#' 
#' The implementation is recursive, so it can handle an arbitrary number of 
#' inputs, albeit inefficiently for large numbers of inputs.
#' 
#' This method is still experimental, but should work in most cases.  If the
#' data sets consist solely of data frames, \code{\link[plyr]{rbind.fill}} is 
#' preferred.
#' 
#' @param x,... Objects to combine.  If the first argument is a list and 
#'   \code{..} is unpopulated, the objects in that list will be combined.
#' @param fill value with which to fill unmatched columns
#' @param out.class the class of the output object.  Defaults to the class of x. 
#'  Note that some output classes are not possible due to R coercion
#'  capabilities, such as converting a character matrix to a Matrix.
#' @return a single object of the same class as the first input, or of class
#'   \code{matrix} if the first object is one dimensional
#' @seealso \code{\link[plyr]{rbind.fill}}
#' @seealso \code{\link[plyr]{rbind.fill.matrix}}
rBind.fill<-function(x,...,fill=NULL,out.class=class(rbind(x,x))[1])
{
  if (is.list(x) && !is.data.frame(x) && missing(...)) {
    Reduce(function (x,y) rBind.fill.internal(x,y,fill,out.class),x)
  }
  else {
    Reduce(function (x,y) rBind.fill.internal(x,y,fill,out.class),list(x,...))
  }
}

rBind.fill.internal<-function(x,y,fill,out.class)
{
  out.class<-force(out.class)
  fillMissing<-is.null(fill)
  if(fillMissing)
    fill<-if(is(x,'Matrix')) 0 else NA
  if (is.null(nrow(x)))
       x<-matrix(x,nrow=1,dimnames=list(NULL,names(x)))
  if (is.null(nrow(y)))
      y<-matrix(y,nrow=1,dimnames=list(NULL,names(y)))
  
  nullNames<-FALSE
  #Cannot currently handle duplicate column names
  if(is.null(colnames(x)))
    colnames(x)<-make.names(colnames(y)[1:ncol(x)],unique = TRUE)
  if(is.null(colnames(y)))
     colnames(y)<-make.names(colnames(x)[1:ncol(y)],unique = TRUE)
  if(is.null(colnames(x)))
  {
    nullNames<-TRUE
    colnames(x)<-1:ncol(x)
    colnames(y)<-1:ncol(y)
  }
  ymiss<-colnames(x)[which(is.na(match(colnames(x),colnames(y))))]
  ybind<-rsparsematrix(nrow=nrow(y),ncol=length(ymiss),0)
  colnames(ybind)<-ymiss
  if(!fillMissing)
    ybind[seq_along(ybind)]<-fill
  xmiss<-colnames(y)[which(is.na(match(colnames(y),colnames(x))))]
  xbind<-rsparsematrix(nrow=nrow(x),ncol=length(xmiss),0)
  colnames(xbind)<-xmiss
  if(!fillMissing)
    xbind[seq_along(xbind)]<-fill
  if (ncol(xbind>0))
    x<-cbind2(x,xbind)
  if(ncol(ybind)>0)
    y<-cbind2(y,ybind)
  y<-as(y,out.class)
  x<-as(x,out.class)
  result<-rbind2(x,y[,order(match(colnames(y),colnames(x)))])
  if(nullNames)
    colnames(result)<-NULL
  return(result)
}

len<-function (data) 
{
  result <- ifelse(is.null(nrow(data)), length(data), nrow(data))
  return(result)
}

setAs('Matrix','data.frame',function (from) as.data.frame(as.matrix(from)),)

setAs('Matrix','data.frame',function (from) as.data.frame(as.matrix(from)))

setAs('data.frame','dgeMatrix', function (from) as(as.matrix(from),'dgeMatrix'))

setAs('data.frame','dgCMatrix', function (from) as(as.matrix(from),'dgCMatrix'))

setAs('matrix','data.frame',function (from) as.data.frame(from))

setAs('vector','data.frame',function (from) data.frame(from))

setMethod("cbind2",c('data.frame','Matrix'),function (x,y) 
{
  y<-as.matrix(y)
  cbind2(x,y)
}
)

setMethod("cbind2",c('Matrix','data.frame'),function (x,y) 
{
  y<-as.matrix(y)
  cbind2(x,y)
}
)
