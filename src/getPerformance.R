
library(boot)
library(ROCR)

crossValidate <- function(x, glm0, nsim=100, k=10)
{
	p0.vect<-vector()
	for (j in 1:nsim)
	{
   		cv0<-cv.glm(data.frame(x), glm0, K=k)
   		p0<-cv0$delta[2]
		p0.vect<-c(p0.vect,p0)
	}
	return(p0.vect)
}

plotROC <- function(glm0, y)
{
	pred <- prediction(glm0$fitted.values, y)
	perf <- performance(pred, "tpr", "fpr")
	plot(perf, colorize=TRUE)
	abline(0,1)
}

getAUC <- function(glm0,y)
{
	pred <- prediction(glm0$fitted.values, y)
	return(performance(pred, 'auc'))
}
