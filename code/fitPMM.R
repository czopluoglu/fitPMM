
# IF YOU ARE NOT INTERESTED IN THE DETAILS OF THE ESTIMATION PROCEDURE AND R ROUTINE, AND YOUR GOAL IS ONLY TO ANALYZE
# YOUR DATA, THEN PLEASE CLOSE THIS FILE AND OPEN THE FILE NAMED "demo_analysis.R".

# YOU HAVE TO SAVE THIS R SCRIPT FILE IN A FOLDER IN YOUR COMPUTER, AND KEEP THE FOLDER PATH IN YOUR MIND
# AS YOU WILL NEED IT FOR THE "demo_analysis.R" SCRIPT.

# NOW, YOU CAN CLOSE THIS FILE AND GO TO THE "demo_analysis.R" SCRIPT.

# "demo_analysis.R" SCRIPT HAS A SAMPLE DATASET AND PROVIDE THE R CODE TO ANALYZE THE SAMPLE DATASET USING
# THE ROUTINE IN THIS R SCRIPT.

# IF YOU HAVE QUESTIONS, PLEASE CONTACT ME:
#	CENGIZ ZOPLUOGLU
#	E-MAIL ADDRESS: c.zopluoglu@miami.edu

############################################################################################################################
############################################################################################################################
############################################################################################################################

# Function to fit finite mixtures of piecewise mixed effects model with unknown random knots


fitPMM  <- function(data,
                    start,
			  time.vec,
			  L,
			  nquad     =10,
                    isnit_max =15,
                    is_st     =0.01,
			  is_nset   =20,
			  fsnit_max =1000,
			  fs_st     =0.001,
			  fs_nset   =3,
			  sel       =2,
			  fix.cov=TRUE,
                    fix.err = TRUE
			  )
{

	
#######################################################################################################################
# INPUT:
#	   data     , N x m matrix, N is the number of observations and m is the number of timepoints
#        start    , N x 5 matrix which stores the individual estimates from the piecewise model
#                   This is used as a core set to generate multiple sets of start values by using random perturbations
#        time.vec , a vector of numbers with length "m" to specify the time metric to be used in the analysis
#	   L        , number of classes in the model
#        nquad    , number of quadrature points to approximate the likelihood value
#        isnit_max, maximum number of iterations for the initial stage estimation
#        is_st    , stopping criteria (difference in two successive log-likelihoods) for the inital stage      
#        is_nset  , number of start value sets for the initial stage estimation
#        fsnit_max, maximum number of iterations for the final stage estimations
#        fs_st    , stopping criteria (difference in two successive log-likelihoods) for the final stage
#        fs_nset  , number of start value sets with the best log-likelihood to be run at the final stage
#        sel      , Choose either BFGS(1) or DFPQ(2) method to update the parameters
#        fix.cov  , fix the variance covariance matrix of random parameters across classes
#        fix.err  , fix the residual variance across classes
#
# RETURN:
#	  
#        PAR	, parameter estimates from the best start value set
#        CP	      , class proportion estimates from the best start value set
#        POSTERIOR, posterior probabilities for the latent class membership
######################################################################################################################

cat("***********************************************************************","\n")
cat("An R Routine to Fit Finite Mixture of Piecewise Mixed-Effect Models with Unknown Random Knots","\n")
cat("","\n")
cat("Maintainer: Cengiz Zopluoglu","\n")
cat("            Research, Measurement, and Evaluation Program","\n")
cat("            Department of Educational and Psychological Studies","\n")
cat("            University of Miami","\n")
cat("","\n")
cat("Email     : c.zopluoglu@miami.edu","\n")
cat("*************************************************************************","\n")
cat("Please use the following citation for any use of this R routine.","\n")
cat("","\n")
cat("     Zopluoglu, C., Harring, J.R., Kohli,N. (2013). An R Routine to fit Piecewise","\n") 
cat("         Mixed-effects Mixture Models with Random Knots. Applied Psychological","\n") 
cat("         Measurement,??,",sprintf("%7s","??_??."),"\n") 
cat("","\n")
cat("*************************************************************************","\n")
cat("","\n")
cat("Processing Date: ",date(),"\n")
cat(" ","\n")

if(L==1) { 
	 fix.cov=FALSE
 	 fix.err = FALSE
}
######################################################################################################################
#
#                                 GENERATE THE STARTING PARAMETER SETS FOR THE INITIAL STAGE ITERATIONS
#
######################################################################################################################


# We follow the ideas in the following paper when generating multiple start value sets
#	 Hipp, J.R. & Bauer, D.J.(2006). Local Solutions in the Estimation of Growth Mixture Models.
#		Psychological Methods, 11(1), 36-53.


	# A start value set is a vector in the following order for one class:

		# B1,B2,B3,K,sigma,B1.var,B2.var,B3.var,K.var,
            # B1B2,B1B3,B1K,B2B3,B2K,B3K
		
			# B1     : Average intercept for the first phase               (Fixed effect)
			# B2     : Average slope for the first phase                   (Fixed effect)  
			# B3     : Average slope for the second phase                  (Fixed effect)
			# K      : Average Knot, change/break point between two phases (Fixed effect)
			# sigma  : Residual Variance (Fixed effect)

			# B1var : Variance of random intercepts 
			# B2var : Variance of random slopes for the first phase
			# B3var : Variance of random slopes for the second phase
			# Kvar  : Variance of random knots

			# B1B2   : Covariance between random B1s and B2s
			# B1B3   : Covariance between random B1s and B3s
			# B1K    : Covariance between random B1s and Ks
			# B2B3   : Covariance between random B2s and B3s
			# B2K    : Covariance between random B2s and Ks
			# B3K    : Covariance between random B3s and Ks


	start <- na.omit(start)  # Remove the individuals in which the model was not able to fitted successfullly
	start <- start[is.finite(rowSums(start)), ]
	par <- colMeans(start,na.rm=TRUE)
	ranef.cov <- cov(start[,1:4],use="pairwise.complete.obs")
      colnames(ranef.cov) <- c("B1.var","B2.var","B3.var","K.var")
	rownames(ranef.cov) <- c("B1.var","B2.var","B3.var","K.var")
	par <- c(par,diag(ranef.cov))
	par <- c(par,ranef.cov[2,1],ranef.cov[3,1],ranef.cov[4,1],
                   ranef.cov[3,2],ranef.cov[4,2],ranef.cov[4,3])
      names(par) <- c("B1","B2","B3","K","Sigma","B1var","B2var","B3var","Kvar",
                      "B1B2","B1B3","B1K", "B2B3","B2K","B3K")

	# Below, we generate a number of start value sets by using random perturbations of the "core" set.
	# See Hipp & Bauer paper for details.

	    B1.lower <- as.numeric(strsplit(summary(start)[2,1],":")[[1]][2])     
	    B1.upper <- as.numeric(strsplit(summary(start)[5,1],":")[[1]][2])       
	    B2.lower <- as.numeric(strsplit(summary(start)[2,2],":")[[1]][2])
	    B2.upper <- as.numeric(strsplit(summary(start)[5,2],":")[[1]][2])
	    B3.lower <- as.numeric(strsplit(summary(start)[2,3],":")[[1]][2])
	    B3.upper <- as.numeric(strsplit(summary(start)[5,3],":")[[1]][2])
	    K.lower  <- as.numeric(strsplit(summary(start)[2,4],":")[[1]][2])
	    K.upper  <- as.numeric(strsplit(summary(start)[5,4],":")[[1]][2])

	    B1var.lower <- .001
	    B2var.lower <- .001
	    B3var.lower <- .001
	    Kvar.lower  <- .001

	    B1var.upper <- var(start[,1],na.rm=TRUE)   
	    B2var.upper <- var(start[,2],na.rm=TRUE)
	    B3var.upper <- var(start[,3],na.rm=TRUE)
	    Kvar.upper  <- var(start[,4],na.rm=TRUE)
	
	    sigma <- mean(start[,5],na.rm=TRUE)
 
	       boundaries <- matrix(c(B1.lower,B1.upper,
        		                  B2.lower,B2.upper,
		                        B3.lower,B3.upper,
	 	                        K.lower,K.upper,
		                        sigma,sigma,
		                        B1var.lower,B1var.upper,
		                        B2var.lower,B2var.upper,
		                        B3var.lower,B3var.upper,
		                        Kvar.lower,Kvar.upper),
		                        nrow=9,ncol=2,byrow=TRUE)

		 rownames(boundaries) <- c("B1","B2","B3","K","sigma","B1var","B2var","B3var","Kvar")
	       colnames(boundaries) <- c("LowerBound","UpperBound")

		start.param.sets <- vector("list",L)

		for(cc in 1:L) {    # Loop over the number of classes
					  # Generate the sets of starting values for each class

	       start.param.sets[[cc]] <- matrix(nrow=is_nset,ncol=15)   

	       colnames(start.param.sets[[cc]]) <- c("B1","B2","B3","K","sigma","B1var","B2var","B3var","Kvar",
                                             "B1B2","B1B3","B1K","B2B3","B2K","B3K")

			for(j in 1:is_nset) {
		        for(i in 1:9) { 
		      	start.param.sets[[cc]][j,i]=runif(1,boundaries[i,1],boundaries[i,2])
			   }
 		      }

		      for(j in 1:is_nset) {

			     	# Boundaries for covariances [ COV = COR*(sqrt(VAR1*VAR2))  ]
				# assumes that upper and lower limits for the correlation among random effects is 0.75 and -0.75

				 B1B2.lower <- -.75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,7])
				 B1B2.upper <-  .75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,7])
				 B1B3.lower <- -.75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,8])
				 B1B3.upper <-  .75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,8])
				 B1K.lower  <- -.75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,9])
				 B1K.upper  <-  .75*sqrt(start.param.sets[[cc]][j,6]*start.param.sets[[cc]][j,9])
			       B2B3.lower <- -.75*sqrt(start.param.sets[[cc]][j,7]*start.param.sets[[cc]][j,8])
				 B2B3.upper <-  .75*sqrt(start.param.sets[[cc]][j,7]*start.param.sets[[cc]][j,8])
			       B2K.lower  <- -.75*sqrt(start.param.sets[[cc]][j,7]*start.param.sets[[cc]][j,9])
				 B2K.upper  <-  .75*sqrt(start.param.sets[[cc]][j,7]*start.param.sets[[cc]][j,9])
			       B3K.lower  <- -.75*sqrt(start.param.sets[[cc]][j,8]*start.param.sets[[cc]][j,9])
				 B3K.upper  <-  .75*sqrt(start.param.sets[[cc]][j,8]*start.param.sets[[cc]][j,9])

					boundaries1 <- matrix(c(B1B2.lower,B1B2.upper,
                        				      B1B3.lower,B1B3.upper,
			                                    B1K.lower,B1K.upper,
				                              B2B3.lower,B2B3.upper,
				                              B2K.lower,B2K.upper,
				                              B3K.lower,B3K.upper),nrow=6,ncol=2,byrow=TRUE)

				      rownames(boundaries1) <- c("B1B2","B1B3","B1K","B2B3","B2K","B3K")
				      colnames(boundaries1) <- c("LowerBound","UpperBound")
    

				 for(i in 10:15) { 
			      	start.param.sets[[cc]][j,i]=runif(1,boundaries1[i-9,1],boundaries1[i-9,2])
				 }
			}
    		}

# Combine the starting parameter sets over all classes in one vector

	if(fix.cov==FALSE & fix.err == FALSE) {

		start.parameters <- matrix(nrow=is_nset,ncol=15*L) 

		for(i in 1:is_nset) {
			comb <- start.param.sets[[1]][i,]
			if(L>=2) {
				for(j in 2:L) { comb <- c(comb,start.param.sets[[j]][i,]) }	
			}
			start.parameters[i,] = comb
		}

		colnames(start.parameters) <- names(comb)
	}

	if(fix.cov==TRUE & fix.err == FALSE) {

		start.parameters <- matrix(nrow=is_nset,ncol=15*L) 

		for(i in 1:is_nset) {
			comb <- start.param.sets[[1]][i,]
			if(L>=2) {
				for(j in 2:L) { comb <- c(comb,start.param.sets[[j]][i,]) }	
			}
			start.parameters[i,] = comb
		}

		colnames(start.parameters) <- names(comb)

		if(L>=2) { 
            	seeq = 1:15
			for(u in 2:L) { seeq=c(seeq,(15*(u-1)+1):(15*u-10)) }
      	}

		start.parameters <- start.parameters[,seeq]

	}

	if(fix.cov==FALSE & fix.err == TRUE) {

		start.parameters <- matrix(nrow=is_nset,ncol=15*L) 

		for(i in 1:is_nset) {
			comb <- start.param.sets[[1]][i,]
			if(L>=2) {
				for(j in 2:L) { comb <- c(comb,start.param.sets[[j]][i,]) }	
			}
			start.parameters[i,] = comb
		}

		colnames(start.parameters) <- names(comb)

		if(L>=2) { 
            	seeq = 1:15
			for(u in 2:L) { seeq=c(seeq,(((15*(u-1)+1):(15*u))[-5])) }
      	}

		start.parameters <- start.parameters[,seeq]

	}

	if(fix.cov==TRUE & fix.err == TRUE) {

		start.parameters <- matrix(nrow=is_nset,ncol=15*L) 

		for(i in 1:is_nset) {
			comb <- start.param.sets[[1]][i,]
			if(L>=2) {
				for(j in 2:L) { comb <- c(comb,start.param.sets[[j]][i,]) }	
			}
			start.parameters[i,] = comb
		}

		colnames(start.parameters) <- names(comb)

		if(L>=2) { 
            	seeq = 1:15
			for(u in 2:L) { seeq=c(seeq,(15*(u-1)+1):(15*u-11)) }
      	}

		start.parameters <- start.parameters[,seeq]

	}



cat(" ","\n")
cat(paste("Multiple sets of start values are written in a file 'StartValueSets.txt' in the following folder:",sep=""),"\n")
cat(getwd(),"\n")
cat(" ","\n")
cat("*************************************************************************","\n")
write.fwf(x=as.data.frame(start.parameters),
            file="StartValueSets.txt",
            rownames=TRUE,colnames=TRUE,sep="   ")


pik <- rep(1/L,L)
pik= pik + runif(L,-.03,+.03)
pik <- pik*(1/sum(pik))

###############################################################################
#
#         LIKELIHOOD FUNCTION for the GIVEN DATA and STARTING PARAMETERS 
#
##################################################################################


likelihood <- function(data,start.param,nquad=nquad,L=L,pis) {

	# data         is the N x t dataset

      # start.param  is a vector of starting parameters with a length Lx15

        # where L is the number latent classes
        # For instance, it's a vector of 30 elements when we want to fit a model
        # with two classes, First 15 are the parameters for the first, 
        # second 15 are the parameters for the second class

     # pis		   is starting values for the class proportions

      
#############################################################################################
	ind.likelihood <- function(par,data,nquad=nquad) {    # START internal function to compute individual 
                                                            # likelihoods for one class for a given 
                                                            # set of parameters and dataset

	  	# "par"     is a vector of starting values for 15 parameters
	      # "data"    is the dataset
	      # "nquad"   is the number of quadrature points

		# The computation of the loglikelihood function is based on three papers. Please see them for details
		#
		# 	du Toit, S.H.C. & Cudeck, R.(2009). Estimation of the Nonlinear Random Coefficient Model
		#		When Some Random Effects are Separable. Psychometrika, 74(1), 65-82.
		# 
		#	Harring, J.(2009). A Nonlinear Mixed Effects Model for Latent Variables. Journal of Educational
		#		and Behavioral Statistics, 34(3), 293-318.
		#
		#	Harring, J.(2012). Finite Mixtures of Nonlinear Mixed Effects Models. Book Chapter in 
		#		Advances in Longitudinal Methods in the Social and Behavioral Sciences by Harring & Hancock.

      
		
      	# Create Gauss-Hermit Quadrature Points and Corresponding Weights #

  		  abx <- as.matrix(ghq(nquad)$zeros)                                        
  		  wts <- as.matrix(ghq(nquad,modified=FALSE)$weights)
              wts <- wts/sqrt(pi)           

 		 # Vector for the time points

    		  x <- as.matrix(time.vec)   
	
 		 # Define model parameters

		   g1 = par[4]                            # Average knot 
	
		   thdelta <- par[5]*diag(length(time.vec))  # Error var-cov matrix
									# Independent error structure with constant error variance on the diagonal

 	         phibb = xpnd(par[c(6,10,11,7,13,8)])    # "xpnd" is a function to create matrix using lower diagonal elements; 
                                                       # phibb is the var-cov matrix among the linear random parameters 
									 # B0, B1, and B2 

		   phibg = rbind(par[12],par[14],par[15]) # phibg is the covariance matrix among the linear and nonlinear random parameters 

		   phigg = par[9]                         # phigg is the var-cov matrix among the nonlinear random parameters 

 	         Beta <- rbind(par[1],par[2],par[3])    # Average values for the intercept and slopes (B0, B1, B2)

   
		    # Compute likelihood values given a quadrature point
		    # Do a loop over quadrature points at the end
		    # The loop over people can be vectorized if there is no missing data for any individual. 
		    # This saves a lot of time. If missing data occurs, the loop is necessary.


		    lk <- matrix(nrow=nrow(data),ncol=nquad)
	
		    for(n in 1:nquad) {  
            
 		        u  = abx[n,1]
		        wk = wts[n,1]

 	     			gi = sqrt(2*phigg)*u                               # transformation of quadrature point following Harring(2009) paper -- Eq.(10)

		        mubg   = Beta  + (phibg%*%solve(phigg)%*%gi)          # Equations are on page 72, just below Eq.(6) in Harring(2009);
                                                                          # "solve()" function computes the inverse of the matrix

	 	        phibgi = phibb - (phibg%*%solve(phigg)%*%t(phibg))    # Equations are on page 72, just below Eq.(6) in Harring(2009)


		 	  # Design matrix for linear spline model, page 70 (du Toit & Cudeck ,2009)

		 	     xdes <- matrix(nrow=ncol(data),ncol=3)
			     xdes[,1]=1
			     xdes[,2]=time.vec
   			     xdes[,3]=time.vec
	      	     xdes[,2]=ifelse(xdes[,2]<=(g1+gi),xdes[,2],(g1+gi))
			     xdes[,3]=ifelse(xdes[,3]>(g1+gi),xdes[,3]-(g1+gi),0)


	      	    muyg  = xdes%*%mubg 			               # Eq.(7) on page 72 in Harring(2009)
		          phiyg = xdes%*%phibgi%*%t(xdes) + thdelta            # Eq.(7) on page 72 in Harring(2009)

			    muyg.matrix <- matrix(as.vector(muyg),nrow=nrow(t(data)),ncol=ncol(t(data)),byrow=FALSE) 

			 # If there is missing data, then use the following loop

				if(sum(is.na(data))!=0) {
	
			         qk <- c()
			   	    for(uuu in 1:nrow(data)) { 
					cases <- which(is.na(data[uuu,])==FALSE)
					qk[uuu] <- (data[uuu,cases]-muyg.matrix[cases,uuu])%*%solve(phiyg[cases,cases])%*%as.matrix(data[uuu,cases]-muyg.matrix[cases,uuu])
				    }

				   cst <- c()
			 	    for(uuu in 1:nrow(data)) { 
					cases <- which(is.na(data[uuu,])==FALSE)
					cst[uuu]=wk*(det(phiyg[cases,cases])^(-.5))
				    }
				}

			 # If there is no missing data, then use the following vectorization to make things faster

				if(sum(is.na(data))==0) {
				   qk <- diag(t(t(data)-muyg.matrix)%*%solve(phiyg)%*%(t(data)-muyg.matrix)) 
		          	   cst = wk*(det(phiyg)^(-.5))
				}

       	  contrib <- cst*exp(-.5*(qk))  
		  lk[,n]=contrib
           
	       }
	
     		 const = (2*pi)^(-(ncol(data)-rowSums(is.na(data)))/2)
		 lik <- const*rowSums(lk)                      # individual likelihoods
    	       return(lik)                                   # returns the individual likelihoods

	}  # end of the internal function for individual likelihoods given the item parameters and data

#########################################################################################################

# Create vectors for starting parameters given the "start.param"

	if(fix.cov==FALSE & fix.err == FALSE) {

	   	locate <- seq(from=1,to=length(start.param),by=15)      # positions for the start values for each class
      	params <- vector("list",L)
	      for(i in 1:length(locate)) { params[[i]] <- start.param[locate[i]:(locate[i]+14)] }
	}

	if(fix.cov==TRUE & fix.err == FALSE) {

	   	params <- vector("list",L)    
     	      params[[1]] <- start.param[1:15]
	
		if(L>=2) {
      	  locate <- seq(from=16,to=length(start.param),by=5)
	        for(u in 2:L) {
      	   params[[u]] = params[[1]]  
	         params[[u]][1:5] <- start.param[(locate[u-1]):(locate[u-1]+4)]
  		  }
		}
	}

	if(fix.cov==FALSE & fix.err == TRUE) {

	   	params <- vector("list",L)    
     	      params[[1]] <- start.param[1:15]
	
		if(L>=2) {
      	  locate <- seq(from=16,to=length(start.param),by=14)
	        for(u in 2:L) {
      	   params[[u]] = params[[1]]  
	         params[[u]][1:4] <- start.param[(locate[u-1]):(locate[u-1]+3)]
		   params[[u]][6:15] <- start.param[(locate[u-1]+4):(locate[u-1]+13)]
  		  }
		}
	}

	if(fix.cov==TRUE & fix.err == TRUE) {

	   	params <- vector("list",L)    
     	      params[[1]] <- start.param[1:15]
	
		if(L>=2) {
      	  locate <- seq(from=16,to=length(start.param),by=4)
	        for(u in 2:L) {
      	   params[[u]] = params[[1]]  
	         params[[u]][1:4] <- start.param[(locate[u-1]):(locate[u-1]+3)]
  		  }
		}
	}


# Compute individual likelihoods for separate sets of parameters (depending on the number of latent classes)

	ind.lik <- vector("list",L)

      for(i in 1:L) { 
		ind.lik[[i]] <- ind.likelihood(par=params[[i]],data=data,nquad=nquad)  
	}

# Compute the posterior probabilities for each individual (appendix, on page 187 in Harring (2012))
# Require this part of R code to update the pis for the next iteration.

	ind.prob <- matrix(nrow=nrow(data),ncol=L)
	for(i in 1:L) { ind.prob[,i] = pis[i]*ind.lik[[i]] }

	post.prob <- matrix(nrow=nrow(data),ncol=L)
	for(i in 1:L) { post.prob[,i] = ind.prob[,i]/rowSums(ind.prob) }

# Compute the complete data loglikelihood ---- Equation 7.8 on page 168 in Harring(2012)

	LL <- sum(log(rowSums(ind.prob)))

return(list(LL=LL,post.prob=post.prob))

} # END OF FUNCTION to COMPUTE LIKELIHOOD VALUE GIVEN the DATA and STARTING PARAMETERS 

##############################################################################################################

####### FUNCTION to compute the first derivatives of LL with respect to each parameter - Gradient  #######

# Gradient vector is computed numerically 
                 
grad <- function(data,start.param,nquad=nquad,L=2,pis) {  # start function

	lnull <- likelihood(data,start.param,nquad=nquad,L=L,pis=pis)$LL
	
	s=1e-5   						#small increment to compute first derivative numerically
      g <- c()
 
	for(j in 1:length(start.param)) {   
	
	  hold=start.param[j]
	  h=s*hold
        start.param[j]=hold+h; 
        lj=likelihood(data,start.param,nquad=nquad,L=L,pis=pis)$LL
        g[j]<-(lj - lnull)/h
        start.param[j]=hold
      }

 return(g)

} # end function 'grad'

###################### FUNCTION to compute the the second derivatives of LL  - HESSIAN  #######################

# Hessian is computed numerically 

nhess <- function(data,start.param,nquad=nquad,L=2,pis) { #start function

	g0 = grad(data,start.param,nquad=nquad,L=L,pis=pis)
	fh <- matrix(nrow=length(start.param),ncol=length(start.param))
	s=1e-5

	for(j in 1:length(start.param)) {

	  hold=start.param[j]
        h=s*hold
        start.param[j]=hold+h
        gj= grad(data,start.param,nquad=nquad,L=L,pis=pis)
        fh[,j]=(gj - g0)/ h
        start.param[j]= hold
	}

	hh = (fh + t(fh))/2 
      r = diag(fh)
	diag(hh) <- r

 return(hh)

} # end function
                 

##################        START INITIAL STAGE ITERATIONS            #######################

RESULTS.INITIAL <- vector("list",is_nset)    # List object to store results for each set of start values

for(set in 1:is_nset) {

if (set==1) { 
cat("","\n")
cat("                              INITIAL STAGE ESTIMATION              ","\n") 
cat("","\n")
cat("The results from the initial stage estimation are saved in the following folder:","\n") 
cat(getwd(),"\n")
cat("","\n")
}
cat("","\n")
cat(" START VALUE SET ",set,":","\n")
cat("","\n")


	####   First Iteration ##########
            
	l <- c()                                    # object to store the loglikelihood values for the data for each iteration
	oldpar <- vector("list",isnit_max)          # list object to store the estimates parameters at each iteration
	g <- vector("list",isnit_max)               # object to store the gradient vector at each iteration  
	H <- vector("list",isnit_max)               # object to store the Hessian matrix at each iteration  
	IH <- vector("list",isnit_max)              # object to store the inverse of the Hessian matrix at each iteration  
	maxg <- c()                                 # object to store the maximum gradient value at each iteration
	posterior <- vector("list",isnit_max)       # object to store individual posterior probabilities at each iteration
	class.prop <- matrix(nrow=isnit_max,ncol=L) # object to store class proportions at each iteration
      d    <- c()                                 # loglikelihood difference between two cycles

	  it=1

	  oldpar[[it]]    = as.numeric(start.parameters[set,])

        # Starting parameters for covariances are generated pairwise 
        # above, so it sometimes makes the starting variance-covariance 
        # matrix of random effects not positive definite
        # Below I check for it, and if necessary, smooth the matrix
        # to make it positive definite
        # I replace the starting values for covariance terms after smoothing 


		est.cov <- vector("list",L)

		if(fix.cov==FALSE & fix.err == FALSE) {

			for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
            	                                          13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
			}
			
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6+15*(u-1),7+15*(u-1),
                                                               8+15*(u-1),9+15*(u-1))]))    # Transform the smoothed correlation matrix back to the covariance matrix
      	            oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
            	                     13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))] <- vech(covv)
		  	}
			}
		}


		if(fix.cov==TRUE & fix.err == FALSE) {

				for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
			}
		
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6,7,8,9)]))                   # Transform the smoothed correlation matrix back to the covariance matrix
                 	 oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

		}

		if(fix.cov==FALSE & fix.err == TRUE) {

				params <- vector("list",L)    
	     		      params[[1]] <- oldpar[[it]][1:15]
	
				if(L>=2) {
      			  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
			        for(u in 2:L) {
      			   params[[u]] = params[[1]]  
			         params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
				   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
		  		  }
				}


			for(u in 1:L) {			
				est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
			}
			
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(params[[u]][c(6,7,8,9)]))    # Transform the smoothed correlation matrix back to the covariance matrix
      	            params[[u]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

				oldpar[[it]][1:15] <- params[[1]]
	
				if(L>=2) {
      			  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
			        for(u in 2:L) { 
			         oldpar[[it]][(locate[u-1]):(locate[u-1]+3)] <- params[[u]][1:4] 
				   oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)] <- params[[u]][6:15] 
		  		  }
				}

		}

		if(fix.cov==TRUE & fix.err == TRUE) {

				for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
			}
		
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6,7,8,9)]))                   # Transform the smoothed correlation matrix back to the covariance matrix
                 	 oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

		}





	  class.prop[it,] = pik
	  LogL            = likelihood(data=data,start.param=oldpar[[it]],nquad=nquad,L=L,pis=pik)

	  l[it]           = LogL$LL
	  posterior[[it]] = LogL$post.prob
	  g[[it]]         = grad(data,oldpar[[it]] ,nquad=nquad,L=L,pis=pik)
	  H[[it]]         = nhess(data,oldpar[[it]] ,nquad=nquad,L=L,pis=pik)
	  IH[[it]]        = solve(H[[it]])
	  maxg[it]        = max(abs(g[[it]]),na.rm=TRUE)
        d[it]           = 5

	  sel = 2     # Choose either BFGS(1) or DFPQ(2) method to update the parameters


	while(it < isnit_max & (d[it] > is_st | d[it] < 0) & is.na(l[it])!=TRUE & is.infinite(l[it])!=TRUE) {   
             
            # This while loop continues until it reaches the maximum number of iterations
      	# or when the difference between loglikelihoods is bigger than "st" defined above
            # or when the estimated var-cov matrices are not positive definite 
               # after reducing the step size up to 1/10000

            # This while loop stop if any of these criteria is not met.
 
		  it <- it+1

		  step.size <- 1/(1:10000)   # These are the step sizes for gradient

	 		if(sel==1){ 

		 	     step=1
      		     ev=-1
				    while(sum(ev<0)!=0 & step <= 10000) {  # check whether the estimated var-cov matrix
                  	                                         # is positive definite for each iteration
                        	                                   # If not, reduce the step size until to get a 
                              	                             # positive definite var-cov matrix of random effects

	                  	oldpar[[it]]= oldpar[[it-1]]-(solve(H[[it-1]])%*%as.matrix(g[[it-1]]*step.size[step]))
					est.cov <- vector("list",L)

					if(fix.cov==FALSE & fix.err == FALSE) {

						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
				                                                      13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
						for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
					}
				
					if(fix.cov==TRUE & fix.err == FALSE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}
				
					if(fix.cov==FALSE & fix.err == TRUE) {
			
						params <- vector("list",L)    
			     		      params[[1]] <- oldpar[[it]][1:15]
	
						if(L>=2) {
      					  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
					        for(u in 2:L) {
      					   params[[u]] = params[[1]]  
			      		   params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
						   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
				  		  }
						}


						for(u in 1:L) {			
						est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
						}
			
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

					if(fix.cov==TRUE & fix.err == TRUE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

				    }
      	      }
	  	
	      	if(sel==2){ 

		        
		 	     step=1
      		     ev=-1
				    while(sum(ev<0)!=0 & step <= 10000) {  # check whether the estimated var-cov matrix
                  	                                         # is positive definite for each iteration
                        	                                   # If not, reduce the step size until to get a 
                              	                             # positive definite var-cov matrix of random effects

	                  	oldpar[[it]]= oldpar[[it-1]]-(IH[[it-1]]%*%as.matrix(g[[it-1]]*step.size[step])) 
					est.cov <- vector("list",L)

					if(fix.cov==FALSE & fix.err == FALSE) {

						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
				                                                      13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
						for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
					}
				
					if(fix.cov==TRUE & fix.err == FALSE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}
				
					if(fix.cov==FALSE & fix.err == TRUE) {
			
						params <- vector("list",L)    
			     		      params[[1]] <- oldpar[[it]][1:15]
	
						if(L>=2) {
      					  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
					        for(u in 2:L) {
      					   params[[u]] = params[[1]]  
			      		   params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
						   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
				  		  }
						}


						for(u in 1:L) {			
						est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
						}
			
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

					if(fix.cov==TRUE & fix.err == TRUE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

				    }
      	      }
    
	        class.prop[it,] <- colMeans(posterior[[it-1]])   # Update class proportions from the previous iteration
                                                               # Page 188
		

		  LogL            = likelihood(data,oldpar[[it]],nquad=nquad,L=L,pis=class.prop[it,])
		    l[it]           = LogL$LL
		    posterior[[it]] = LogL$post.prob

		if(is.na(l[it])!=TRUE & is.infinite(l[it])!=TRUE) {

		  g[[it]]         = grad(data,oldpar[[it]],nquad=nquad,L=L,pis=class.prop[it,])
      	  maxg[it]        = max(g[[it]],na.rm=TRUE)
	        gradiff  = g[[it]] - g[[it-1]]
		  pardiff  = oldpar[[it]]-oldpar[[it-1]]
              d[it]    = l[it]-l[it-1]

	          # Broyden-Fletcher-Goldfarb-Shanno Quasi-Newton Method for updating the Hessian matrix #

      	        Hstep  = ((H[[it-1]]%*%pardiff%*%t(pardiff)%*%H[[it-1]])/as.numeric((t(pardiff)%*%H[[it-1]]%*%pardiff))) -
            	           ((gradiff%*%t(gradiff))/as.numeric(t(pardiff)%*%gradiff))
	              H[[it]]= H[[it-1]] - Hstep

      	    # Davidon-Fletcher-Powell Quasi-Newton Method for updating the inverse of the Hessian matrix #

	      	  iHstep  = ((pardiff%*%t(pardiff))/as.numeric(t(pardiff)%*%gradiff)) - 
                  	     ((IH[[it-1]]%*%gradiff%*%t(gradiff)%*%IH[[it-1]])/as.numeric(t(gradiff)%*%IH[[it-1]]%*%gradiff))  

	              IH[[it]]  = IH[[it-1]] + iHstep  
           
	    }  

	  	if(it==2) { 
		   cat(sprintf("%12s %10s %16s %10s %12s","Iteration","LL","Max.Gradient","LL.Change","Step Size"),"\n")
		   cat(sprintf("%8.0f %16.3f %11.3f %9s %12.0f",it-1,l[it-1],maxg[it-1],"-",step),"\n")
		  }   
		   cat(sprintf("%8.0f %16.3f %11.3f %9.3f %12.0f",it,l[it],maxg[it],d[it],step),"\n")
	}

	output <- as.data.frame(matrix(nrow=it,ncol=4+L))
	nnn    <- c("iter","LL","maxg","change")
	for(u in 1:L) { nnn <- c(nnn,paste("Class",u,".Prop",sep="")) }
	colnames(output) <- nnn
	for(ij in 1:it) {
		output[ij,1]=ij
		output[ij,2]=l[ij]
		output[ij,3]=maxg[ij]
 		output[ij,4]=d[ij]
		for(u in 1:L) {
			output[ij,4+u]=class.prop[ij,u]
		}

	}

  out <- as.data.frame(output)
  write.fwf(x=out,
            file=paste("SET",set,"_ITERATIONS.txt",sep=""),
            rownames=FALSE,sep="   ")
      
RESULTS.INITIAL[[set]] <- output

}
	
max.LL <- c()
for(set in 1:is_nset) { max.LL[set]= na.omit(RESULTS.INITIAL[[set]])[nrow(na.omit(RESULTS.INITIAL[[set]])),2] }
max.LL <- cbind(1:is_nset,max.LL)
max.LL <- max.LL[order(max.LL[,2],decreasing=TRUE),]
if(length(which(max.LL[,2]>0))!=0) { max.LL <- max.LL[-which(max.LL[,2]>0),] }
best  <- max.LL[1:fs_nset,1]

cat("","\n")
cat("**************************************************************************","\n") 
cat("**************************************************************************","\n") 
cat("","\n")
cat(" Initial stage estimation is completed.","\n") 
cat("","\n")
cat(paste(" The Best ",fs_nset," Start Value Sets:",sep=""),"\n") 
cat("","\n")
cat("          SET      LL","\n") 
for(i in 1:fs_nset){
lim=abs(round(max.LL[i,2]/1000,0))+1
cat("     ",sprintf(paste("%6.0f %",lim+6,".3f",sep=""),max.LL[i,1],max.LL[i,2]),"\n")
}
cat(" ","\n")
cat("**************************************************************************","\n") 
##########               FINAL STAGE ITERATION               #############
cat("","\n")
cat("                  FINAL STAGE ESTIMATION              ","\n")
cat("","\n")
cat("The results from the final stage estimation are saved in the following folder:","\n") 
cat(getwd(),"\n")
cat("","\n")

OUT       <- vector("list",length(best))
PAR       <- vector("list",length(best))
POSTERIOR <- vector("list",length(best))
CPHIST    <- vector("list",length(best))

for(best.set in 1:length(best)) {

	cat("","\n")
	cat("START VALUE SET ",best[best.set],":","\n")

	####   First Iteration ##########
               
	l <- c()                              # object to store the loglikelihood values for the data for each iteration
	oldpar <- vector("list",fsnit_max)          # list object to store the estimates parameters at each iteration
	g <- vector("list",fsnit_max)               # object to store the gradient vector at each iteration  
	H <- vector("list",fsnit_max)               # object to store the Hessian matrix at each iteration  
	IH <- vector("list",fsnit_max)              # object to store the inverse of the Hessian matrix at each iteration  
	maxg <- c()                           # object to store the maximum gradient value at each iteration
	posterior <- vector("list",fsnit_max)       # object to store individual posterior probabilities at each iteration
	class.prop <- matrix(nrow=fsnit_max,ncol=L) # object to store class proportions at each iteration
      d    <- c()                           # loglikelihood difference between two cycles

	  it=1

	  oldpar[[it]]    = as.numeric(start.parameters[best[best.set],])

	   	est.cov <- vector("list",L)

		if(fix.cov==FALSE & fix.err == FALSE) {

			for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
            	                                          13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
			}
			
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6+15*(u-1),7+15*(u-1),
                                                               8+15*(u-1),9+15*(u-1))]))    # Transform the smoothed correlation matrix back to the covariance matrix
      	            oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
            	                     13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))] <- vech(covv)
		  	}
			}
		}


		if(fix.cov==TRUE & fix.err == FALSE) {

				for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
			}
		
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6,7,8,9)]))                   # Transform the smoothed correlation matrix back to the covariance matrix
                 	 oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

		}

		if(fix.cov==FALSE & fix.err == TRUE) {

				params <- vector("list",L)    
	     		      params[[1]] <- oldpar[[it]][1:15]
	
				if(L>=2) {
      			  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
			        for(u in 2:L) {
      			   params[[u]] = params[[1]]  
			         params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
				   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
		  		  }
				}


			for(u in 1:L) {			
				est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
			}
			
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(params[[u]][c(6,7,8,9)]))    # Transform the smoothed correlation matrix back to the covariance matrix
      	            params[[u]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

				oldpar[[it]][1:15] <- params[[1]]
	
				if(L>=2) {
      			  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
			        for(u in 2:L) { 
			         oldpar[[it]][(locate[u-1]):(locate[u-1]+3)] <- params[[u]][1:4] 
				   oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)] <- params[[u]][6:15] 
		  		  }
				}

		}

		if(fix.cov==TRUE & fix.err == TRUE) {

				for(u in 1:L) {			
				est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
			}
		
			eig <- vector("list",L)
			for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }

			for(u in 1:L) {

			  if(length(which(eig[[u]]<0))!=0) {

				corr <- cov2cor(est.cov[[u]])                                             # transform the covariance matrix to the correlation matrix
				corr <- cor.smooth(corr)                                                  # Smooth the correlation matrix to make it positive definite if it's not positive definite.
				covv  <- cor2cov(corr,sd=sqrt(oldpar[[it]][c(6,7,8,9)]))                   # Transform the smoothed correlation matrix back to the covariance matrix
                 	 oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)] <- vech(covv)
		  	}
			}

		}



	  class.prop[it,] = pik
	  LogL            = likelihood(data,oldpar[[it]] ,nquad=nquad,L=L,pis=pik)

	  l[it]           = LogL$LL
	  posterior[[it]] = LogL$post.prob
	  g[[it]]         = grad(data,oldpar[[it]] ,nquad=nquad,L=L,pis=pik)
	  H[[it]]         = nhess(data,oldpar[[it]] ,nquad=nquad,L=L,pis=pik)
	  IH[[it]]        = solve(H[[it]])
	  maxg[it]        = max(abs(g[[it]]),na.rm=TRUE)
        d[it]           = 5


	while(it < fsnit_max & (d[it] > fs_st | d[it] < 0) & is.na(l[it])!=TRUE & is.infinite(l[it])!=TRUE) {   
             
            # This while loop continues until it reaches the maximum number of iterations
      	# or when the difference between loglikelihoods is bigger than .1
            # or when the estimated var-cov matrices are not positive definite 
               # after reducing the step size 1/10000
 
		  it <- it+1

			step.size <- 1/(1:10000)   # These are the step sizes for gradient

	 		if(sel==1){ 

		 	     step=1
      		     ev=-1
				    while(sum(ev<0)!=0 & step <= 10000) {  # check whether the estimated var-cov matrix
                  	                                         # is positive definite for each iteration
                        	                                   # If not, reduce the step size until to get a 
                              	                             # positive definite var-cov matrix of random effects

	                  	oldpar[[it]]= oldpar[[it-1]]-(solve(H[[it-1]])%*%as.matrix(g[[it-1]]*step.size[step]))
					est.cov <- vector("list",L)

					if(fix.cov==FALSE & fix.err == FALSE) {

						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
				                                                      13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
						for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
					}
				
					if(fix.cov==TRUE & fix.err == FALSE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}
				
					if(fix.cov==FALSE & fix.err == TRUE) {
			
						params <- vector("list",L)    
			     		      params[[1]] <- oldpar[[it]][1:15]
	
						if(L>=2) {
      					  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
					        for(u in 2:L) {
      					   params[[u]] = params[[1]]  
			      		   params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
						   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
				  		  }
						}


						for(u in 1:L) {			
						est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
						}
			
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

					if(fix.cov==TRUE & fix.err == TRUE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

				    }
      	      }
	  	
	      	if(sel==2){ 

		        
		 	     step=1
      		     ev=-1
				    while(sum(ev<0)!=0 & step <= 10000) {  # check whether the estimated var-cov matrix
                  	                                         # is positive definite for each iteration
                        	                                   # If not, reduce the step size until to get a 
                              	                             # positive definite var-cov matrix of random effects

	                  	oldpar[[it]]= oldpar[[it-1]]-(IH[[it-1]]%*%as.matrix(g[[it-1]]*step.size[step])) 
					est.cov <- vector("list",L)

					if(fix.cov==FALSE & fix.err == FALSE) {

						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6+15*(u-1),10+15*(u-1),11+15*(u-1),12+15*(u-1),7+15*(u-1),
				                                                      13+15*(u-1),14+15*(u-1),8+15*(u-1),15+15*(u-1),9+15*(u-1))])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
						for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
					}
				
					if(fix.cov==TRUE & fix.err == FALSE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}
				
					if(fix.cov==FALSE & fix.err == TRUE) {
			
						params <- vector("list",L)    
			     		      params[[1]] <- oldpar[[it]][1:15]
	
						if(L>=2) {
      					  locate <- seq(from=16,to=length(oldpar[[it]]),by=14)
					        for(u in 2:L) {
      					   params[[u]] = params[[1]]  
			      		   params[[u]][1:4] <- oldpar[[it]][(locate[u-1]):(locate[u-1]+3)]
						   params[[u]][6:15] <- oldpar[[it]][(locate[u-1]+4):(locate[u-1]+13)]
				  		  }
						}


						for(u in 1:L) {			
						est.cov[[u]] <-  xpnd(params[[u]][c(6,10,11,12,7,13,14,8,15,9)])
						}
			
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }	

						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

					if(fix.cov==TRUE & fix.err == TRUE) {
						for(u in 1:L) {			
							est.cov[[u]] <- xpnd(oldpar[[it]][c(6,10,11,12,7,13,14,8,15,9)])
						}
		
						eig <- vector("list",L)
						for(u in 1:L) { eig[[u]] = eigen(est.cov[[u]])$values }
						ev <- eig[[1]]
						if(L>=2) { 
							for(u in 2:L) { ev <- c(ev,eig[[u]]) }
						}
            	      		step=step+1
						
					}

				    } 
      	      }

	        class.prop[it,] <- colMeans(posterior[[it-1]])   # Update class proportions from the previous iteration
                                                               # Page 188
		

		  LogL            = likelihood(data,oldpar[[it]],nquad=nquad,L=L,pis=class.prop[it,])
		    l[it]           = LogL$LL
		    posterior[[it]] = LogL$post.prob

		if(is.na(l[it])!=TRUE & is.infinite(l[it])!=TRUE) {

		  g[[it]]         = grad(data,oldpar[[it]],nquad=nquad,L=L,pis=class.prop[it,])
      	  maxg[it]        = max(g[[it]],na.rm=TRUE)
	        gradiff  = g[[it]] - g[[it-1]]
		  pardiff  = oldpar[[it]]-oldpar[[it-1]]
              d[it]    = l[it]-l[it-1]

	          # Broyden-Fletcher-Goldfarb-Shanno Quasi-Newton Method for updating the Hessian matrix #

      	        Hstep  = ((H[[it-1]]%*%pardiff%*%t(pardiff)%*%H[[it-1]])/as.numeric((t(pardiff)%*%H[[it-1]]%*%pardiff))) -
            	           ((gradiff%*%t(gradiff))/as.numeric(t(pardiff)%*%gradiff))
	              H[[it]]= H[[it-1]] - Hstep

      	    # Davidon-Fletcher-Powell Quasi-Newton Method for updating the inverse of the Hessian matrix #

	      	  iHstep  = ((pardiff%*%t(pardiff))/as.numeric(t(pardiff)%*%gradiff)) - 
                  	     ((IH[[it-1]]%*%gradiff%*%t(gradiff)%*%IH[[it-1]])/as.numeric(t(gradiff)%*%IH[[it-1]]%*%gradiff))  

	              IH[[it]]  = IH[[it-1]] + iHstep  
           
	    }  

	  	if(it==2) { 
		   cat(sprintf("%12s %10s %16s %10s %12s","Iteration","LL","Max.Gradient","LL.Change","Step Size"),"\n")
		   cat(sprintf("%8.0f %16.3f %11.3f %9s %12.0f",it-1,l[it-1],maxg[it-1],"-",step),"\n")
		  }   
		   cat(sprintf("%8.0f %16.3f %11.3f %9.3f %12.0f",it,l[it],maxg[it],d[it],step),"\n")
	}


	output <- as.data.frame(matrix(nrow=it,ncol=4+L))
	nnn    <- c("iter","LL","maxg","change")
	for(u in 1:L) { nnn <- c(nnn,paste("Class",u,".Prop",sep="")) }
	colnames(output) <- nnn
	for(ij in 1:it) {
		output[ij,1]=ij
		output[ij,2]=l[ij]
		output[ij,3]=maxg[ij]
 		output[ij,4]=d[ij]
		for(u in 1:L) {
			output[ij,4+u]=class.prop[ij,u]
		}

	}

	out <- as.data.frame(output)
	OUT[[best.set]] <- out

	#  setwd("\\\\ad.umn.edu\\CEHD\\Projects\\EdPsy\\NKohli\\Cengiz\\Amanda\\Math\\Two-classPME_withmissing")

  	write.fwf(x=out,
            file=paste("FINAL_STAGE_ITERATION_SET",best[best.set],".txt",sep=""),
            rownames=FALSE,sep="   ")
 
	Hess   <- nhess(data,oldpar[[it]],nquad=nquad,L=L,pis=pik)       # Hessian matrix at the last step 
	infor  <- -Hess                                                  # Inverse of the negative of Hessian
	SE     <- sqrt(diag(solve(infor))) 
	param  <- cbind(oldpar[[it]],SE)


	if(fix.cov==FALSE & fix.err == FALSE) {
		ccc <- rep("1",15)
		for(u in 2:L) { ccc <- c(ccc,rep(paste(u,sep=""),15))}
	
		rnames <- c("B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)",
                  "Residual Variance","B1 variance","B2 variance","B3 variance","K variance",
                  "B1B2 Covariance","B1B3 Covariance","B1K Covariance","B2B3 Covariance","B2K Covariance","B3K Covariance")
		rnames <- rep(rnames,3)
	}

	if(fix.cov==TRUE & fix.err == FALSE){
		ccc <- rep("1",15)
		for(u in 2:L) { ccc <- c(ccc,rep(paste(u,sep=""),5))}
	
		rnames <- c("B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)",
                  "Residual Variance","B1 variance","B2 variance","B3 variance","K variance",
                  "B1B2 Covariance","B1B3 Covariance","B1K Covariance","B2B3 Covariance","B2K Covariance","B3K Covariance")

		for(u in 2:L) { rnames <- c(rnames,"B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)","Residual Variance")}
	}


	if(fix.cov==FALSE & fix.err == TRUE){
		ccc <- rep("1",15)
		for(u in 2:L) { ccc <- c(ccc,rep(paste(u,sep=""),14))}
	
		rnames <- c("B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)",
                  "Residual Variance","B1 variance","B2 variance","B3 variance","K variance",
                  "B1B2 Covariance","B1B3 Covariance","B1K Covariance","B2B3 Covariance","B2K Covariance","B3K Covariance")

		for(u in 2:L) { rnames <- c(rnames,rnames[-5])}
	}

	if(fix.cov==TRUE & fix.err == TRUE){
		ccc <- rep("1",15)
		for(u in 2:L) { ccc <- c(ccc,rep(paste(u,sep=""),4))}
	
		rnames <- c("B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)",
                  "Residual Variance","B1 variance","B2 variance","B3 variance","K variance",
                  "B1B2 Covariance","B1B3 Covariance","B1K Covariance","B2B3 Covariance","B2K Covariance","B3K Covariance")

		for(u in 2:L) { rnames <- c(rnames,"B1 (Intercept)","B2 (Slope 1)","B3 (Slope 2)","K (Knot)")}
	}


	
	parameters <- as.data.frame(cbind(ccc,rnames,round(param,6)))
	colnames(parameters) <- c("Class","Parameter","Est","SE")
	PAR[[best.set]] <- parameters

      write.fwf(x=parameters,
            file=paste("FINAL_PARAMETER_EST_SET",best[best.set],".txt",sep=""),
            rownames=TRUE,sep="   ")

	CPHIST[[best.set]] <- as.data.frame(na.omit(class.prop))
	write.fwf(x=as.data.frame(na.omit(class.prop)),
                file=paste("CLASS_PROP_EST_HIST_SET",best[best.set],".txt",sep=""),
                rownames=TRUE,sep="   ")

	POSTERIOR[[best.set]] <- as.data.frame(round(posterior[[it]],4))
      write.fwf(x=as.data.frame(round(posterior[[it]],4)),
                file=paste("POSTERIOR_PROB_EST_SET",best[best.set],".txt",sep=""),
                rownames=TRUE,sep="   ")
}
cat("**************************************************************************","\n") 

best.out <- c()
for(i in 1:length(best)) { best.out[i] = round(na.omit(OUT[[i]])[nrow(na.omit(OUT[[i]])),2],3) }
best.out <- cbind(best,best.out)
if(fs_nset>=2) { best.out <- best.out[order(best.out[,2],decreasing=TRUE),]}

cat("","\n")
cat(" Final stage estimation is completed.","\n")
cat("","\n")
cat(" Random Starts Results Ranked from the Best to the Worst Loglikelihood Values","\n")
cat("","\n")
cat("   START VALUE SET        LL  ","\n")
for(i in 1:length(best)){
lim=abs(round(best.out[i,2]/1000,0))+1
cat("     ",sprintf(paste("%6.0f %",lim+15,".3f",sep=""),best.out[i,1],best.out[i,2]),"\n")
}
cat("","\n")

if(sum(best.out==max(best.out))<2) {
cat("WARNING: The best loglikelihood was not replicated. The final parameter","\n")
cat("estimates may not be trustworthy due to local maxima. Increase the number","\n")
cat("of start value sets.","\n")
}

if(sum(best.out==max(best.out))>=2) {
cat("WARNING: The best loglikelihood has been replicated. You may want to re-run","\n")
cat("one or two more times to check that the best loglikelihood is still obtained and","\n")
cat("replicated","\n")
}

cat("","\n")
cat(paste("        MODEL FIT INFORMATION FOR START VALUE SET ",best.out[1,1],sep=""),"\n")
cat("","\n")
cat("Number of Parameters:",nrow(PAR[[which(best==best.out[1,1])]]),"\n")
cat("Loglikelihood Value :",best.out[1,2],"\n")
cat("                AIC :",-2*best.out[1,2]+(nrow(PAR[[which(best==best.out[1,1])]])),"\n")
cat("                BIC :",-2*best.out[1,2]+(nrow(PAR[[which(best==best.out[1,1])]]))*log(nrow(data)),"\n")
cat("","\n")

if(L>=2){

	cat("Final Class Proportions:","\n")
	cat("","\n")
	for(i in 1:L){
		CM <- as.numeric(CPHIST[[which(best==best.out[1,1])]][nrow(CPHIST[[which(best==best.out[1,1])]]),])
		cat(paste("    Class ",i,"   :",sep=""),CM[i],"\n")
	}
}

cat("","\n")
cat("        PARAMETER ESTIMATES","\n")
cat("","\n")

par <- PAR[[which(best==best.out[1,1])]]

cat("    Class     Parameter        Estimate      SE","\n")
for(i in 1:nrow(par)) {
cat(sprintf("%-10s %-18s %8s %8s",paste("  ",par[i,1],sep="    "),par[i,2],par[i,3],par[i,4]),"\n")
}
cat("","\n")
if(sum(par[,4]=="NaN")>0) {
cat("WARNING: The standard errors are not available due to the negative values at the diagonal","\n")
cat("of the inverse of the Hessian matrix at the last iteration.","\n")
}
cat("","\n")

aaa.post <- as.data.frame(POSTERIOR[[which(best==best.out[1,1])]])
aaa.name <- "P1"
if(ncol(aaa.post)>=2){
for(i in 2:ncol(aaa.post)) { aaa.name <- c(aaa.name, paste("P",i,sep=""))}
}
colnames(aaa.post) <- aaa.name
data2 <- cbind(data,aaa.post)
data2$Class <- NA
for(i in 1:nrow(data2)) { data2[i,]$Class <- which.max(aaa.post[i,])}

# Plotted Observed and Fitted Curves

class.means <- matrix(nrow=L,ncol=ncol(data))
for(i in 1:L) { 
	class.means[i,]=colMeans(data2[which(data2$Class==i),1:ncol(data)],na.rm=TRUE)
}

plot(time.vec,class.means[1,],ylim=c(0,max(data)),xlab="Time",ylab="Outcome",pch=0,cex=1)
if(L>=2){
for(i in 2:L) {
	 points(time.vec,class.means[i,],ylim=c(0,max(data)),xlab="Time",ylab="Outcome",pch=i-1,cex=1)
}
}

ltitle <- "Observed Mean for Estimated Class 1 Members"
if(L>=2){
	for(i in 2:L) { ltitle <- c(ltitle,paste("Observed Mean for Estimated Class ",i," Members",sep=""))}
}

legend("topleft",ltitle,pch=0:(L-1),cex=.75)

plotPME <- function (B1,B2,B3,K,t,LTY) {
 a <- seq(from=min(t),to=max(t),by=.001)
 p <- c()
 for(i in 1:length(a)) { p[i] <- ifelse(a[i]<K,B1 + B2*a[i],B1 + B2*K + B3*(a[i]-K)) } 
 points(a,p,type="l",lty=LTY)
}

nclass <- L
for(i in 1:nclass) {
	B1=as.numeric(as.character(par[which(par[,1]==i),][1,3]))
	B2=as.numeric(as.character(par[which(par[,1]==i),][2,3]))
	B3=as.numeric(as.character(par[which(par[,1]==i),][3,3]))
	K =as.numeric(as.character(par[which(par[,1]==i),][4,3]))
      t =time.vec
	plotPME(B1,B2,B3,K,t,LTY=i)
}

ltitle <- "Fitted Curve for Class 1"
if(L>=2){
	for(i in 2:L) { ltitle <- c(ltitle,paste("Fitted Curve for Class ",i,sep=""))}
}

legend("bottomright",ltitle,lty=1:L,cex=.75)

return(list(PAR       = par,
		CP        = as.numeric(CPHIST[[which(best==best.out[1,1])]][nrow(CPHIST[[which(best==best.out[1,1])]]),]),
		POSTERIOR = data2
		)
	 )

}













