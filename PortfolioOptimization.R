library(data.table)
library(openxlsx)
library(plotly)
library(quantmod)
library(nloptr)
library(quadprog)


################################Optimal Investing Starts########################################
# Notes: To calculate the short sale vs no short sale, I used two method here.
# An easier way is to simply change the lb in the noshortsale function change the constraints.


# Start by choosing 15 equities
tickers <- c("SBUX", "MSFT", "JPM", "BAC",'AAPL',
             "XOM", "F", "MA", "IBM", 'MCD',
             "WFC", "PG", "GE", "BA", 'WMT')

# Loop to download stock data from yahoo
# The date range is 2010-01-01 to 2022-12-31
# Here, I only take the sixth column, which is the adjusted close to perform analysis
stock_data <- list()
for(ticker in tickers) {
  data <- getSymbols(ticker, src = "yahoo", from = "2010-01-01", to = "2020-12-31", auto.assign = FALSE)
  stock_data[[ticker]] <- data[,6]
}

# Combine all into a data frame and check if there is NA
combined_stock_data <- do.call(merge, stock_data)
sum(is.na(combined_stock_data))

# Calculating the stock returns
# The frequency of returns is daily
# Used diff with log returns here and transform the returns back to normal
log_returns <- diff(log(combined_stock_data))[-1]
returns     <- exp(log_returns)-1


# So there are two ways that we can go in calculating the covariance matrix.
# Showing one example here:
# To calculate the covariance of SBUX and MSFT, the following formula is used:
sum((returns[,1]-mean(returns[,1]))*(returns[,2]-mean(returns[,2])))/length(returns[,1])

# Here I used the cov function, which easily help calculate the covariance matrix
# Calculating the return and covariance matrix and annualized by multiply 252 trading days
avg_returns <- colMeans(returns)
mu <- avg_returns * 252

covariance_matrix <- cov(returns)
cov  <- covariance_matrix * 252

sd_returns <- apply(returns,2,sd)*sqrt(252)
cor(returns)

library(corrplot)
# Visualizing the heatmap of annualized covariance
corrplot(cov, method = "color", 
         type = "lower", order = "hclust",addCoef.col = "gray")

# Visualizing the correlation helps knowing more about the data
corr_matrix <- cor(returns)
corrplot(corr_matrix, method = "color", 
         type = "lower", order = "hclust",addCoef.col = "gray")

# The efficient frontier is essentially the most optimal risk-reward option for investor to invest in.
# In this case, it will consist of the 15 stocks above and creates a frontier. The investor will be able
# to earn the best return at each risk by investing in the portfolio of 15 assets.
# Before creating frontier, I plotted the mu and sd of each stock.
names(mu) = tickers
sd.vec = sqrt(diag(cov)) 
cex.val = 1
plot(sd.vec, mu,  ylim=c(0, 0.50), xlim=c(0, 0.50), 
     ylab=expression(mu[p]), xlab=expression(sigma[p]), 
     pch=16, col="blue", cex=cex.val, cex.lab=1)      
text(sd.vec, mu, labels=names(mu), pos=4, cex = cex.val) 


############################# Used Nloptr to create no short-sale MVF ##############################
# Setting n as the length of returns, which is 15 in this case.
n = length(mu)

# Created a function for no short sale
noshortsale <- function(y){
# Objective function
    obj <- function(w) {
            n <- length(mu)
            w <- w / sum(w) # ensure weights sum to 1
            ret <- sum(mu * w) # portfolio return
            var <- t(w) %*% cov %*% w # portfolio variance
            penalty <- (ret - target_return)^2 # quadratic penalty on deviation from target return
            obj <- var + penalty # objective function
            grad <- 2 * (cov %*% w) + 2 * (ret - target_return) * mu # gradient of objective function
            return(list(objective = obj, gradient = grad))
            }

    # Setting up the constraints
    eval_g <- function(w) {
               return(sum(w)-1)
            }

    eval_jac_g <- function( w ) {
                return( rep(-1,n) )
    # this step returns the vector of first derivatives of eval_g(x) w.r.t. x, which is a vector of -1.
             }

    # Setting the variables
     w0 <- rep(1/n,n) #Created a starting weight
     lb <- rep(0, n) # Setting the lower bound to be 0
     ub <- rep(Inf, n) # Setting the upper bound
     target_return <- y # the y is the target return, which is used in the no short sale function

     # Setting the optimization
     opt <- nloptr(x0 = w0,
              eval_f = obj,
              lb = lb,
              ub = ub,
              eval_g_eq = eval_g,              # tells it that the equality constraint is eval_g
              eval_jac_g_eq = eval_jac_g,
              opts = list("algorithm" = "NLOPT_LD_SLSQP",         # opts is a list of controls that tells nloptr what to do.  "algorithm" is the method nloptr uses to minizmie the objective function
                          "xtol_rel"=1.0e-10,                     # xtol_rel is a tolerance level
                          "print_level" = 0,                      # print_level tells nloptr how much to print out while it is running
                          "check_derivatives" = TRUE,             # this tells nloptr to check the derivatives (using numeric derivatives)
                          "check_derivatives_print" = "all"))

     # Here the expected return is calculated after optimizing the weights with minimum variance
      expected_return <- sum(opt$solution*mu)
      sd <-  sqrt(t(opt$solution) %*% cov %*% opt$solution)
      
      # The function will return the expected return and its standard deviation accordingly.
      return(c(expected_return,sd))

}

# Creating a simple dataframe to store the ER and SD for no short sale
df_noss <- data.frame(0,0)
amt <- seq(.05,1,by = .05)
for (i in amt){
df_noss <- rbind(df_noss,noshortsale(i))
}
df_noss <- df_noss[-1,]
colnames(df_noss) <- c('ER','SD')

# Plotting the No Short Sale Mean Variance Frontier
plot(df_noss$SD,df_noss$ER,ylim=c(0, 0.45), xlim=c(0, 0.45), type ='l')
############################ End: Used Nloptr to create no short-sale MVF  ######################


############################ Start: Used quadprog to create short sale MVF ######################
# Can also use above "noshortsale" function, just need to change lb to -Inf
# Optimization with short sales allow
# Here we used the quadratic to find the weights, portfolio's ER and SD
# Setting a target return list
target_returns <- seq(0.05, 0.45, by = 0.02)

# Create a matrix with row according to the 15 tickers
# The columns are equal to the length of the target returns
portfolio_weights <- matrix(NA, nrow = length(mu), ncol = length(target_returns))

# Create a for loop to run the optimization and store the weights into the portfolio_weights matrix
for (i in seq_along(target_returns)) {
  # Define the objective function and constraints
  dvec <- rep(0, length(mu))
  Dmat <- cov
  Amat <- cbind(mu, rep(1, length(mu)))
  bvec <- c(target_returns[i], 1)

  # Solve the optimization problem
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  
  # Store the weights for this target return
  portfolio_weights[, i] <- sol$solution
}

# Calculate the portfolio returns and standard deviations for each set of weights
portfolio_returns <- mu %*% portfolio_weights
portfolio_stddevs <- sqrt(diag(t(portfolio_weights) %*% cov %*% portfolio_weights))

# Plot the EF with short sale
par(mfrow = c(1, 1))
plot(portfolio_stddevs, portfolio_returns,ylim=c(0, 0.45), xlim=c(0, 0.45), type = "l", lwd = 2, xlab = "Portfolio Standard Deviation", ylab = "Portfolio Return")
# Plot the EF with no short sale using points
points(df_noss$SD,df_noss$ER,ylim=c(0, 0.45), xlim=c(0, 0.45))

# This is to plot the individual stocks returns and standard deviation
points(sd.vec, mu,  ylim=c(0, 0.50), xlim=c(0, 0.50), 
     ylab=expression(mu[p]), xlab=expression(sigma[p]), 
     pch=16, col="blue", cex=cex.val, cex.lab=1)   
text(sd.vec, mu, labels=names(mu), pos=4, cex = cex.val) 
################################ Optimal Investing End ########################################

# Efficient frontier- it is the best combination risk-reward that investor can create using all the risky assets
# in this case, the 15 different stocks. Instead of investing in each stock, investor can invest in different
# weights on each stock and create the best return under each risk(sd). The weights
# are essentially the combination that form on the efficient frontiers.

# The no short sale constraint pull the frontier down. Because of the constraint, there is limitation on
# how we can assign the weights on each stock. If we have short sale, we can short some stock and
# long in the others.

# If we are calculating covariance by hand, then it will be a lot of works!
# However, there is a simple function in R, which is cov(), this creates a cov matrix easily.
# With this, we do not need to worry about how many equities we are using.






###Part 2#####
################################ Linear Regression Starts ################################
# 2a
# Reading both files into R
infotxt <- fread('info.txt')
performance_df <- fread('performance.csv')

# Perform left join on two tables
Merged_df <- infotxt[performance_df, on ="ID" ]
Final_Merged_df <- Merged_df[,c("Date","Performance","Name")]

#2b
# Write the final dataframe to excel
fwrite(Final_Merged_df, "Final_Merged_df.xlsx")

#2c
# Regressing StockCZ against Index
# Preparing data
Index_data <- Final_Merged_df[Final_Merged_df$Name == 'Index', c(1:2)]
StockCZ_data <- Final_Merged_df[Final_Merged_df$Name == 'Stock CZ', c(1:2)]
Combined <- Index_data[StockCZ_data, on = 'Date']
colnames(Combined) <- c('Date', 'Index', 'Stock_CZ')

# Linear Regression
model <- lm(Stock_CZ ~ Index, data = Combined)

# Looking at the p value(0.00) and t stat(68.295) of coefficient, we can tell that the it is significant.
# The beta is 0.8069
# Essentially, this means that the stock is less volatile than market.
# If the market drop by 1%, the stock drop by 0.8%
summary(model)

# data frame with fitted and actual returns
fitted <- model$fitted
actual <- Combined$Stock_CZ
returns <- data.frame(fitted, actual)

# create a scatter plot of fitted vs. actual values
plot_ly(returns, x = ~fitted, y = ~actual, type = 'scatter', mode = 'markers')

################################ Linear Regression Ends ################################





