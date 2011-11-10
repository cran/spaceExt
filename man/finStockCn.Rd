\name{finStocksCn}
          \docType{data}
          \alias{finStocksCn}
          \title{Continuously Compounded Daily Returns for Financial Sector of Chinese Stock Market}
          \description{
            This data set contains the continuously compounded daily returns for financial sector of Chinese stock market.
			The time ranges from January 4th, 2011 to July 27th, 2011. It totally consists of 34 stocks: 3 insurance stocks, 
			15 security stocks and 16 bank stocks.
          }
          \usage{finStocksCn}
          \format{A list containing two components
		  \enumerate{
		  \item $names, a vector of length 34 containing the stock names. 
		  \item $returns, a 155-by-34 matrix, each column corresponding to a stock and each row to a trading day. 
		  The stock name of its i(th) column is contained in the i(th) component of finStocksCn$names.
		  }
		  }
\keyword{datasets}