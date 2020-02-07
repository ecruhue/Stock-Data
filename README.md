# Stock-Data
Stock price data is diverse and situational, and it is unlikely that any single model will be uniformly best across all industries or contexts. Based on our results, ARIMA GARCH methods are better for Consumer Discretionary and Financial industries, and LSTM models are better for Healthcare and Industrials. Specifically, we found Cumulative Year (CumYr) ARIMA GARCH performs best for the Consumer Discretionary industry, year by year (YrByYr) ARIMA GARCH performs best for Financials, YrByYr multivariate LSTM performs best for Healthcare, and YrByYr univariate LSTM performs best for Industrials. Overall, the LSTM models with YrByYr under multivariate condition perform better than LSTM models under other conditions. 
