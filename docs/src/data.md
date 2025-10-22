# Data
We have included a market dataset that we use for examples and activities in the equity domain. This dataset holds the daily open, high, low, close, and volume data for a selection of stocks between 2014 and 2024. 


```@docs
VLQuantitativeFinancePackage.MyTrainingMarketDataSet
VLQuantitativeFinancePackage.MyTestingMarketDataSet
```

## Options Data
We have also included an options chain dataset that we use for examples and activities in the options domain. This dataset holds a snapshot of options data for a selection of stocks on a specific date. Currently, we only provide data for AMD options as of October 22, 2025. The dataset includes information such as strike prices, expiration dates, bid and ask prices, and implied volatilities.

```@docs
VLQuantitativeFinancePackage.MyOptionsChainDataSet
```