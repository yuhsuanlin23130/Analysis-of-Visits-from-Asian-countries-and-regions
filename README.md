#### 說明
對亞洲部分國家地區(日本、香港、澳門、馬來西亞、新加坡、印度)來台人次，採用以下三種方法來建構時間序列分析模型： 1. Smoothing by Centered Moving Average 2. Smoothing by Linear Regression Model 3. Regression Model by Indicator Variables，並對於各國家地區選擇最佳模型來預測其未來一年(民國108年5月到109年4月)每月來台人次，並分析其季節性指數。 
資料來源為<a href="https://stat.taiwan.net.tw/inboundSearch">交通部觀光局觀光統計資料庫</a>。
 
#### 執行環境
R + Rstudio

#### 目錄結構說明
* model資料夾底下包含所分析國家地區之原始資料(\*.csv)與模型建構及分析(\*.Rmd, \.nb)
* result.xlsx為各國家地區之預測數據、季節性指標、Error Metrics結果彙整。
