import pandas as pd
import sys


df = pd.read_csv("output/gradientSeries.csv")
df['modS'] = df['mod'].shift(1)
df['phaseS'] = df['phase'].shift(1)

df = df.drop([i for i in range(3000)])

print(df.shape[0])
c = [i for i in range(0,df.shape[0])]

df = df.assign(code=c)

df.to_csv("gradientSeries_Shift.csv",index=False)
