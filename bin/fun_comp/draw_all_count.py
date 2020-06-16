import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("all_count.txt", header=None,index_col=0) 
df.columns=['pir', 'hyb','clan']

print(df)

plt.figure()
df.iloc[0:6,:].plot(kind='bar',colormap='cividis')
plt.savefig("human_count.png")
plt.figure()
df.iloc[6:,:].plot(kind='bar',colormap='cividis')
plt.savefig("elegan_count.png")
