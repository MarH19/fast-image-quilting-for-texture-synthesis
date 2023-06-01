import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('ggplot')

# Read the CSV file into a dataframe
base = pd.read_csv('base/benchmark_numblocks_basic.csv')
opt1 = pd.read_csv('opt1/benchmark_numblocks_opt1.csv')
opt2 = pd.read_csv('opt2/benchmark_numblocks_opt2.csv')
opt3 = pd.read_csv('opt3/benchmark_numblocks_opt3.csv')

# Plot the data
#plt.plot(df['seconds'], df['flops'])
#plt.xlabel('Seconds')
#plt.ylabel('Flops')
#plt.title('Flops vs Seconds')
#plt.show()





# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(6)]

ax.plot(opt2['numblocks'],opt2['flops']/opt2['cycles'], marker='o', linewidth=2, color="darkred", label='Optimization 2')

ax.plot(opt3['numblocks'],opt3['flops']/opt3['cycles'], marker='8', linewidth=2, color="darkgreen", label='Optimization 3')



plt.xticks(ticks=base['numblocks'])
ax.set_xlabel('Number of blocks in the output', fontsize=14, labelpad=15)


ax.set_ylabel('Performance [intops/cycles]', fontsize=14, labelpad=15)#, rotation=0, loc='top')


plt.legend(facecolor="white", fontsize=12)

ax.set_title('Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz', fontsize=16)

ax.set_ylim(0, None)


ax.grid(axis='x')
#plt.yscale('log') 
plt.show()



