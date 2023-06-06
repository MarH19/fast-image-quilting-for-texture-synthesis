import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('ggplot')

# Read the CSV file into a dataframe
base = pd.read_csv('../base/benchmark_blocksize_basic.csv')
opt1 = pd.read_csv('../opt1/benchmark_blocksize_opt1.csv')
opt2 = pd.read_csv('../opt2/benchmark_blocksize_opt2.csv')
opt3 = pd.read_csv('../opt3/benchmark_blocksize_opt3.csv')

# Plot the data
#plt.plot(df['seconds'], df['flops'])
#plt.xlabel('Seconds')
#plt.ylabel('Flops')
#plt.title('Flops vs Seconds')
#plt.show()





# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(7)]


ax.plot(n,base['seconds'],marker='s', linewidth=2, color="mediumblue", label='Base implementation')

ax.plot(n,opt1['seconds'],marker='^', linewidth=2, color="dimgray", label='Optimization 1')

ax.plot(n,opt2['seconds'], marker='o', linewidth=2, color="darkred", label='Optimization 2')

ax.plot(n,opt3['seconds'], marker='8', linewidth=2, color="darkgreen", label='Optimization 3')



plt.xticks(ticks=n, rotation=0)
ax.set_xlabel('(blocksize, overlap) with fixed input size of 239 pixels', fontsize=14, labelpad=15)
xticks = list(zip(base['blocksize'], base['overlap']))
xtick_labels = [f'({x[0]}, {x[1]})' for x in xticks]
plt.xticks(range(len(xticks)), xtick_labels)

ax.set_ylabel('Runtime [s]', fontsize=14, labelpad=15)


#plt.legend(facecolor="white", fontsize=12, loc='center right')

ax.set_title('Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz', fontsize=16)



ax.grid(axis='x')
plt.yscale('log') 
plt.show()




