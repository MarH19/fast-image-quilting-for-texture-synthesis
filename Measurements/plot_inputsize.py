import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('ggplot')

# Read the CSV file into a dataframe
base = pd.read_csv('base/benchmark_inputsize_basic.csv')
opt1 = pd.read_csv('opt1/benchmark_inputsize_opt1.csv')
opt2 = pd.read_csv('opt2/benchmark_inputsize_opt2.csv')
opt3 = pd.read_csv('opt3/benchmark_inputsize_opt3.csv')

# Plot the data
#plt.plot(df['seconds'], df['flops'])
#plt.xlabel('Seconds')
#plt.ylabel('Flops')
#plt.title('Flops vs Seconds')
#plt.show()





# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(12)]


ax.plot(base['inputsize'],base['seconds'],marker='s', linewidth=2, color="mediumblue", label='Base implementation')

ax.plot(opt1['inputsize'],opt1['seconds'],marker='^', linewidth=2, color="dimgray", label='Optimization 1')

ax.plot(opt2['inputsize'],opt2['seconds'], marker='o', linewidth=2, color="darkred", label='Optimization 2')

ax.plot(opt3['inputsize'],opt3['seconds'], marker='8', linewidth=2, color="darkgreen", label='Optimization 3')



ax.set_xlabel('Input image size', fontsize=14, labelpad=15)
xticks = [50, 100, 150, 200, 250]
plt.xticks(xticks)
plt.xlim(50, None) 


ax.set_ylabel('Runtime [s]', fontsize=14, labelpad=15)#, rotation=0, loc='top')


plt.legend(facecolor="white", fontsize=12)

ax.set_title('AMD Ryzen 7 5800U @ 1.90GHz', fontsize=16)


ax.grid(axis='x')
plt.yscale('log') 
plt.show()








