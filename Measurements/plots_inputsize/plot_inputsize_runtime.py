import matplotlib.pyplot as plt
import pandas as pd
import os

plt.style.use('ggplot')

# Base implementation
base_O3_novec = pd.read_csv('../base/benchmark_inputsize_basic_O3_novec.csv')
base_O3 = pd.read_csv('../base/benchmark_inputsize_basic_O3.csv')
base_O3_fastmath = pd.read_csv('../base/benchmark_inputsize_basic_O3_fastmath.csv')

# Optimization 1
opt1_O3_novec = pd.read_csv('../opt1/benchmark_inputsize_opt1_O3_novec.csv')
opt1_O3 = pd.read_csv('../opt1/benchmark_inputsize_opt1_O3.csv')
opt1_O3_fastmath = pd.read_csv('../opt1/benchmark_inputsize_opt1_O3_fastmath.csv')

# Optimization 2
opt2_O3_novec = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3_novec.csv')
opt2_O3 = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3.csv')
opt2_O3_fastmath = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3_fastmath.csv')

# Optimization 3
opt3_O3_novec = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3_novec.csv')
opt3_O3 = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3.csv')
opt3_O3_fastmath = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3_fastmath.csv')




# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(12)]

#Base implementation
#ax.plot(base_O3_novec['inputsize'],base_O3_novec['seconds'],marker='s', markersize=12, linewidth=4, color="#87CEEB", label='Base -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(base_O3['inputsize'],base_O3['seconds'],marker='s', markersize=12, linewidth=4, color="#1E90FF", label='Base -O3 -mfma')
#ax.plot(base_O3_fastmath['inputsize'],base_O3_fastmath['seconds'],marker='s', markersize=12, linewidth=4, color="#0000FF", label='Base -O3 -mfma -ffast-math -march=native')

#Optimization 1
#ax.plot(opt1_O3_novec['inputsize'],opt1_O3_novec['seconds'],marker='s', markersize=12, linewidth=4, color="#FFA07A", label='Opt1 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(opt1_O3['inputsize'],opt1_O3['seconds'],marker='s', markersize=12, linewidth=4, color="#FF4500", label='Opt1 -O3 -mfma')
#ax.plot(opt1_O3_fastmath['inputsize'],opt1_O3_fastmath['seconds'],marker='s', markersize=12, linewidth=4, color="#FF0000", label='Opt1 -O3 -mfma -ffast-math -march=native')

#Optimization 2
#ax.plot(opt2_O3_novec['inputsize'],opt2_O3_novec['seconds'],marker='s', markersize=12, linewidth=4, color="#7CFC00", label='Opt2 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(opt2_O3['inputsize'],opt2_O3['seconds'],marker='s', markersize=12, linewidth=4, color="#00A36C", label='Opt2 -O3 -mfma')
#ax.plot(opt2_O3_fastmath['inputsize'],opt2_O3_fastmath['seconds'],marker='s', markersize=12, linewidth=4, color="#355E3B", label='Opt2 -O3 -mfma -ffast-math -march=native')

#Optimization 3
#ax.plot(opt3_O3_novec['inputsize'],opt3_O3_novec['seconds'],marker='s', markersize=12, linewidth=4, color="#DA70D6", label='Opt3 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(opt3_O3['inputsize'],opt3_O3['seconds'],marker='s', markersize=12, linewidth=4, color="#8A2BE2", label='Opt3 -O3 -mfma')
#ax.plot(opt3_O3_fastmath['inputsize'],opt3_O3_fastmath['seconds'],marker='s', markersize=12, linewidth=4, color="#800080", label='Opt3 -O3 -mfma -ffast-math -march=native')





ax.set_xlabel('Input image size', fontsize=20, color='black')
xticks = [50, 100, 150, 200, 250]
plt.xticks(xticks, fontsize=14, color='black')
plt.xlim(50, None) 


plt.yticks(fontsize=14, color='black')
ax.text(0.0, 1.02, 'Runtime [s]', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize=19)


ax.set_title('Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz', fontsize=22, color='black', pad=30)


ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_linewidth(2)


plt.legend()



ax.grid(axis='x', linewidth=0)
ax.grid(axis='y', linewidth=2)

ax.tick_params(axis='both', width=2, color='black')
ax.tick_params(axis='both', which='both', width=2, color='black')

#plt.yscale('log') 


# store the plot
directory = 'generated_plots' 

if not os.path.exists(directory):
    os.makedirs(directory)

plt.savefig('generated_plots/plot_inputsize_runtime.pdf', bbox_inches='tight')
plt.show()








