import matplotlib.pyplot as plt
import pandas as pd
import os

plt.style.use('ggplot')

# Optimization 2
#opt2_O3_novec = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3_novec.csv')
opt2_O3 = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3.csv')
#opt2_O3_fastmath = pd.read_csv('../opt2/benchmark_inputsize_opt2_O3_fastmath.csv')

# Optimization 3
#opt3_O3_novec = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3_novec.csv')
opt3_O3 = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3.csv')
#opt3_O3_fastmath = pd.read_csv('../opt3/benchmark_inputsize_opt3_O3_fastmath.csv')





# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(12)]

#Optimization 3
#ax.plot(opt3_O3_novec['inputsize'],opt3_O3_novec['flops']/opt3_O3_novec['cycles'],marker='s', markersize=12, linewidth=4, color="#DA70D6", label='Opt3 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(opt3_O3['inputsize'],opt3_O3['flops']/opt3_O3['cycles'],marker='s', markersize=12, linewidth=4, color="#8A2BE2", label='Opt3 -O3 -mfma')
#ax.plot(opt3_O3_fastmath['inputsize'],opt3_O3_fastmath['flops']/opt3_O3_fastmath['cycles'],marker='s', markersize=12, linewidth=4, color="#800080", label='Opt3 -O3 -mfma -ffast-math -march=native')

#Optimization 2
#ax.plot(opt2_O3_novec['inputsize'],opt2_O3_novec['flops']/opt2_O3_novec['cycles'],marker='s', markersize=12, linewidth=4, color="#7CFC00", label='Opt2 -O3 -mfma -fno-tree-vectorize')
ax.plot(opt2_O3['inputsize'],opt2_O3['flops']/opt2_O3['cycles'],marker='s', markersize=12, linewidth=4, color="#00A36C", label='Opt2 -O3 -mfma')
#ax.plot(opt2_O3_fastmath['inputsize'],opt2_O3_fastmath['flops']/opt2_O3_fastmath['cycles'],marker='s', markersize=12, linewidth=4, color="#355E3B", label='Opt2 -O3 -mfma -ffast-math -march=native')





ax.set_xlabel('Input image size', fontsize=20, color='black')
xticks = [50, 100, 150, 200, 250]
plt.xticks(xticks, fontsize=14, color='black')
plt.xlim(50, None) 



plt.yticks(fontsize=14, color='black')
ax.text(0.0, 1.02, 'Performance [intops/cycles]', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize=19)
ax.set_ylim(0, None)



ax.set_title('Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz', fontsize=22, color='black', pad=30)


plt.legend()



ax.spines['bottom'].set_color('black')
ax.spines['bottom'].set_linewidth(2)

ax.grid(axis='x', linewidth=0)
ax.grid(axis='y', linewidth=2)

ax.tick_params(axis='both', width=2, color='black')




# store the plot
directory = 'generated_plots' 

if not os.path.exists(directory):
    os.makedirs(directory)

plt.savefig('generated_plots/plot_inputsize_perf.pdf', bbox_inches='tight')
plt.show()







