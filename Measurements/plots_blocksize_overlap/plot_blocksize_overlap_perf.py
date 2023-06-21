import matplotlib.pyplot as plt
import pandas as pd
import os

plt.style.use('ggplot')


# Base implementation
base_O3_novec = pd.read_csv('../base/benchmark_blocksize_basic_O3_novec.csv')
#base_O3 = pd.read_csv('../base/benchmark_blocksize_basic_O3.csv')
#base_O3_fastmath = pd.read_csv('../base/benchmark_blocksize_basic_O3_fastmath.csv')



# Optimization 2
#opt2_O3_novec = pd.read_csv('../opt2/benchmark_blocksize_opt2_O3_novec.csv')
opt2_O3 = pd.read_csv('../opt2/benchmark_blocksize_opt2_O3.csv')
#opt2_O3_fastmath = pd.read_csv('../opt2/benchmark_blocksize_opt2_O3_fastmath.csv')

# Optimization 3
#opt3_O3_novec = pd.read_csv('../opt3/benchmark_blocksize_opt3_O3_novec.csv')
opt3_O3 = pd.read_csv('../opt3/benchmark_blocksize_opt3_O3.csv')
#opt3_O3_fastmath = pd.read_csv('../opt3/benchmark_blocksize_opt3_O3_fastmath.csv')





# Creating figure and axis objects using subplots()
fig, ax = plt.subplots(figsize=[10, 7])

n = [i for i in range(7)]

#Optimization 3
#ax.plot(n,opt3_O3_novec['flops']/opt3_O3_novec['cycles'],marker='s', markersize=12, linewidth=4, color="#DA70D6", label='Opt3 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(n,opt3_O3['flops']/opt3_O3['cycles'],marker='s', markersize=12, linewidth=4, color="#8A2BE2", label='Opt3 -O3 -mfma')
#ax.plot(n,opt3_O3_fastmath['flops']/opt3_O3_fastmath['cycles'],marker='s', markersize=12, linewidth=4, color="#800080", label='Opt3 -O3 -mfma -ffast-math -march=native')

#Optimization 2
#ax.plot(n,opt2_O3_novec['flops']/opt2_O3_novec['cycles'],marker='s', markersize=12, linewidth=4, color="#7CFC00", label='Opt2 -O3 -mfma -fno-tree-vectorize -fno-slp-vectorize')
ax.plot(n,opt2_O3['flops']/opt2_O3['cycles'],marker='s', markersize=12, linewidth=4, color="#00A36C", label='Opt2 -O3 -mfma')
#ax.plot(n,opt2_O3_fastmath['flops']/opt2_O3_fastmath['cycles'],marker='s', markersize=12, linewidth=4, color="#355E3B", label='Opt2 -O3 -mfma -ffast-math -march=native')


plt.xticks(ticks=n, rotation=0)
ax.set_xlabel('(blocksize, overlap) with fixed input width of 239 pixels', fontsize=20, color='black')
xticks = list(zip(base_O3_novec['blocksize'], base_O3_novec['overlap']))
xtick_labels = [f'({x[0]}, {x[1]})' for x in xticks]
plt.xticks(range(len(xticks)), xtick_labels, fontsize=14, color='black')



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

plt.savefig('generated_plots/plot_blocksize_overlap_perf.pdf', bbox_inches='tight')
plt.show()





