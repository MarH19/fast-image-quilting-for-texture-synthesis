import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

fig, ax = plt.subplots(figsize=[10, 7])


opt3_O3 = pd.read_csv('benchmark_numblocks_opt3_O3.csv')
opt2_O3 = pd.read_csv('benchmark_numblocks_opt2_O3.csv')





# Define peak performance and memory bandwidth
peak_perf = 2  # flops/cycle
bandwidth = 18.57  # bytes/cycle

# Create data points for roofline plot
xs = np.logspace(-5, 16, num=1000, base=2)  # operational intensity (FLOPs/byte)

ys = np.minimum(xs * bandwidth, peak_perf)    # Scalar performance peak
ys2 = np.minimum(xs * bandwidth, 48) # Vetor performance peak


def nr_bytes(numblocks, bs, ov):
    input_bytes = (239**2) * 3 * 2
    errormatrix_bytes = ((239-bs+1)**2)*4
    output_bytes = ((numblocks * bs - (bs-1) * ov)**2) * 3 * 2
    errormatrix_dp_bytes = (ov*bs)*4
    return input_bytes + errormatrix_bytes + output_bytes


def nr_bytes2(numblocks, bs, ov):
    input_bytes = (239**2) * 3 * 4
    errormatrix_bytes = ((239-bs+1)**2)*4
    output_bytes = ((numblocks * bs - (bs-1) * ov)**2) * 3 * 4
    errormatrix_dp_bytes = (ov*bs)*4
    return input_bytes + errormatrix_bytes + output_bytes


# Plot the roofline and data points
plt.xscale('log', base=2)
plt.yscale('log', base=2)
ax.plot(xs, ys, 'k-', linewidth=3, markersize=12)
ax.plot(xs, ys2, 'k--', linewidth=3, markersize=12)
ax.plot(opt3_O3['flops']/opt3_O3['numblocks'].apply(nr_bytes,args=(48,8)), opt3_O3['flops']/opt3_O3['cycles'], marker='o', linewidth=1.5, markersize=8)
ax.plot(opt2_O3['flops']/opt2_O3['numblocks'].apply(nr_bytes2,args=(48,8)), opt2_O3['flops']/opt2_O3['cycles'], marker='o', color="blue", linewidth=1.5, markersize=8)


plt.xlabel('Operational Intensity [intops/byte]', fontsize=19, color='black')

ax.text(0.0, 1.02, 'Performance [intops/cycle]', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize=19)



ax.set_title('Roofline Plot', fontsize=22, color='black', pad=30)



# Set the x-ticks to be evenly spaced on a log scale
xticks = [2**i for i in range(-5, 16, 2)]
plt.xticks(xticks, [str(x) for x in xticks], fontsize=14, color='black')



# Set the y-ticks to be evenly spaced on a log scale
yticks = [2**i for i in range(-3, 6)]
plt.yticks(yticks, [str(y) for y in yticks], fontsize=14, color='black')



#ax.spines['bottom'].set_color('black')
#ax.spines['bottom'].set_linewidth(1)
ax.grid(axis='x', linewidth=1)
ax.grid(axis='y', linewidth=1)
ax.tick_params(axis='both', width=1, color='black')

plt.savefig('roofline.pdf', bbox_inches='tight')

plt.show()
