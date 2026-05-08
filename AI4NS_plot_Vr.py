import matplotlib.pyplot as plt
import numpy as np

# Font size settings
LABEL_SIZE = 22
TICK_SIZE = 16

fig, axes = plt.subplots(2, 6, figsize=(24, 10), sharex=True, sharey=True)
axes = axes.flatten()

# Define labels
labels = ['DFT', 'FF']

for i in range(1, 13):
    filename = f"split_{i}.txt"
    ax = axes[i-1]
    
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]
        
        # Plot lines - passing only the specific string label for each
        line1, = ax.plot(x, data[:, 1], color='blue', lw=2.5, label=labels)
        line3, = ax.plot(x, data[:, 3], color='red', lw=2.5, label=labels)
        
        # Set X-axis range and specific ticks
        ax.set_xlim(2, 13)
        ax.set_xticks([2, 4, 6, 8, 10, 12])
        ax.set_ylim(-2, 2)
        
        ax.tick_params(axis='both', which='major', labelsize=TICK_SIZE)
        ax.grid(True, linestyle='--', alpha=0.4)
        
    except (FileNotFoundError, IOError):
        continue

# Add a single legend at the top
fig.legend([line1, line3], labels, loc='upper center', 
           bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=LABEL_SIZE, frameon=False)

# Common labels
fig.text(0.5, 0.02, 'Distance ($\AA$)', ha='center', fontsize=LABEL_SIZE)
fig.text(0.01, 0.5, 'Energy (kcal/mol)', va='center', rotation='vertical', fontsize=LABEL_SIZE)

plt.tight_layout(rect=[0.03, 0.05, 1, 0.95])

# Save as PDF
plt.savefig("my_plots.pdf", format='pdf', bbox_inches='tight', dpi=300)
print("Figure saved successfully as my_plots.pdf")

