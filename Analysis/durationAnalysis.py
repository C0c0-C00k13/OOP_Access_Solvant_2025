import matplotlib.pyplot as plt

# Data lists
duration_2c8r = [3.3847103118896484, 2.61596417427063, 2.3511838912963867, 2.209131956100464, 2.2690656185150146]
duration_2oe4 = [25.170685052871704, 27.05111598968506, 25.61398482322693, 25.91768741607666, 27.69119167327881]
duration_1c26 = [1.4570181369781494, 1.2581472396850586, 1.2174153327941895, 1.3090920448303223, 1.318755865097046]
duration_1bj5 = [340.2541913986206, 618.5669720172882, 695.0584111213684, 758.2242970466614, 654.4952387809753]
duration_6pwf = [4965.598150968552, 3889.2021148204803, 4451.032356023788, 3647.5049815177917, 4869.671242713928]

lists = [duration_2c8r, duration_2oe4, duration_1c26, duration_1bj5, duration_6pwf]
colors = ['red', 'blue', 'green', 'purple', 'orange']
labels = ['2c8r', '2oe4', '1c26', '1bj5', '6pwf']

# X-axis positions (indices)
positions = list(range(1,6))

# Second X-axis reference (linked number for each list)
reference_x = [10, 20, 30, 40, 50]  # You can change these values to reflect your context

# --- Plot 1: Values vs Index ---
plt.figure(figsize=(10, 5))
for i, lst in enumerate(lists):
    plt.scatter(positions, lst, label=labels[i], color=colors[i])
plt.title('Scatter Plot: Values vs Position')
plt.xlabel('Number of Run')
plt.ylabel('Time (s)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 2: Values vs Associated Number ---
# Create boxplot
plt.figure(figsize=(10, 6))
plt.boxplot(lists, patch_artist=True,
            boxprops=dict(facecolor='lightblue', color='blue'),
            medianprops=dict(color='red'),
            whiskerprops=dict(color='blue'),
            capprops=dict(color='blue'),
            flierprops=dict(markerfacecolor='orange', marker='o', markersize=5))

plt.title("Duration of the calculation of ASA of the 5 proteins")
plt.xlabel("Proteins")
plt.ylabel("Duration (s)")
plt.xticks([1, 2, 3, 4, 5],labels)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
