
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# # Read CSV file
# df = pd.read_csv(r'/Users/dkp116/Desktop/Master York/Diss/code/Intergral Check/build/output.csv')

# # Configure plot style
# sns.set_style("whitegrid")
# plt.figure(figsize=(12, 6))

# # Define barrier level
# barrier_level = 80  # Adjust this value as needed

# # Split data at the jump points (time = 0.1 and 0.8 years)
# jump_time1 = 0.1
# jump_time2 = 0.8

# # Split data into three segments
# df_part1 = df[df['Time (Years)'] <= jump_time1]
# df_part2 = df[(df['Time (Years)'] > jump_time1) & (df['Time (Years)'] <= jump_time2)]
# df_part3 = df[df['Time (Years)'] > jump_time2]

# # Plot stock price trajectory with discontinuities
# plt.plot(df_part1['Time (Years)'], df_part1['Stock Price'], 
#          color='royalblue', linewidth=1.5, linestyle='-', label='Stock Price')
# plt.plot(df_part2['Time (Years)'], df_part2['Stock Price'], 
#          color='royalblue', linewidth=1.5, linestyle='-')
# plt.plot(df_part3['Time (Years)'], df_part3['Stock Price'], 
#          color='royalblue', linewidth=1.5, linestyle='-')

# # Get prices at the jump discontinuities
# # First jump (0.1)
# last_price_part1 = df_part1['Stock Price'].iloc[-1]
# first_price_part2 = df_part2['Stock Price'].iloc[0]

# # Second jump (0.8)
# last_price_part2 = df_part2['Stock Price'].iloc[-1]
# first_price_part3 = df_part3['Stock Price'].iloc[0]

# # Add vertical dotted lines at both jumps
# plt.plot([jump_time1, jump_time1], [last_price_part1, first_price_part2],
#          linestyle=':', color='blue', linewidth=1.5, label='Jump')
# plt.plot([jump_time2, jump_time2], [last_price_part2, first_price_part3],
#          linestyle=':', color='blue', linewidth=1.5)

# # Add barrier level line
# plt.axhline(y=barrier_level, color='red', linestyle='--', linewidth=1.5, label='Barrier Level')

# # Add initial price line (optional)
# plt.plot(df['Time (Years)'], [df['Stock Price'].iloc[0]]*len(df), 
#          '--', color='gray', label='Initial Price')

# # Format plot
# plt.title("Merton Jump Diffussion Stock Simulation with Crossing of Barrier in Brownian Bridge", fontsize=14)
# plt.xlabel("Time", fontsize=12)
# plt.ylabel("Stock Price", fontsize=12)
# plt.xlim(0, df['Time (Years)'].max())
# plt.legend()

# # Save and show plot
# plt.tight_layout()
# plt.savefig('stock_simulation_jump.png', dpi=300)
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read CSV file
df = pd.read_csv(r'/Users/dkp116/Desktop/Master York/Diss/code/Intergral Check/build/output.csv')

# Set style and context
sns.set_context("talk")  # Larger font for presentations/dissertations
sns.set_style("whitegrid")
palette = sns.color_palette("deep")

plt.figure(figsize=(14, 7))

# Parameters
barrier_level = 80
jump_time1 = 0.1
jump_time2 = 0.6

# Split data
df_part1 = df[df['Time (Years)'] <= jump_time1]
df_part2 = df[(df['Time (Years)'] > jump_time1) & (df['Time (Years)'] <= jump_time2)]
df_part3 = df[df['Time (Years)'] > jump_time2]

# Plot stock trajectory
plt.plot(df_part1['Time (Years)'], df_part1['Stock Price'], 
         color=palette[0], linewidth=2.2, label='Stock Price')
plt.plot(df_part2['Time (Years)'], df_part2['Stock Price'], 
         color=palette[0], linewidth=2.2)
plt.plot(df_part3['Time (Years)'], df_part3['Stock Price'], 
         color=palette[0], linewidth=2.2)

# Jump prices
last_price_part1 = df_part1['Stock Price'].iloc[-1]
first_price_part2 = df_part2['Stock Price'].iloc[0]
last_price_part2 = df_part2['Stock Price'].iloc[-1]
first_price_part3 = df_part3['Stock Price'].iloc[0]

# Jump lines and markers
# Jump at 0.1
plt.plot([jump_time1, jump_time1], [last_price_part1, first_price_part2],
         linestyle=':', color=palette[1], linewidth=2, label='Jump Discontinuity')
plt.scatter(jump_time1, last_price_part1, s=30, facecolors='white', edgecolors=palette[1], linewidths=2, zorder=5)
plt.scatter(jump_time1, first_price_part2, s=30, color=palette[1], zorder=5)

# Jump at 0.8
plt.plot([jump_time2, jump_time2], [last_price_part2, first_price_part3],
         linestyle=':', color=palette[1], linewidth=2)
plt.scatter(jump_time2, last_price_part2, s=30, facecolors='white', edgecolors=palette[1], linewidths=2, zorder=5)
plt.scatter(jump_time2, first_price_part3, s=30, color=palette[1], zorder=5)

# Barrier line in red
plt.axhline(y=barrier_level, color='red', linestyle='--', linewidth=2, label='Barrier Level')

# Initial price line
initial_price = df['Stock Price'].iloc[0]
plt.axhline(y=initial_price, color='gray', linestyle='--', linewidth=1.5, label='Initial Price')

# Title and axis labels
plt.title("Stock Path under Merton Jump Diffusion\nwith No Crossing", 
          fontsize=16, weight='bold')
plt.xlabel("Time", fontsize=14)
plt.ylabel("Stock Price", fontsize=14)
plt.xlim(0, df['Time (Years)'].max())

# Legend outside the plot
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., fontsize=12)

# Save and show
plt.tight_layout()
plt.savefig('stock_simulation_jump.png', dpi=300, bbox_inches='tight')
plt.show()
