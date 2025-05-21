# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# # Read CSV file
# df = pd.read_csv(r'/Users/dkp116/Desktop/Master York/Diss/code/Intergral Check/build/output.csv')

# # Configure plot style
# sns.set_style("whitegrid")
# plt.figure(figsize=(12, 6))

# # Plot stock price trajectory
# plt.plot(df['Time (Years)'], df['Stock Price'], 
#          color='royalblue', linewidth=1.5)

# # Format plot
# plt.title("Stock Price Simulation (Geometric Brownian Motion)", fontsize=14)
# plt.xlabel("Time (Years)", fontsize=12)
# plt.ylabel("Stock Price ($)", fontsize=12)
# plt.xlim(0, df['Time (Years)'].max())
# plt.gca().xaxis.set_major_locator(plt.MaxNLocator(10))  # Clean x-axis ticks

# # Add technical indicators (optional)
# plt.plot(df['Time (Years)'], [df['Stock Price'].iloc[0]]*len(df), 
#          '--', color='gray', label='Initial Price')
# plt.legend()

# # Save and show plot
# plt.tight_layout()
# plt.savefig('stock_simulation.png', dpi=300)
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read CSV file
df = pd.read_csv(r'/Users/dkp116/Desktop/Master York/Diss/code/Intergral Check/build/output.csv')

# Configure plot style
sns.set_style("whitegrid")
plt.figure(figsize=(12, 6))

# Define barrier level
barrier_level = 80  # Adjust this value as needed

# Split data at the jump point (time = 0.1 years)
jump_time = 0.1  # Matches the C++ loop condition (100*T where T=1/1000)
df_before = df[df['Time (Years)'] <= jump_time]
df_after = df[df['Time (Years)'] > jump_time]

# Get prices at the jump discontinuity
last_price_before = df_before['Stock Price'].iloc[-1]
first_price_after = df_after['Stock Price'].iloc[0]

# Plot stock price trajectory with discontinuity
plt.plot(df_before['Time (Years)'], df_before['Stock Price'], 
         color='royalblue', linewidth=1.5, linestyle='-', label='Stock Price')
plt.plot(df_after['Time (Years)'], df_after['Stock Price'], 
         color='royalblue', linewidth=1.5, linestyle='-')

# Add vertical dotted line at the jump discontinuity
plt.plot([jump_time, jump_time], [last_price_before, first_price_after],
         linestyle=':', color='royalblue', linewidth=1.5)

# Add barrier level line
plt.axhline(y=barrier_level, color='red', linestyle='--', linewidth=1.5, label='Barrier Level')

# Add initial price line (optional)
plt.plot(df['Time (Years)'], [df['Stock Price'].iloc[0]]*len(df), 
         '--', color='gray', label='Initial Price')

# Format plot
plt.title("Stock Price Simulation with Jump Discontinuity", fontsize=14)
plt.xlabel("Time (Years)", fontsize=12)
plt.ylabel("Stock Price ($)", fontsize=12)
plt.xlim(0, df['Time (Years)'].max())
plt.legend()

# Save and show plot
plt.tight_layout()
plt.savefig('stock_simulation_jump.png', dpi=300)
plt.show()