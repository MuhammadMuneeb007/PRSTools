import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Sample data

df = pd.read_csv("ExecutionTime.csv",sep=",")
df['Total_Time_Normalized'] = (df['Total_Time'] - df['Total_Time'].min()) / (df['Total_Time'].max() - df['Total_Time'].min())
df['Total_Memory_Normalized'] = (df['Total_Memory'] - df['Total_Memory'].min()) / (df['Total_Memory'].max() - df['Total_Memory'].min())

# Set seaborn style
sns.set(style="whitegrid")

# Create figure and axes
fig, ax = plt.subplots(figsize=(12, 6))

# Create bar width and x locations
bar_width = 0.4
x = np.arange(len(df['Method']))

# Plot bars for time and memory (normalized)
ax.bar(x - bar_width/2, df['Total_Time_Normalized'], width=bar_width, color='blue', label='Normalized Time')
ax.bar(x + bar_width/2, df['Total_Memory_Normalized'], width=bar_width, color='orange', label='Normalized Memory')

# Set x-axis labels with vertical orientation
ax.set_xticks(x)
ax.set_xticklabels(df['Method'], rotation=90, ha="center")

# Set y-axis to log scale
ax.set_yscale('log')

# Add labels and title
ax.set_xlabel('Method')
ax.set_ylabel('Normalized Value (Log Scale)')
ax.set_title('Comparison of Normalized Time and Memory Usage for Different Methods')

# Add legend
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

exit(0)


def convert_time_unit(df, time_col):
    # Convert to minutes
    df[time_col] = df[time_col] / 60
    unit = "minutes"
    
    # If all values are >= 1, convert to hours
    if df[time_col].min() >= 1:
        df[time_col] = df[time_col] / 60
        unit = "hours"
        
        # If all values are >= 1, convert to days
        if df[time_col].min() >= 1:
            df[time_col] = df[time_col] / 24
            unit = "days"
    
    return df, unit

def convert_memory_unit(df, memory_col):
    # Convert to MB
    df[memory_col] = df[memory_col] / 1024
    unit = "MB"
    print(df)
    # If all values are >= 1 MB, convert to GB
    if df[memory_col].min() >= 1:
        df[memory_col] = df[memory_col] / 1024
        unit = "GB"
        print(df)
        # If all values are >= 1 GB, convert to TB
        #if df[memory_col].min() >= 1:
        #    df[memory_col] = df[memory_col] / 1024
        #    unit = "TB"
          
    return df, unit

 

def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else: 
        alignment = "left"
    return rotation, alignment

def add_labels(angles, values, labels, offset, ax):
    
    # This is the space between the end of the bar and the label
    padding = 4
    
    # Iterate over angles, values, and labels, to add all of them.
    for angle, value, label, in zip(angles, values, labels):
        angle = angle
        
        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle, 
            y=value + padding, 
            s=label, 
            ha=alignment, 
            va="center", 
            rotation=rotation, 
            rotation_mode="anchor",size=5
        ) 


results = pd.read_csv("ExecutionTime.csv",sep=",")
print(results.head())

# Convert Total_Time to the appropriate unit
results, time_unit = convert_time_unit(results, 'Total_Time')
print(time_unit)

# Convert Total_Memory to the appropriate unit
results, memory_unit = convert_memory_unit(results, 'Total_Memory')
print(memory_unit)

def groupbasedonalgorithm(results,col):
    results = results.sort_values(by=['Phenotype','Method'])
    print(results)
    OFFSET = np.pi / 2
    ANGLES = np.linspace(0, 2 * np.pi, len(results), endpoint=False)
    
    VALUES = results[col].values
    
    results[col+"A"] = results[col]/1
    results[col+"A"] = results[col+"A"].round(2)
    
    #LABELS = results["Total_TimeA"].astype(str) +" : "+ results["Phenotype"]
    if "Time" in col:
        LABELS = results[col+"A"].astype(str) +" : "+ time_unit
    else:
        LABELS = results[col+"A"].astype(str) +" : "+ memory_unit

    LABELS = LABELS.values
    GROUP = results["Method"].values

    grouped = results.groupby('Method')
    print(grouped)

    GROUPS_SIZE = []
    unique = []
    for name, group in grouped:
        GROUPS_SIZE.append(len(group))
        unique.append(name)



    PAD = 3
    ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
    ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
    WIDTH = (2 * np.pi) / len(ANGLES)
    offset = 0
    IDXS = []
    for size in GROUPS_SIZE:
        IDXS += list(range(offset + PAD, offset + size + PAD))
        offset += size + PAD
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "polar"})
    ax.set_theta_offset(OFFSET)
    ax.set_ylim(-100, 100)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    COLORS = [f"C{i}" for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]
    print(VALUES)
    print(ANGLES[IDXS])
    ax.bar(ANGLES[IDXS], VALUES, width=WIDTH, color=COLORS, edgecolor="white", linewidth=2)

    add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)

    offset = 0 
    for group, size in zip(unique, GROUPS_SIZE):

        x1 = np.linspace(ANGLES[offset + PAD], ANGLES[offset + size + PAD - 1], num=50)
        ax.plot(x1, [-5] * 50, color="#333333")
        
        ax.text(
            np.mean(x1), -35, "      "+ str(group)+": "+str(size), color="#333333", fontsize=3, 
            fontweight="bold", ha="center", va="center"
        )
        
        x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=50)
        ax.plot(x2, [20] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [40] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [60] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [80] * 50, color="#bebebe", lw=0.8)
        
        offset += size + PAD

    plt.tight_layout() 
    #plt.savefig('plot3.png',dpi=1000)
    plt.show()    

    pass



groupbasedonalgorithm(results,"Total_Time")