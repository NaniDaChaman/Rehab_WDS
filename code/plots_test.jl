using Plots
using DataFrames

# Create a sample DataFrame with binary variables
df = DataFrame(
    A = [1, 0, 1, 0, 1],
    B = [0, 1, 0, 1, 0],
    C = [1, 1, 0, 0, 1]
)

# Convert the DataFrame to a matrix
matrix = Matrix(df)

# Create the heatmap
heatmap(z=matrix)
xlabel!("Variables")
ylabel!("Observations")
title!("Heatmap of Binary Variables")

