import pandas as pd

# Load your file
# df = pd.read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/3. Pos_annoated_exogenous.csv")
df = pd.read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/4. Neg_annoated_exogenous.csv") 

# List of IDs to remove (both with and without _r1/_r2)
remove_ids = [
    "BI00707_r1", "BI00707_r2", "BI00709_r1", "BI00709_r2",
    "BI00712_r1", "BI00712_r2", "BI00715_r1", "BI00715_r2",
    "BI00724_r1", "BI00724_r2", "BI00727_r1", "BI00727_r2",
    "BI00730_r1", "BI00730_r2", "BI00741_r1", "BI00741_r2",
    "BI00759_r1", "BI00759_r2", "BI00746_r1", "BI00746_r2",
    "BI00753_r1", "BI00753_r2",
    "BI00707", "BI00709", "BI00712", "BI00715",
    "BI00724", "BI00727", "BI00730", "BI00741",
    "BI00759", "BI00746", "BI00753"
]

# Remove columns that contain any of those IDs
df = df.loc[:, ~df.columns.isin(remove_ids)]

# Remove rows where any of those IDs appear in any cell
mask = df.apply(lambda x: x.astype(str).str.contains('|'.join(remove_ids)).any(), axis=1)
df = df.loc[~mask]

# Save cleaned file
# output_path = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Pos_annoated_exogenous_CLEAN.csv"
output_path = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Neg_annoated_exogenous_CLEAN.csv"

df.to_csv(output_path, index=False)

print(f"Final shape: {df.shape}")
