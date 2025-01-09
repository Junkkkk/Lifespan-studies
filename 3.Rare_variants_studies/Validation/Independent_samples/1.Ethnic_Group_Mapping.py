import pandas as pd

# Load the dataset containing all samples, including multi-ancestry
df = pd.read_csv('pheno_covar_487149_total.csv')

# Filter out samples based on genetic ethnic grouping: Keep only non-Caucasian samples (Genotype_ethnic == 0)
df = df[df['Genotype_ethnic'] == 0]

# Remove samples with "Do not know" (-1) and "Prefer not to answer" (-3) for self-reported ethnicity
df = df[df['Self_reported_ethnic'] > 0]

# Define a mapping for self-reported ethnic groups to broader categories
ethnic_mapping = {
    1.0: "White", 1001.0: "White", 1002.0: "White", 1003.0: "White",  # White group
    2.0: "Mixed", 2001.0: "Mixed", 2002.0: "Mixed", 2003.0: "Mixed", 2004.0: "Mixed",  # Mixed group
    3.0: "Asian", 5.0: "Asian", 3001.0: "Asian", 3002.0: "Asian", 3003.0: "Asian", 3004.0: "Asian",  # Asian group
    4.0: "Black", 4001.0: "Black", 4002.0: "Black", 4003.0: "Black",  # Black group
    6.0: "Other",  # Other ethnicities
}

# Map self-reported ethnicities to the broader categories defined in ethnic_mapping
df['ethnic_group'] = df['Self_reported_ethnic'].map(ethnic_mapping)

# Create dummy variables for the ethnic groups (excluding the first category as the reference group)
ethnic_dummies = pd.get_dummies(df['ethnic_group'], drop_first=True)

# Concatenate the dummy variables with the original dataframe
df = pd.concat([df, ethnic_dummies], axis=1)

# Save the filtered and processed dataset to a new CSV file
df.to_csv('pheno_covar_73281_independent_samples.csv', index=False)