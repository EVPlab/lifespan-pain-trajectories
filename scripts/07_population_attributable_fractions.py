# ============================================================
# Population Attributable Fraction (PAF) Analysis
# Combines age- and sex-specific RRs with GBD exposure prevalence
# to estimate country- and region-level attributable pain burden.
# ============================================================

import pandas as pd
import numpy as np
import re
import os
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# CONFIGURATION
# ============================================================

# --- Paths: GBD inputs ---
POP_PATH       = "/Users/Patty/Desktop/LifespanPain/PAF/GBD_Estimates/GBD_Pop/Pop_AllRegions.csv"
OB_ADULT_PATH  = "/Users/Patty/Desktop/LifespanPain/PAF/GBD_Estimates/GBD_Obesity/Adult/IHME_GLOBAL_OVERWEIGHT_OBESITY_PREVALENCE_1990_2050_AGES_25_125/IHME_GLOBAL_OVERWEIGHT_OBESITY_PREVALENCE_1990_2050_AGES_25_125_OB_Y2025M02D06.CSV"
OB_ADOL_PATH   = "/Users/Patty/Desktop/LifespanPain/PAF/GBD_Estimates/GBD_Obesity/Adolescent/IHME_GLOBAL_OBESITY_ADOLESCENT_1990_2050_AGES_5_24/IHME_GLOBAL_OBESITY_ADOLESCENT_1990_2050_AGES_5_24_OB_Y2025M02D06.CSV"
SMOKE_PATH     = "/Users/Patty/Desktop/LifespanPain/PAF/GBD_Estimates/GBD_Smoking/GBD_Smoking_2015/IHME_GBD_2015_SMOKING_PREVALENCE_1980_2015_Y2017M04D05.CSV"

# --- Paths: Study-derived inputs ---
RR_PATH        = "/Users/Patty/Desktop/LifespanPain/PAF/Spline_RR/SplineRR_UNADJ_AgeSpecific_SexStratified.csv"
BASELINE_DIR   = "/Users/Patty/Desktop/LifespanPain/Results/Global_trajectories/Stats/"

# --- Paths: Intermediate cleaned files ---
CLEAN_DIR      = "/Users/Patty/Desktop/LifespanPain/PAF/GBD_Estimates/"
POP_CLEAN      = os.path.join(CLEAN_DIR, "Population_Cleaned.csv")
OB_CLEAN       = os.path.join(CLEAN_DIR, "Obesity_Master_Cleaned.csv")
SMOKE_CLEAN    = os.path.join(CLEAN_DIR, "Smoking_Master_Cleaned.csv")
RR_PIVOT_PATH  = "/Users/Patty/Desktop/LifespanPain/PAF/Spline_RR/RR_Pivot_Master.csv"
BASELINE_CLEAN = os.path.join(BASELINE_DIR, "Baseline_Pain_Cleaned.csv")

# --- Paths: Outputs ---
RESULTS_DIR    = "/Users/Patty/Desktop/LifespanPain/Results/"
FIGURES_DIR    = "/Users/Patty/Desktop/TMP_Traj/"
OUTPUT_MAIN    = os.path.join(RESULTS_DIR, "Global_Impact_Map_Data_SexStrat.csv")
OUTPUT_RATES   = os.path.join(RESULTS_DIR, "Global_Pain_Burden_Rates_100k.csv")
TABLE6_PATH    = os.path.join(RESULTS_DIR, "Supplementary_Table/Appendix_Table6_Regional_PAF.csv")

# --- Country name crosswalk (GBD exposure names → GBD population names) ---
NAME_CROSSWALK = {
    "Bolivia": "Bolivia (Plurinational State of)",
    "Brunei": "Brunei Darussalam",
    "Cape Verde": "Cabo Verde",
    "Cote d'Ivoire": "Côte d'Ivoire",
    "Czech Republic": "Czechia",
    "Federated States of Micronesia": "Micronesia (Federated States of)",
    "Iran": "Iran (Islamic Republic of)",
    "Laos": "Lao People's Democratic Republic",
    "Macedonia": "North Macedonia",
    "Moldova": "Republic of Moldova",
    "North Korea": "Democratic People's Republic of Korea",
    "Russia": "Russian Federation",
    "South Korea": "Republic of Korea",
    "Swaziland": "Eswatini",
    "Syria": "Syrian Arab Republic",
    "Tanzania": "United Republic of Tanzania",
    "The Bahamas": "Bahamas",
    "The Gambia": "Gambia",
    "Turkey": "Türkiye",
    "United States": "United States of America",
    "Venezuela": "Venezuela (Bolivarian Republic of)",
    "Vietnam": "Viet Nam",
    "Virgin Islands, U.S.": "United States Virgin Islands",
    "Taiwan (Province of China)": "Taiwan (Province of China)",
}

# --- Natural Earth → GBD name fixes (for choropleths) ---
NE_NAME_MAP = {
    'Tanzania': 'United Republic of Tanzania',
    'Dem. Rep. Congo': 'Democratic Republic of the Congo',
    'Russia': 'Russian Federation',
    'Bolivia': 'Bolivia (Plurinational State of)',
    'Venezuela': 'Venezuela (Bolivarian Republic of)',
    "Cote d'Ivoire": "Côte d'Ivoire",
    'Laos': "Lao People's Democratic Republic",
    'Vietnam': 'Viet Nam',
    'North Korea': "Democratic People's Republic of Korea",
    'South Korea': 'Republic of Korea',
    'Iran': 'Iran (Islamic Republic of)',
    'Syria': 'Syrian Arab Republic',
    'Moldova': 'Republic of Moldova',
    'Brunei': 'Brunei Darussalam',
    'Bosnia and Herz.': 'Bosnia and Herzegovina',
    'Dominican Rep.': 'Dominican Republic',
    'Central African Rep.': 'Central African Republic',
    'S. Sudan': 'South Sudan',
    'Solomon Is.': 'Solomon Islands',
    'Eq. Guinea': 'Equatorial Guinea',
    'eSwatini': 'Eswatini',
}

# --- GBD aggregate locations to exclude from country-level analysis ---
GBD_AGGREGATES = {
    'Andean Latin America', 'Australasia', 'Caribbean', 'Central Asia',
    'Central Europe', 'Central Europe, Eastern Europe, and Central Asia',
    'Central Latin America', 'Central Sub-Saharan Africa', 'East Asia',
    'Eastern Europe', 'Eastern Sub-Saharan Africa', 'Global', 'High-income',
    'High-income Asia Pacific', 'High-income North America',
    'Latin America and Caribbean', 'North Africa and Middle East', 'Oceania',
    'South Asia', 'Southeast Asia', 'Southeast Asia, East Asia, and Oceania',
    'Southern Latin America', 'Southern Sub-Saharan Africa',
    'Sub-Saharan Africa', 'Tropical Latin America', 'Western Europe',
    'Western Sub-Saharan Africa',
}

# --- SDG Region Mapping ---
SDG_REGION_MAP = {
    # Sub-Saharan Africa
    'Angola': 'Sub-Saharan Africa', 'Benin': 'Sub-Saharan Africa',
    'Botswana': 'Sub-Saharan Africa', 'Burkina Faso': 'Sub-Saharan Africa',
    'Burundi': 'Sub-Saharan Africa', 'Cabo Verde': 'Sub-Saharan Africa',
    'Cameroon': 'Sub-Saharan Africa', 'Central African Republic': 'Sub-Saharan Africa',
    'Chad': 'Sub-Saharan Africa', 'Comoros': 'Sub-Saharan Africa',
    'Congo': 'Sub-Saharan Africa', 'Democratic Republic of the Congo': 'Sub-Saharan Africa',
    "Côte d'Ivoire": 'Sub-Saharan Africa', 'Djibouti': 'Sub-Saharan Africa',
    'Equatorial Guinea': 'Sub-Saharan Africa', 'Eritrea': 'Sub-Saharan Africa',
    'Eswatini': 'Sub-Saharan Africa', 'Ethiopia': 'Sub-Saharan Africa',
    'Gabon': 'Sub-Saharan Africa', 'Gambia': 'Sub-Saharan Africa',
    'Ghana': 'Sub-Saharan Africa', 'Guinea': 'Sub-Saharan Africa',
    'Guinea-Bissau': 'Sub-Saharan Africa', 'Kenya': 'Sub-Saharan Africa',
    'Lesotho': 'Sub-Saharan Africa', 'Liberia': 'Sub-Saharan Africa',
    'Madagascar': 'Sub-Saharan Africa', 'Malawi': 'Sub-Saharan Africa',
    'Mali': 'Sub-Saharan Africa', 'Mauritania': 'Sub-Saharan Africa',
    'Mauritius': 'Sub-Saharan Africa', 'Mozambique': 'Sub-Saharan Africa',
    'Namibia': 'Sub-Saharan Africa', 'Niger': 'Sub-Saharan Africa',
    'Nigeria': 'Sub-Saharan Africa', 'Rwanda': 'Sub-Saharan Africa',
    'Sao Tome and Principe': 'Sub-Saharan Africa', 'Senegal': 'Sub-Saharan Africa',
    'Seychelles': 'Sub-Saharan Africa', 'Sierra Leone': 'Sub-Saharan Africa',
    'Somalia': 'Sub-Saharan Africa', 'South Africa': 'Sub-Saharan Africa',
    'South Sudan': 'Sub-Saharan Africa', 'Sudan': 'Sub-Saharan Africa',
    'Togo': 'Sub-Saharan Africa', 'Uganda': 'Sub-Saharan Africa',
    'United Republic of Tanzania': 'Sub-Saharan Africa',
    'Zambia': 'Sub-Saharan Africa', 'Zimbabwe': 'Sub-Saharan Africa',
    # Northern Africa and Western Asia
    'Algeria': 'Northern Africa and Western Asia', 'Bahrain': 'Northern Africa and Western Asia',
    'Egypt': 'Northern Africa and Western Asia', 'Iraq': 'Northern Africa and Western Asia',
    'Iran (Islamic Republic of)': 'Northern Africa and Western Asia',
    'Israel': 'Northern Africa and Western Asia', 'Jordan': 'Northern Africa and Western Asia',
    'Kuwait': 'Northern Africa and Western Asia', 'Lebanon': 'Northern Africa and Western Asia',
    'Libya': 'Northern Africa and Western Asia', 'Morocco': 'Northern Africa and Western Asia',
    'Oman': 'Northern Africa and Western Asia', 'Palestine': 'Northern Africa and Western Asia',
    'Qatar': 'Northern Africa and Western Asia', 'Saudi Arabia': 'Northern Africa and Western Asia',
    'Syrian Arab Republic': 'Northern Africa and Western Asia',
    'Tunisia': 'Northern Africa and Western Asia', 'Turkey': 'Northern Africa and Western Asia',
    'United Arab Emirates': 'Northern Africa and Western Asia',
    'Yemen': 'Northern Africa and Western Asia', 'Cyprus': 'Northern Africa and Western Asia',
    # Central and Southern Asia
    'Afghanistan': 'Central and Southern Asia', 'Bangladesh': 'Central and Southern Asia',
    'Bhutan': 'Central and Southern Asia', 'India': 'Central and Southern Asia',
    'Kazakhstan': 'Central and Southern Asia', 'Kyrgyzstan': 'Central and Southern Asia',
    'Maldives': 'Central and Southern Asia', 'Nepal': 'Central and Southern Asia',
    'Pakistan': 'Central and Southern Asia', 'Sri Lanka': 'Central and Southern Asia',
    'Tajikistan': 'Central and Southern Asia', 'Turkmenistan': 'Central and Southern Asia',
    'Uzbekistan': 'Central and Southern Asia',
    # Eastern and South-Eastern Asia
    'Brunei Darussalam': 'Eastern and South-Eastern Asia',
    'Cambodia': 'Eastern and South-Eastern Asia', 'China': 'Eastern and South-Eastern Asia',
    "Democratic People's Republic of Korea": 'Eastern and South-Eastern Asia',
    'Indonesia': 'Eastern and South-Eastern Asia', 'Japan': 'Eastern and South-Eastern Asia',
    "Lao People's Democratic Republic": 'Eastern and South-Eastern Asia',
    'Malaysia': 'Eastern and South-Eastern Asia', 'Mongolia': 'Eastern and South-Eastern Asia',
    'Myanmar': 'Eastern and South-Eastern Asia', 'Philippines': 'Eastern and South-Eastern Asia',
    'Republic of Korea': 'Eastern and South-Eastern Asia',
    'Singapore': 'Eastern and South-Eastern Asia', 'Thailand': 'Eastern and South-Eastern Asia',
    'Timor-Leste': 'Eastern and South-Eastern Asia', 'Viet Nam': 'Eastern and South-Eastern Asia',
    # Latin America and the Caribbean
    'Antigua and Barbuda': 'Latin America and the Caribbean',
    'Argentina': 'Latin America and the Caribbean', 'Bahamas': 'Latin America and the Caribbean',
    'Barbados': 'Latin America and the Caribbean', 'Belize': 'Latin America and the Caribbean',
    'Bolivia (Plurinational State of)': 'Latin America and the Caribbean',
    'Brazil': 'Latin America and the Caribbean', 'Chile': 'Latin America and the Caribbean',
    'Colombia': 'Latin America and the Caribbean', 'Costa Rica': 'Latin America and the Caribbean',
    'Cuba': 'Latin America and the Caribbean', 'Dominica': 'Latin America and the Caribbean',
    'Dominican Republic': 'Latin America and the Caribbean',
    'Ecuador': 'Latin America and the Caribbean', 'El Salvador': 'Latin America and the Caribbean',
    'Grenada': 'Latin America and the Caribbean', 'Guatemala': 'Latin America and the Caribbean',
    'Guyana': 'Latin America and the Caribbean', 'Haiti': 'Latin America and the Caribbean',
    'Honduras': 'Latin America and the Caribbean', 'Jamaica': 'Latin America and the Caribbean',
    'Mexico': 'Latin America and the Caribbean', 'Nicaragua': 'Latin America and the Caribbean',
    'Panama': 'Latin America and the Caribbean', 'Paraguay': 'Latin America and the Caribbean',
    'Peru': 'Latin America and the Caribbean', 'Puerto Rico': 'Latin America and the Caribbean',
    'Saint Lucia': 'Latin America and the Caribbean',
    'Saint Vincent and the Grenadines': 'Latin America and the Caribbean',
    'Suriname': 'Latin America and the Caribbean',
    'Trinidad and Tobago': 'Latin America and the Caribbean',
    'Uruguay': 'Latin America and the Caribbean',
    'Venezuela (Bolivarian Republic of)': 'Latin America and the Caribbean',
    'United States Virgin Islands': 'Latin America and the Caribbean',
    # North America
    'Canada': 'North America', 'Greenland': 'North America',
    'United States of America': 'North America',
    # Western Europe
    'Andorra': 'Western Europe', 'Austria': 'Western Europe',
    'Belgium': 'Western Europe', 'Denmark': 'Western Europe',
    'Finland': 'Western Europe', 'France': 'Western Europe',
    'Germany': 'Western Europe', 'Greece': 'Western Europe',
    'Iceland': 'Western Europe', 'Ireland': 'Western Europe',
    'Italy': 'Western Europe', 'Luxembourg': 'Western Europe',
    'Malta': 'Western Europe', 'Netherlands': 'Western Europe',
    'Norway': 'Western Europe', 'Portugal': 'Western Europe',
    'Spain': 'Western Europe', 'Sweden': 'Western Europe',
    'Switzerland': 'Western Europe', 'United Kingdom': 'Western Europe',
    # Eastern Europe
    'Albania': 'Eastern Europe', 'Belarus': 'Eastern Europe',
    'Bosnia and Herzegovina': 'Eastern Europe', 'Bulgaria': 'Eastern Europe',
    'Croatia': 'Eastern Europe', 'Czechia': 'Eastern Europe',
    'Estonia': 'Eastern Europe', 'Hungary': 'Eastern Europe',
    'Latvia': 'Eastern Europe', 'Lithuania': 'Eastern Europe',
    'Montenegro': 'Eastern Europe', 'North Macedonia': 'Eastern Europe',
    'Poland': 'Eastern Europe', 'Republic of Moldova': 'Eastern Europe',
    'Romania': 'Eastern Europe', 'Russian Federation': 'Eastern Europe',
    'Serbia': 'Eastern Europe', 'Slovakia': 'Eastern Europe',
    'Slovenia': 'Eastern Europe', 'Ukraine': 'Eastern Europe',
    'Armenia': 'Eastern Europe', 'Azerbaijan': 'Eastern Europe',
    'Georgia': 'Eastern Europe',
    # Oceania
    'American Samoa': 'Oceania', 'Australia': 'Oceania', 'Bermuda': 'Oceania',
    'Fiji': 'Oceania', 'Guam': 'Oceania', 'Kiribati': 'Oceania',
    'Marshall Islands': 'Oceania', 'Micronesia (Federated States of)': 'Oceania',
    'New Zealand': 'Oceania', 'Northern Mariana Islands': 'Oceania',
    'Papua New Guinea': 'Oceania', 'Samoa': 'Oceania', 'Solomon Islands': 'Oceania',
    'Tonga': 'Oceania', 'Vanuatu': 'Oceania',
}

REGION_SHORT = {
    'Sub-Saharan Africa': 'SS Africa',
    'Northern Africa and Western Asia': 'N Africa & W Asia',
    'Central and Southern Asia': 'C & S Asia',
    'Eastern and South-Eastern Asia': 'E & SE Asia',
    'Latin America and the Caribbean': 'Latin America',
    'North America': 'N America',
    'Western Europe': 'W Europe',
    'Eastern Europe': 'E Europe',
    'Oceania': 'Oceania',
}

# --- Baseline prevalence file list ---
BASELINE_FILES = [
    "ElbowPain_Global_PastMonth_Results.csv", "FacialPain_Global_PastMonth_Results.csv",
    "BackPain_Global_PastMonth_Results.csv", "StomachAbdominalPain_Global_PastMonth_Results.csv",
    "Pain_Global_PastMonth_Results.csv", "FootPain_Global_PastMonth_Results.csv",
    "KneePain_Global_PastMonth_Results.csv", "HipPain_Global_PastMonth_Results.csv",
    "HandPain_Global_PastMonth_Results.csv", "Headache_Global_PastMonth_Results.csv",
    "JointPain_Global_PastMonth_Results.csv", "NeckShoulderPain_Global_PastMonth_Results.csv",
    "ChestPain_Global_PastMonth_Results.csv",
]

# --- Pain site display names ---
SITE_NAME_MAP = {
    'BackPain': 'Back', 'KneePain': 'Knee', 'Headache': 'Headache',
    'HipPain': 'Hip', 'NeckShoulderPain': 'Neck/Shoulder',
    'HandPain': 'Hand/Wrist', 'FootPain': 'Foot/Ankle',
    'StomachAbdominalPain': 'Stomach/Abdominal', 'ElbowPain': 'Elbow',
    'ChestPain': 'Chest', 'FacialPain': 'Facial',
    'Pain': 'Any Bodily Pain', 'JointPain': 'Joint Pain',
}

# ============================================================
# SECTION 1: CLEAN GBD POPULATION
# ============================================================
print("=" * 60)
print("SECTION 1: Cleaning GBD Population data")
print("=" * 60)

df_pop = pd.read_csv(POP_PATH)
df_pop = df_pop[df_pop['sex'].isin(['Male', 'Female'])]

keep_ages = {
    '5-9 years': 5, '10-14 years': 10, '15-19 years': 15, '20-24 years': 20,
    '25-29 years': 25, '30-34 years': 30, '35-39 years': 35, '40-44 years': 40,
    '45-49 years': 45, '50-54 years': 50, '55-59 years': 55, '60-64 years': 60,
    '65-69 years': 65, '70-74 years': 70, '75-79 years': 75, '80+ years': 80,
}
df_pop = df_pop[df_pop['age'].isin(keep_ages.keys())]
df_pop['age_start'] = df_pop['age'].map(keep_ages)

pop_final = df_pop[['location', 'age_start', 'sex', 'val']].copy()
pop_final.columns = ['location_name', 'age_start', 'sex', 'pop_count']
pop_final.to_csv(POP_CLEAN, index=False)
print(f"  Saved: {POP_CLEAN}")

# ============================================================
# SECTION 2: CLEAN GBD OBESITY PREVALENCE
# ============================================================
print("\n" + "=" * 60)
print("SECTION 2: Cleaning GBD Obesity prevalence")
print("=" * 60)

def clean_obesity_file(path):
    df = pd.read_csv(path)
    df = df[(df['year_id'] == 2021) & (df['sex_id'].isin([1, 2])) & (df['metric'] == 'obesity')]
    df['sex'] = df['sex_id'].map({1: 'Male', 2: 'Female'})

    def map_age(age_str):
        age_str = str(age_str).lower()
        if any(x in age_str for x in ['80', '85', '90', '95']):
            return 80
        nums = re.findall(r'\d+', age_str)
        return int(nums[0]) if nums else None

    df['age_start'] = df['age_group_name'].apply(map_age)
    return df.groupby(['location_name', 'age_start', 'sex'])['mean_prev'].mean().reset_index()

df_ob = pd.concat([clean_obesity_file(OB_ADOL_PATH), clean_obesity_file(OB_ADULT_PATH)], axis=0)
df_ob['location_name'] = df_ob['location_name'].replace(NAME_CROSSWALK)
df_ob.to_csv(OB_CLEAN, index=False)
print(f"  Saved: {OB_CLEAN}")

# ============================================================
# SECTION 3: CLEAN GBD SMOKING PREVALENCE
# ============================================================
print("\n" + "=" * 60)
print("SECTION 3: Cleaning GBD Smoking prevalence")
print("=" * 60)

df_smoke = pd.read_csv(SMOKE_PATH)
df_smoke = df_smoke[
    (df_smoke['year_id'] == 2015) &
    (df_smoke['metric'] == 'Percent') &
    (df_smoke['sex_id'].isin([1, 2]))
]
df_smoke['sex'] = df_smoke['sex_id'].map({1: 'Male', 2: 'Female'})

def map_smoke_age(age_str):
    age_str = str(age_str).lower()
    if '80' in age_str:
        return 80
    nums = re.findall(r'\d+', age_str)
    return int(nums[0]) if nums else None

df_smoke['age_start'] = df_smoke['age_group_name'].apply(map_smoke_age)
df_smoke_final = df_smoke[['location_name', 'age_start', 'sex', 'mean']].copy()
df_smoke_final.rename(columns={'mean': 'smoke_prev'}, inplace=True)

# Backfill age 5 with 0 prevalence (no smoking data for young children)
locations = df_smoke_final['location_name'].unique()
backfill = pd.concat([
    pd.DataFrame({'location_name': locations, 'age_start': 5, 'sex': 'Male', 'smoke_prev': 0.0}),
    pd.DataFrame({'location_name': locations, 'age_start': 5, 'sex': 'Female', 'smoke_prev': 0.0}),
])
df_smoke_final = pd.concat([df_smoke_final, backfill], axis=0).sort_values(
    ['location_name', 'sex', 'age_start'])

df_smoke_final['location_name'] = df_smoke_final['location_name'].replace(NAME_CROSSWALK)
df_smoke_final.to_csv(SMOKE_CLEAN, index=False)
print(f"  Saved: {SMOKE_CLEAN}")

# ============================================================
# SECTION 4: CLEAN AND PIVOT RISK RATIOS
# ============================================================
print("\n" + "=" * 60)
print("SECTION 4: Cleaning and pivoting Risk Ratios")
print("=" * 60)

df_rr = pd.read_csv(RR_PATH)
df_rr['sex'] = df_rr['SexStratum'].astype(str).replace({'1': 'Male', '2': 'Female'})

contrast_map = {
    'mean(Ex-Smoker) / mean(Non-Smoker)': 'Ex-Smoker',
    'mean(Smoker) / mean(Non-Smoker)': 'Smoker',
    'mean(Obese) / mean(Normal)': 'Obese',
    'mean(Overweight) / mean(Normal)': 'Overweight',
    'mean(Lower Middle) / mean(Upper)': 'Lower Middle',
    'mean(Lower) / mean(Upper)': 'Lower',
    'mean(Middle) / mean(Upper)': 'Middle',
    'mean(Upper Middle) / mean(Upper)': 'Upper Middle',
}
df_rr['risk_category'] = df_rr['contrast'].map(contrast_map)

# Backfill ages 5 and 10 from age 15 (no exposure RR data for youngest children)
backfill_rows = []
for site in df_rr['Pain_Site'].unique():
    for risk in df_rr['risk_category'].unique():
        for s in ['Male', 'Female']:
            val_15 = df_rr[
                (df_rr['Pain_Site'] == site) &
                (df_rr['risk_category'] == risk) &
                (df_rr['sex'] == s) &
                (df_rr['Age'] == 15)
            ]
            if not val_15.empty:
                for new_age in [5, 10]:
                    new_row = val_15.iloc[0].copy()
                    new_row['Age'] = new_age
                    backfill_rows.append(new_row)

df_rr = pd.concat([df_rr, pd.DataFrame(backfill_rows)], ignore_index=True)

# Cap at 80+ to match GBD age bins
df_rr['Age'] = df_rr['Age'].apply(lambda x: 80 if x >= 80 else x)

# Pivot to wide format
df_pivot = df_rr[['Pain_Site', 'Age', 'sex', 'risk_category', 'RR']].pivot_table(
    index=['Pain_Site', 'Age', 'sex'],
    columns='risk_category',
    values='RR',
).reset_index()
df_pivot.columns = [c.replace(' ', '_') for c in df_pivot.columns]

df_pivot.to_csv(RR_PIVOT_PATH, index=False)
print(f"  Saved: {RR_PIVOT_PATH}")

# ============================================================
# SECTION 5: CLEAN BASELINE PAIN PREVALENCE
# ============================================================
print("\n" + "=" * 60)
print("SECTION 5: Cleaning baseline pain prevalence from trajectory models")
print("=" * 60)

def map_to_gbd_age(age):
    if age >= 80:
        return 80
    return int((age // 5) * 5)

all_sites_data = []
for filename in BASELINE_FILES:
    full_path = os.path.join(BASELINE_DIR, filename)
    if not os.path.exists(full_path):
        continue

    site_name = filename.split('_')[0]
    df = pd.read_csv(full_path)
    df['Pain_Site'] = site_name

    # Keep only Male/Female (exclude "Both" to avoid double-counting)
    df = df[df['Sex'].isin(['Male', 'Female'])]

    df['age_start'] = df['Age'].apply(map_to_gbd_age)
    df_binned = df.groupby(['Pain_Site', 'age_start', 'Sex'])['Predicted_Prob'].mean().reset_index()
    df_binned.rename(columns={'Sex': 'sex', 'Predicted_Prob': 'Base_Pain_Prev'}, inplace=True)
    all_sites_data.append(df_binned)

master_baseline = pd.concat(all_sites_data, axis=0)
master_baseline.to_csv(BASELINE_CLEAN, index=False)
print(f"  Saved: {BASELINE_CLEAN}")

# ============================================================
# SECTION 6: COMPUTE POPULATION ATTRIBUTABLE FRACTIONS
# ============================================================
print("\n" + "=" * 60)
print("SECTION 6: Computing PAFs")
print("=" * 60)

# Load cleaned inputs
df_pop   = pd.read_csv(POP_CLEAN)
df_ob    = pd.read_csv(OB_CLEAN).rename(columns={'mean_prev': 'obesity_prev'})
df_smoke = pd.read_csv(SMOKE_CLEAN)
df_rr    = pd.read_csv(RR_PIVOT_PATH).rename(columns={'Age': 'age_start'})
df_base  = pd.read_csv(BASELINE_CLEAN)

# Merge country-level GBD exposure and population
df_countries = pd.merge(df_pop, df_ob, on=['location_name', 'age_start', 'sex'], how='inner')
df_countries = pd.merge(df_countries, df_smoke, on=['location_name', 'age_start', 'sex'], how='inner')

# Merge global study parameters (RRs and baseline prevalence)
df_params = pd.merge(df_rr, df_base, on=['Pain_Site', 'age_start', 'sex'], how='inner')

# Grand merge: broadcast 11-site parameters to every country
df_final = pd.merge(df_countries, df_params, on=['age_start', 'sex'], how='inner')

# --- A. INCOME PAF: Partial polytomous (Lower quintile vs population mean) ---
p_q = 0.20
df_final['Den_Income'] = (
    p_q * 1.0 +
    p_q * df_final['Upper_Middle'] +
    p_q * df_final['Middle'] +
    p_q * df_final['Lower_Middle'] +
    p_q * df_final['Lower']
)
df_final['PAF_Income'] = (p_q * (df_final['Lower'] - 1.0)) / df_final['Den_Income']

# --- B. OBESITY & SMOKING PAF: Standard binary formula ---
df_final['PAF_Obese'] = (
    (df_final['obesity_prev'] * (df_final['Obese'] - 1)) /
    (1 + df_final['obesity_prev'] * (df_final['Obese'] - 1))
)
df_final['PAF_Smoke'] = (
    (df_final['smoke_prev'] * (df_final['Smoker'] - 1)) /
    (1 + df_final['smoke_prev'] * (df_final['Smoker'] - 1))
)

# --- C. Floor at 0 (clean numeric noise) ---
for col in ['PAF_Income', 'PAF_Obese', 'PAF_Smoke']:
    df_final[col] = df_final[col].clip(lower=0)

# --- D. COMBINED PAF: Multiplicative independence ---
df_final['PAF_Joint'] = 1 - (
    (1 - df_final['PAF_Obese']) *
    (1 - df_final['PAF_Smoke']) *
    (1 - df_final['PAF_Income'])
)

# --- E. Proportional segment decomposition (non-overlapping, sum to joint) ---
paf_sum = df_final['PAF_Obese'] + df_final['PAF_Smoke'] + df_final['PAF_Income']
paf_sum_safe = paf_sum.replace(0, np.nan)

df_final['PAF_Obese_Seg']  = (df_final['PAF_Obese']  / paf_sum_safe) * df_final['PAF_Joint']
df_final['PAF_Smoke_Seg']  = (df_final['PAF_Smoke']  / paf_sum_safe) * df_final['PAF_Joint']
df_final['PAF_Income_Seg'] = (df_final['PAF_Income'] / paf_sum_safe) * df_final['PAF_Joint']

df_final[['PAF_Obese_Seg', 'PAF_Smoke_Seg', 'PAF_Income_Seg']] = \
    df_final[['PAF_Obese_Seg', 'PAF_Smoke_Seg', 'PAF_Income_Seg']].fillna(0)

# --- F. Absolute case counts ---
df_final['Total_Pain_Cases'] = df_final['pop_count'] * df_final['Base_Pain_Prev']

# Individual (may overlap)
df_final['Cases_Joint']   = df_final['Total_Pain_Cases'] * df_final['PAF_Joint']
df_final['Cases_Obesity'] = df_final['Total_Pain_Cases'] * df_final['PAF_Obese']
df_final['Cases_Smoking'] = df_final['Total_Pain_Cases'] * df_final['PAF_Smoke']
df_final['Cases_Income']  = df_final['Total_Pain_Cases'] * df_final['PAF_Income']

# Segments (non-overlapping, sum to joint)
df_final['Cases_Obesity_Seg'] = df_final['Total_Pain_Cases'] * df_final['PAF_Obese_Seg']
df_final['Cases_Smoke_Seg']   = df_final['Total_Pain_Cases'] * df_final['PAF_Smoke_Seg']
df_final['Cases_Income_Seg']  = df_final['Total_Pain_Cases'] * df_final['PAF_Income_Seg']

# ============================================================
# SECTION 7: AGGREGATE TO COUNTRY × SITE LEVEL
# ============================================================
print("\n" + "=" * 60)
print("SECTION 7: Aggregating to country × site level")
print("=" * 60)

map_summary = df_final.groupby(['location_name', 'Pain_Site']).agg({
    'pop_count': 'sum',
    'Total_Pain_Cases': 'sum',
    'Cases_Joint': 'sum',
    'Cases_Obesity': 'sum',
    'Cases_Smoking': 'sum',
    'Cases_Income': 'sum',
    'Cases_Obesity_Seg': 'sum',
    'Cases_Smoke_Seg': 'sum',
    'Cases_Income_Seg': 'sum',
}).reset_index()

# PAF percentages
map_summary['PAF_Joint_Pct']   = 100 * map_summary['Cases_Joint']   / map_summary['Total_Pain_Cases']
map_summary['PAF_Obesity_Pct'] = 100 * map_summary['Cases_Obesity'] / map_summary['Total_Pain_Cases']
map_summary['PAF_Smoking_Pct'] = 100 * map_summary['Cases_Smoking'] / map_summary['Total_Pain_Cases']
map_summary['PAF_Income_Pct']  = 100 * map_summary['Cases_Income']  / map_summary['Total_Pain_Cases']

# Segment PAF percentages (for stacked bars)
map_summary['Seg_Obesity_Pct'] = 100 * map_summary['Cases_Obesity_Seg'] / map_summary['Total_Pain_Cases']
map_summary['Seg_Smoke_Pct']   = 100 * map_summary['Cases_Smoke_Seg']   / map_summary['Total_Pain_Cases']
map_summary['Seg_Income_Pct']  = 100 * map_summary['Cases_Income_Seg']  / map_summary['Total_Pain_Cases']

# Rates per 100k
map_summary['Joint_Cases_per_100k']   = 100000 * map_summary['Cases_Joint']   / map_summary['pop_count']
map_summary['Obesity_Cases_per_100k'] = 100000 * map_summary['Cases_Obesity'] / map_summary['pop_count']
map_summary['Smoking_Cases_per_100k'] = 100000 * map_summary['Cases_Smoking'] / map_summary['pop_count']
map_summary['Income_Cases_per_100k']  = 100000 * map_summary['Cases_Income']  / map_summary['pop_count']

# Save
map_summary.to_csv(OUTPUT_MAIN, index=False)
map_summary.to_csv(OUTPUT_RATES, index=False)

print(f"  Saved: {OUTPUT_MAIN}")
print(f"  Countries: {map_summary['location_name'].nunique()}")
print(f"\n  Top 10 by Joint PAF:")
print(map_summary.sort_values('PAF_Joint_Pct', ascending=False)[
    ['location_name', 'Pain_Site', 'PAF_Joint_Pct']].head(10).to_string(index=False))

# ============================================================
# SECTION 8: GLOBAL SITE-SPECIFIC PAFs
# ============================================================
print("\n" + "=" * 60)
print("SECTION 8: Global site-specific PAFs")
print("=" * 60)

global_paf = map_summary.groupby('Pain_Site').agg(
    Total_Pain_Cases=('Total_Pain_Cases', 'sum'),
    Cases_Obesity=('Cases_Obesity', 'sum'),
    Cases_Smoking=('Cases_Smoking', 'sum'),
    Cases_Income=('Cases_Income', 'sum'),
    Cases_Joint=('Cases_Joint', 'sum'),
).reset_index()

global_paf['PAF_Obesity_Pct'] = 100 * global_paf['Cases_Obesity'] / global_paf['Total_Pain_Cases']
global_paf['PAF_Smoking_Pct'] = 100 * global_paf['Cases_Smoking'] / global_paf['Total_Pain_Cases']
global_paf['PAF_Income_Pct']  = 100 * global_paf['Cases_Income']  / global_paf['Total_Pain_Cases']
global_paf['PAF_Joint_Pct']   = 100 * global_paf['Cases_Joint']   / global_paf['Total_Pain_Cases']

print(global_paf[['Pain_Site', 'PAF_Obesity_Pct', 'PAF_Smoking_Pct',
                   'PAF_Income_Pct', 'PAF_Joint_Pct']].round(1).to_string(index=False))
