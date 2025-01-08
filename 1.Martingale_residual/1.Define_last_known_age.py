import pandas as pd

#For survival subjects,

y=pd.read_csv('ukb_Year_of_birth.csv')
m=pd.read_csv('ukb_Month_of_birth.csv')

ym=pd.merge(y,m,on='eid')
ym['Birth'] = pd.to_datetime(ym[['34-0.0', '52-0.0']].astype(str).apply('-'.join, axis=1))

# the latest censoring date (November 30, 2022)
ym['Current'] = pd.to_datetime('2022-11-30')

ym['age_years'] = ym['Current'].dt.year - ym['Birth'].dt.year
ym['age_months'] = ym['Current'].dt.month - ym['Birth'].dt.month

ym['age'] = ym['age_years']
ym.loc[ym['age_months'] < 0, 'age'] -= 1
ym.loc[ym['age_months'] < 0, 'age_months'] += 12

ym['last_known_age'] = ym['age'] + ym['age_months']/12