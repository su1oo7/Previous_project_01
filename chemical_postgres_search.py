print('patent_explorer_v2')
'''
patent_explorer ver2
using postgresql 
'''
import pandas as pd
import numpy as np
import time
import json 
import warnings
warnings.filterwarnings(action='ignore')

import argparse
import pandas.io.sql as psql
import psycopg2,time
# print(time.asctime(),"rdkit=",rdBase.rdkitVersion)

'''
from rdkit import Chem
from rdkit import DataStructs
from rdkit import rdBase
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
'''

parser = argparse.ArgumentParser(description='Patent Explorer')
parser.add_argument('--input', type=str, default='example.tsv', help='query tsv file')
parser.add_argument('--output', type=str, default='search_result.json', help='output file name')
parser.add_argument('--pickle', default=None)
parser.add_argument('--drop', default=None)
parser.add_argument('--options', type=str)
args = parser.parse_args()

'''
example:

python patent_explorer_psql_proto.py --input '13-input' --output 'search_result.json' --options '{"code":"JP,WO","cutoff":0.85,"num_save":500, "method":'Similarity'}'

'''
PATH = './'
# PATH = '/home/adli/UI/patent_explorer/'

def connect_psql():
    user = # postgresql DB user ID
    password =  # password
    host_product = # 서버 IP
    dbname =  # DB Name
    port= # port number
    
    product_connection_string = "dbname={dbname} user={user} host={host} password={password} port={port}".format(dbname=dbname,
                                    user=user,
                                    host=host_product,
                                    password=password,
                                    port=port) 
    
    try:
        conn = psycopg2.connect(product_connection_string)
        print(f"Success to connect to the database '{dbname}'")
        return conn
    except:
        print("Unable to connect to the database")
        exit()
        
        
def get_option(args):
    
    options = eval(args.options)

    print('[Query]')
    query_df = pd.read_csv(args.input, sep='\t')
    for e, [q, query] in enumerate(zip(query_df['id'].tolist(),  query_df['smiles'].tolist())):
        print(f"({e+1}) {q}: {query}")
    print(f'Total :{e+1}') 
    print()
    options['id'] = query_df['id'].tolist()
    options['smiles'] = query_df['smiles'].tolist()
    
    return options



if __name__ == "__main__":
    conn = connect_psql()
    option = get_option(args)
    
    name = option['id']
    query = option['smiles']
    threshold = option['cutoff']
    num_input = len(name)
    limit = option['num_save']

    output = pd.DataFrame(columns=  ['query_id','query_smiles','num_patent','num_comp','patent','result'])
    start = time.time()
    for i in range(num_input):
        print(option['method'] )
        if option['method'] == 'similarity':
            print('Start Similarity Search')
            df = psql.read_sql(f"select db.code, db.smiles, db.date, res.similarity from db inner join \
                                (select * from get_mfp_neighbors('{query[i]}',{threshold},{int(limit)})) as res \
                                 on db.smiles=res.smiles;",conn)
            df = df.sort_values('similarity',ascending = False).reset_index(drop = True)
#         elif option['method'] == 'substructure':
        else:
            print('Start Substruture Search')
            df = psql.read_sql(f"select db.code, db.smiles, db.date from db where db.smiles in \
                                (select smiles from mol_mfp where m@>cast(mol_to_smarts(mol_from_smiles('{query[i]}'))as text)::qmol \
                                order by smiles limit {int(limit)});",conn)

        df['source'] = 'SureChEMBL' 
        df['country'] = df['code'].apply(lambda x:x[:2])
        df = df[df['country'].isin(option['code'].split(","))].reset_index(drop = True)
        
        num = df.groupby('smiles')['code'].nunique()
        df['num'] = df['smiles'].apply(lambda x: num.loc[x].item())
        kigdom = df.groupby('smiles')['country'].apply(lambda x: list(np.unique(x)))
        df['contries'] = df['smiles'].apply(lambda x: kigdom.loc[x])
    
        if option['method'] == 'similarity':
            df.columns = ['patent_code','smiles','date','similarity','source','each_country_code','num','country_code']
            df = df[['smiles','similarity','country_code','patent_code','date','source','num']]
            df['similarity'] = df['similarity'].apply(lambda x: round(x,3))
      
        else:
            df.columns = ['patent_code','smiles','date','source','each_country_code','num','country_code']
            df = df[['smiles','country_code','patent_code','date','source','num']]
        
        
        #code = df['patent_code'].unique().tolist()
        code = pd.DataFrame(df.groupby('smiles')['patent_code'].apply(lambda x: list(np.unique(x)))).reset_index().to_dict(orient='record')
        df = df.drop_duplicates('smiles',keep='first').reset_index(drop=True)
        df['method'] = option['method']
        
        output.loc[i,'query_id'] = name[i]
        output.loc[i,'query_smiles'] = query[i]
        output.loc[i,'num_patent'] = num.sum()
        output.loc[i,'num_comp'] = df['smiles'].nunique()
        output.loc[i,'patent'] = [code]
        output.loc[i,'result'] = df.to_dict(orient = 'records')
        
        if i == 0: 
            print(f'>> {i+1}st compound done')
        elif i == 1: 
            print(f'>> {i+1}nd compound done')
        elif i == 2: 
            print(f'>> {i+1}rd compound done')
        else:
            print(f'>> {i+1}th compound done')
            
    end = time.time()  
    print(f'Searched patent for {round(end-start)} second')

    output = output.to_dict('list')
    with open(args.output, 'w') as f:
        json.dump(output, f)

    conn.close()

'''
prompt

Success to connect to the database 'patent'
[Query]
(1) Adavosertib: O=C1N(N(C2=NC(NC3=CC=C(C=C3)N4CCN(CC4)C)=NC=C21)C5=NC(C(O)(C)C)=CC=C5)CC=C
(2) Afatinib: CN(C)CC=CC(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)OC4CCOC4
(3) Alectinib: N#CC1=CC2=C(C3=C(N2)C(C)(C4=CC(N5CCC(CC5)N6CCOCC6)=C(C=C4C3=O)CC)C)C=C1
(4) Alisertib: COC1=C(C(=CC=C1)F)C2=NCC3=CN=C(N=C3C4=C2C=C(C=C4)Cl)NC5=CC(=C(C=C5)C(=O)O)OC
(5) Alvocidib: CN1CCC(C(C1)O)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=CC=C4Cl)O)O
(6) AT9283: C1CC1NC(=O)NC2=C(NN=C2)C3=NC4=C(N3)C=C(C=C4)CN5CCOCC5
Total :6

>> 1st compound done
>> 2nd compound done
>> 3rd compound done
>> 4th compound done
>> 5th compound done
>> 6th compound done
Searched patent for 307 second
Done!


'''

    
    
    