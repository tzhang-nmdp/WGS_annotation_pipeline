
#// ...existing code...
import sys
import os
import pandas as pd
import chardet

if len(sys.argv) < 4:
    raise SystemExit("Usage: python annovar_db_table.py <input_file> <db> <db_type>")

input_path = str(sys.argv[1])
input_dir = os.path.dirname(input_path)
input_base = os.path.basename(input_path)
output_file = input_path + '.db.csv'    # or use input_base + '.db.csv' if you prefer

db = str(sys.argv[2])
db_type = str(sys.argv[3])  # 'generic' or 'region' or 'gene'

# Build db file path(s) using the base name
if db_type == 'generic':
    db_file = os.path.join(input_dir,'annovarg', f'{input_base}.{db}.hg38_generic_dropped')
elif db_type == 'region':
    db_file = os.path.join(input_dir, 'annovarr', f'{input_base}.{db}.hg38_{db}_sv')
elif db_type == 'gene':
    db_file_all = os.path.join(input_dir, 'annovargg', f'{input_base}.{db}.variant_function')
    db_file_exon = os.path.join(input_dir, 'annovargg', f'{input_base}.{db}.exonic_variant_function')
    # prefer exon file if present
    if os.path.exists(db_file_exon):
        db_file = db_file_exon
    elif os.path.exists(db_file_all):
        db_file = db_file_all
    else:
        raise FileNotFoundError(f"Neither {db_file_exon} nor {db_file_all} found")
elif db_type == 'cadd':
    db_file = os.path.join(input_dir, 'cadd', f'{input_base}.{db}.score_all')
else:
    raise ValueError(f"Unknown db_type: {db_type}")

if not os.path.exists(input_path):
    raise FileNotFoundError(f"Input file not found: {input_path}")
if not os.path.exists(db_file):
    raise FileNotFoundError(f"DB file not found: {db_file}")

# Read original file and extract the ID column (6th column, index 5)
orig_df = pd.read_csv(input_path, sep='\t', header=0, dtype=str, na_values = 'NAN') #, low_memory=False)
# orig_df.columns=[ 'Chr', 'Pos', 'End', 'Ref', 'Alt', 'ID']
if orig_df.shape[1] <= 5:
    raise ValueError("Input file has fewer than 6 columns; cannot extract ID from column index 5")
if db in orig_df.columns:
    orig_df = orig_df.drop(columns=[db])

id_series = orig_df['ID'].astype(str).str.strip()
# id_df = pd.DataFrame({'ID': id_series})
if os.path.getsize(db_file) == 0:
    orig_df[db] = 'NAN'
else:
    # Heuristic to pick ID and annotation columns:
    # original script used db_df.iloc[:, [7,1]] for generic/region; keep that for those types
    if db_type in ('generic', 'region'):
        print(db_file)
        # Check if db_file is empty
        # Read db file and extract two columns (ID and annotation)
        with open(db_file, 'rb') as f:
            try:
                enc = chardet.detect(f.read())  # read first 10000 bytes
            except Exception as e:
                print(f"Warning: Error detecting encoding: {e}")
                enc = {'encoding': 'utf-8'}  # fallback to utf-8
            encoding = enc['encoding']
            db_df_raw = pd.read_csv(db_file, sep='\t', header=None, dtype=str, na_values = 'NAN', encoding=encoding,  low_memory=True)

    # if db_df_raw.shape[1] <= 7:
    #     raise ValueError(f"DB file {db_file} does not have expected columns (needs index 7)")
        db_df = db_df_raw.iloc[:, [7, 1]].copy()
        db_df.columns = ['ID', db]
        # Create mapping and add column to original dataframe
        mapping = db_df.set_index('ID')[db]
        orig_df[db] = id_series.map(mapping).fillna('NAN')
    elif db_type == 'gene':
        print(db_file_all)        
        # gene: fall back to first column as ID and second column (if present) as annotation,
        # otherwise use the last column as annotation
        db_df_all_raw = pd.read_csv(db_file_all, sep='\t', header=None, dtype=str, na_values = 'NAN')
        db_df_exon_raw = pd.read_csv(db_file_exon, sep='\t', header=None, dtype=str, na_values = 'NAN')
        id_col = 7
        db_df_all_df = db_df_all_raw.iloc[:, [id_col, 0, 1]].copy()
        db_df_all_df.iloc[:,1] = db_df_all_df.iloc[:,[1,2]].apply(lambda x: '&'.join(x.astype(str)), axis=1)
        db_df_all_df = db_df_all_df.iloc[:, [0,1]]
        db_df_all_df.columns = ['ID', f'{db}_all_func']

        db_file_exon_df = db_df_exon_raw.iloc[:, [id_col+1, 1, 2]].copy()
        db_file_exon_df.iloc[:,1] = db_file_exon_df.iloc[:,1].str.replace(' ', '_')
        db_file_exon_df.iloc[:,1] = db_file_exon_df.iloc[:,[1,2]].apply(lambda x: '&'.join(x.fillna('NAN').astype(str)), axis=1)
        db_file_exon_df = db_file_exon_df.iloc[:, [0,1]]
        db_file_exon_df.columns = ['ID', f'{db}_exon_func']

        # Create mapping and add column to original dataframe
        mapping = db_df_all_df.set_index('ID')[f'{db}_all_func']
        orig_df[f'{db}_all_func'] = id_series.map(mapping).fillna('NAN')
        mapping = db_file_exon_df.set_index('ID')[f'{db}_exon_func']
        orig_df[f'{db}_exon_func'] = id_series.map(mapping).fillna('NAN')
    elif db_type == 'cadd':
        print(db_file)
        # Read db file and extract two columns (ID and annotation)
        db_df_raw = pd.read_csv(db_file, sep='\t', header=None, dtype=str, na_values = 'NAN')

        if db_df_raw.shape[1] <= 1:
            raise ValueError(f"DB file {db_file} does not have expected columns (needs index 2)")

        db_df = db_df_raw.copy()
        db_df.columns  = ['ID', db]
        db_df['ID'] = db_df['ID'].apply(lambda x: 'info:chr' + x if isinstance(x, str) else 'NAN')
        mapping = db_df.set_index('ID')[db]
        orig_df[db] = id_series.map(mapping).fillna('NAN')
    else:
        raise ValueError(f"Unknown db_type: {db_type}")

    # Write output (tab-separated, no index)
orig_df.to_csv(output_file, sep='\t', index=False, na_rep="NAN",  header=True)
    # ...existing code...
