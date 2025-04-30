import polars as pl
import os
from tqdm import tqdm

def iterative_update(df, mapping_files, col, on_col):
    # Assume df is your main DataFrame and mapping_files is a list of your mapping file paths
    assert col in ['pos', 'pos19', 'SNP']
    df = df.with_columns(pl.lit(None).alias(col).cast(pl.Int64))
    def update_col(df, mapping_df, col, on_col):
        if 'pos19' in mapping_df.columns:
            mapping_df = mapping_df.with_columns(pl.col('pos19').cast(pl.Int64).alias('pos19'))
        if 'pos' in mapping_df.columns:
            mapping_df = mapping_df.with_columns(pl.col('pos').cast(pl.Int64).alias('pos'))
        if 'chr' in mapping_df.columns:
            mapping_df = mapping_df.with_columns(pl.col('chr').cast(pl.Int32).alias('chr'))
        # df = df.with_columns(
        #     (pl.min_horizontal("REF", "ALT") + pl.max_horizontal("REF", "ALT")).alias("allele_cat")
        # )
        # Select rows where POS is missing
        missing_pos_df = df.filter(pl.col(col).is_null())
        # If there are no missing POS values, return the original DataFrame
        if missing_pos_df.height == 0:
            return df, True
        df_no_missing = df.filter(pl.col(col).is_not_null())
        mapping_df = mapping_df.select([col] + on_col)
        # Perform the join only on the subset with missing POS
        # Since all POS in mapping are not null, no need to filter the mapping_df
        update_df = missing_pos_df.drop(col).join(mapping_df, on=on_col, how='left').select(df.columns)
        # Update the original DataFrame
        # Use an anti join to find rows in the original df not in the subset with missing POS
        # Concatenate the non-missing part with the updated missing part
        updated_df = pl.concat([df_no_missing, update_df])
        finished = update_df.filter(pl.col('pos').is_null()).height == 0
        return updated_df, finished
    finished = False
    for mapping_file in tqdm(mapping_files):
        print(f'mapping using {mapping_file}')
        if finished:
            break
        if mapping_file.endswith('.parquet'):
            mapping_df = pl.read_parquet(mapping_file)  # Load the mapping file
        elif mapping_file.endswith('.tsv'):
            mapping_df = pl.read_csv(mapping_file, separator='\t') 
        elif mapping_file.endswith('.csv'):
            mapping_df = pl.read_csv(mapping_file)
        else:
            print(f"couldn't open {mapping_file}, skipping")
            continue
        df, finished = update_col(df, mapping_df, col, on_col)  # Update the main DataFrame
    # Now df should have the missing POS values filled in where possible
    return df



# Build
build = 37
if build == 37:
    mapping_files = ['dbSNP/filtered_snps_GRCh38/' + x for x in os.listdir('dbSNP/filtered_snps_GRCh38')]
    col = 'pos'
    on_col=['SNP', 'chr', 'allele_cat']    
elif build == 38:
    mapping_files = ['dbSNP/filtered_snps_GRCh37/' + x for x in os.listdir('dbSNP/filtered_snps_GRCh37')]
    col = 'pos19'
    on_col=['SNP', 'chr', 'allele_cat']

# Load
# df = pl.read_csv(f'data/{in_file_name}', separator='\t', null_values=['NA'], dtypes={'n': pl.Float64})
file_path = ''
df = pl.read_parquet(file_path)

#########################  File specific cleaning (to ensure column types) #########################
# SNP need to be in string
for dt in zip(df.columns, df.dtypes):
    print(dt)

# chr need to be in Int32
df = df.with_columns(pl.col('chr').cast(pl.Int32).alias('chr'))

# pos need to be in Int64
df = df.with_columns(pl.col('pos19').cast(pl.Int64).alias('pos19'))

# create allele_cat
df = df.with_columns(
            (pl.min_horizontal("effect_allele", "other_allele") + pl.max_horizontal("effect_allele", "other_allele")).alias("allele_cat")
        )

# delete INDEL SNP
valid_bases = {'A', 'C', 'G', 'T'}
df = df.filter(
            (pl.col("effect_allele").is_in(valid_bases)) & (pl.col("other_allele").is_in(valid_bases))
        )

# num_cases/controls/samples need to be Int64
df = df.with_columns(
    pl.col("num_cases")
    .fill_nan(None)  # Convert NaN to Null
    .fill_null(0)    # Replace Nulls with 0 (or another default)
    .round(0)
    .cast(pl.Int64)
)
df = df.with_columns(
    pl.col("num_controls")
    .fill_nan(None)  # Convert NaN to Null
    .fill_null(0)    # Replace Nulls with 0 (or another default)
    .round(0)
    .cast(pl.Int64)
)
df = df.with_columns(
    pl.col("num_samples")
    .fill_nan(None)  # Convert NaN to Null
    .fill_null(0)    # Replace Nulls with 0 (or another default)
    .round(0)
    .cast(pl.Int64)
)

######################### iteratively update #########################
if __name__ == "__main__":
    df = iterative_update(df, mapping_files, col=col, on_col=on_col)
    original_pos = 'pos19' if build == 37 else 'pos'
    df = df.sort(['chr', original_pos])
    df = df.select(['SNP', 
                    'chr', 
                    'pos',
                    'pos19', 
                    'effect_allele', 
                    'other_allele', 
                    'eaf', 
                    'beta', 
                    'se', 
                    'pval', 
                    'mlogp', 
                    'num_samples', 
                    'num_cases', 
                    'num_controls'])
    df.write_parquet(file_path)