from Sequoia import Athena, MRPipeline, Harmoniser, PostMRProcessor
from Sequoia import functions as f
import polars as pl

wald_ratio = False
exposure_name = 'CXCR5 '

# need to provide chr
instrument  = f.create_instrument_df(SNP='rs11530564 ', chr=3, pos=118728005 , effect_allele='A', other_allele='G', file_name=exposure_name)
# instrument_2  = f.create_instrument_df(SNP='rs117437913', chr=20, pos=50554637, effect_allele='G', other_allele='A', file_name=exposure_name)
# instrument = pl.concat([instrument_1, instrument_2])

# def reverse_orientation(df: pl.DataFrame) -> pl.DataFrame:
#     df = df.with_columns((-1*pl.col('beta')).alias('beta'))
#     return df

# instrument = reverse_orientation(instrument)

outcome_df = Athena(table_list=[
                                # 'finngen',
                                # 'olink', 
                                # 'soma', 
                                # 'toxicology', 
                                # 'nightingale', 
                                # 'metabolon', 
                                # 'panukbb', 
                                # 'nightingale_ukb', 
                                # 'global_biobank',
                                'autoimmune'
                                ]).query_instrument(instrument, exposure_name=exposure_name)

# outcome_df = pl.read_csv('results/rs4680536_Athena_results.csv')
# outcome_df = outcome_df.with_columns(
#     pl.when(pl.col("effect_allele") == "G")
#     .then(pl.col("beta"))
#     .otherwise(-1*pl.col("beta"))
#     .alias("beta")
# ).with_columns(pl.lit('G').alias('effect_allele'),
#     pl.lit('A').alias('other_allele'))


# original_instrument = instrument
# instrument = instrument.with_columns(
#     pl.lit(1).alias('beta'),
#     pl.lit(1).alias('se'),
#     pl.lit(0.5).alias('pval')
# )

if wald_ratio:
    harmonised_df = Harmoniser().harmonise(exposure_df=instrument, outcome_df=outcome_df)
    # fit and save MR results
    MR_result = MRPipeline('wald_ratio').fit(harmonised_df)
    MR_result.write_csv(f'results/{exposure_name}_MR_results.csv')
    PostMRProcessor().process(MR_result)
else:
    PostMRProcessor().process(outcome_df)
