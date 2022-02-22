import pandas as pd

def reduce_taxontable(df, taxonomic_level, samples, meta_data_to_test):

    if meta_data_to_test == '':
        ## extract the relevant data
        df = df[[taxonomic_level] + samples]
        ## define an aggregation function to combine multiple hit of one taxonimic level
        aggregation_functions = {}
        ## define samples functions
        for sample in samples:
            ## 'sum' will calculate the sum of p/a data
            aggregation_functions[sample] = 'sum'
        ## define taxon level function
        aggregation_functions[taxonomic_level] = 'first'
        ## create condensed dataframe
        df = df.groupby(df[taxonomic_level]).aggregate(aggregation_functions)
        if '' in df.index:
            df = df.drop('')

        return df

    else:
        ## collect metadata for each taxon
        taxon_metadata_dict = {}
        for row in df[[taxonomic_level, meta_data_to_test]].values.tolist():
            taxon = row[0]
            metadata = row[1]
            taxon_metadata_dict[taxon] = metadata

        ## extract the relevant data
        df = df[[taxonomic_level] + samples]
        ## define an aggregation function to combine multiple hit of one taxonimic level
        aggregation_functions = {}
        ## define samples functions
        for sample in samples:
            ## 'sum' will calculate the sum of p/a data
            aggregation_functions[sample] = 'sum'
        ## define taxon level function
        aggregation_functions[taxonomic_level] = 'first'
        ## create condensed dataframe
        df = df.groupby(df[taxonomic_level]).aggregate(aggregation_functions)
        if '' in df.index:
            df = df.drop('')

        ## add metadata back to taxa
        metadata_col = []
        for taxon in df[taxonomic_level].values.tolist():
            metadata_col.append(taxon_metadata_dict[taxon])
        df[meta_data_to_test] = metadata_col

        return df

def strip_metadata(df):
    " Strip the additional metadata from the dataframe to produce a clean TaxonTable "

    ## position of seq == acts as separator
    seq_loc = df.columns.get_loc("seq")

    ## samples are always assending the seq column
    samples = df.columns.tolist()[seq_loc+1:]

    ## standard columns
    standard_columns = ['ID', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Similarity', 'Status', 'seq']

    ## extract legacy TTT frame
    df = df[standard_columns + samples]

    return df

def collect_metadata(df):
    " Collect the additional metadata from the dataframe to be able to place it back later "

    ## position of status and seq == acts as separator
    seq_loc = df.columns.get_loc("seq")
    status = df.columns.get_loc("Status")

    ## metadata is always located between status and seq
    metadata = df.columns.tolist()[status+1:seq_loc]

    ## collect IDs and metadata
    df = df[['ID'] + metadata]

    return df

def add_metadata(df, metadata_df):
    " Collect the additional metadata from the dataframe to add them back to a df "

    for metadata in  metadata_df.columns.tolist()[1:]:
        status_loc = df.columns.get_loc("Status") + 1
        ## loope through IDs and find all that still remain in the altered Taxon table (e.g. after filtering less OTUs might be present)
        sub_df = pd.DataFrame([i for i in metadata_df[['ID', metadata]].values.tolist() if i[0] in df['ID'].values.tolist()], columns = ['ID', metadata])
        ## insert sub_df containing the respective metadata after the Status column
        df.insert(status_loc, str(metadata), sub_df[metadata].values.tolist(), True)

    return df

#
