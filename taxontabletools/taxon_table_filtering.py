import pandas as pd
from pandas import DataFrame
from pathlib import Path
import PySimpleGUI as sg
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import webbrowser
from taxontabletools.taxontable_manipulation import strip_metadata
from taxontabletools.taxontable_manipulation import collect_metadata
from taxontabletools.taxontable_manipulation import add_metadata

## TaXon table filtering
#1.a create mask for user input
def create_taxa_mask(TaXon_table_xlsx, mask, taxon_mask):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    data = pd.read_excel(TaXon_table_xlsx)

    available_taxa = set(data[mask].values.tolist())
    available_taxa = [x for x in available_taxa if str(x) != 'nan']
    available_taxa = sorted(list(available_taxa))

    if taxon_mask != "":
        available_taxa = [taxon for taxon in available_taxa if taxon_mask in taxon]

    return available_taxa

# 1.b filter taxa from taxon table
def taxon_filter(TaXon_table_xlsx, filtered_taxa, mask, appendix_name, threshold, path_to_outdirs, taxon_filter_method):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)

    # convert taxa to exclude to a list if only one taxon is given (which is then string)
    if type(filtered_taxa) == str:
        filtered_taxa = [filtered_taxa]

    if taxon_filter_method == "keep":
        available_taxa = set(TaXon_table_df[mask].values.tolist())
        available_taxa = [x for x in available_taxa if str(x) != '']
        available_taxa = sorted(list(available_taxa))
        filtered_taxa = list(set(available_taxa) - set(filtered_taxa))

        # check for taxa to filter
        mask_position = list(TaXon_table_df.columns).index(mask)
        df_columns = TaXon_table_df.columns

        rows_to_keep = []

        df_rows = TaXon_table_df.values.tolist()
        for row in df_rows:
            taxon_to_evaluate = row[mask_position]
            if taxon_to_evaluate not in filtered_taxa:
                if str(taxon_to_evaluate) != '':
                    rows_to_keep.append(row)


    else:
        # check for taxa to filter
        mask_position = list(TaXon_table_df.columns).index(mask)
        df_columns = TaXon_table_df.columns

        rows_to_keep = []

        df_rows = TaXon_table_df.values.tolist()
        for row in df_rows:
            taxon_to_evaluate = row[mask_position]
            if taxon_to_evaluate not in filtered_taxa:
                rows_to_keep.append(row)

    df_out = pd.DataFrame(rows_to_keep)

    similarity_position = list(df_columns).index("Similarity")
    threshold = int(threshold)

    filtered_rows = []

    for index, row in df_out.iterrows():
        similarity = list(row)[similarity_position]
        if similarity != 'No Match':
            if int(similarity) >= threshold:
                filtered_rows.append(list(row))

    df_out = pd.DataFrame(filtered_rows)

    if df_out.empty:
        sg.PopupError('Filter theshold were to harsh: Nothing to print', title="Error", keep_on_top=True)

    else:
        df_out.columns = df_columns

        # write output file
        file_name = TaXon_table_file.stem
        output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")

        threshold_output = "Similarity threshold = " + str(threshold)

        filtered_taxa.append(threshold_output)
        df_filtered_taxa = pd.DataFrame(filtered_taxa)
        df_filtered_taxa.columns = ['Filter criteria']

        ## open write
        writer = pd.ExcelWriter(output_name, engine = 'xlsxwriter')

        ## add already existing metadata back to the df
        if len(TaXon_table_df_metadata.columns) != 1:
            df_out = add_metadata(df_out, TaXon_table_df_metadata)

        ## write dataframe
        df_out.to_excel(writer, sheet_name = 'TaXon table', index=False)
        df_filtered_taxa.to_excel(writer, sheet_name = 'Filter criteria', index=False)
        writer.save()

        ## print results for the user
        n_old_OTUs = len(TaXon_table_df["ID"].values.tolist())
        n_remaining_OTUs = len(df_out["ID"].values.tolist())
        diff_abs = n_old_OTUs - n_remaining_OTUs
        diff_rel = round(100 - n_remaining_OTUs / n_old_OTUs * 100, 2)

        ## finish script
        closing_text = "Removed " + str(diff_abs) + " OTUs (" + str(diff_rel) + "%).\n\n" + "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        log_text = str(diff_abs) + " OTUs ; " + str(diff_rel) + "%"
        ttt_log("taxon filter", "processing", TaXon_table_file.name, output_name.name, log_text, path_to_outdirs)

# 2.a create mask for user input
def create_sample_mask(TaXon_table_xlsx, sample_mask):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = strip_metadata(TaXon_table_df)

    available_samples = TaXon_table_df.columns.tolist()[10:]

    if sample_mask != "":
        available_samples = [sample for sample in available_samples if sample_mask in sample]

    return available_samples

# 2.b filter samples from taxon list
def filter_samples(TaXon_table_xlsx, selected_samples, appendix_name, path_to_outdirs, sample_filter_method):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)
    n_old_OTUs = len(TaXon_table_df["ID"].values.tolist())

    if type(selected_samples) == str:
        selected_samples = [selected_samples]

    if sample_filter_method == "exclude":
        for sample in selected_samples:
            TaXon_table_df = TaXon_table_df.drop(sample, axis=1)
    else:
        available_samples = TaXon_table_df.columns.tolist()[10:]
        for sample in available_samples:
            if sample not in selected_samples:
                TaXon_table_df = TaXon_table_df.drop(sample, axis=1)

    header = TaXon_table_df.columns.values.tolist()
    row_filter_list = []
    for row in TaXon_table_df.values.tolist():
        reads = set(row[10:])
        if reads != {0}:
            row_filter_list.append(row)

    TaXon_table_df = pd.DataFrame(row_filter_list)
    TaXon_table_df.columns = header

    file_name = TaXon_table_file.stem
    output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")

    ## add already existing metadata back to the df
    if len(TaXon_table_df_metadata.columns) != 1:
        TaXon_table_df = add_metadata(TaXon_table_df, TaXon_table_df_metadata)

    ## write file
    TaXon_table_df.to_excel(output_name, sheet_name = 'TaXon table', index=False)

    ## print results for the user
    n_remaining_OTUs = len(TaXon_table_df["ID"].values.tolist())
    diff_abs = n_old_OTUs - n_remaining_OTUs
    diff_rel = round(100 - n_remaining_OTUs / n_old_OTUs * 100, 2)

    ## finish script
    closing_text = "Removed " + str(diff_abs) + " OTUs (" + str(diff_rel) + "%).\n\n" + "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    log_text = str(diff_abs) + " OTUs ; " + str(diff_rel) + "%"
    from taxontabletools.create_log import ttt_log
    ttt_log("sample filter", "processing", TaXon_table_file.name, output_name.name, log_text, path_to_outdirs)

# 3 read-based filter
def read_filter(TaXon_table_xlsx, path_to_outdirs, read_filter_method, read_filter_treshold):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)
    samples = TaXon_table_df.columns.tolist()[10:]

    if read_filter_method == "absolute_filtering":

        ## transform dataframe to array and apply filter threshold
        a = np.array(TaXon_table_df[samples].values.tolist())
        TaXon_table_df[samples] = np.where(a < int(read_filter_treshold), 0, a).tolist()

        ## remove OTUs that have 0 reads after filtering
        row_filter_list = []
        for row in TaXon_table_df.values.tolist():
            reads = set(row[10:])
            if reads != {0}:
                row_filter_list.append(row)
        TaXon_table_df_filtered = pd.DataFrame(row_filter_list)
        TaXon_table_df_filtered.columns = TaXon_table_df.columns.tolist()

        ## save filtered dataframe to file
        file_name = TaXon_table_file.stem
        output_name = Path(str(path_to_outdirs) + "/TaXon_tables/" + file_name + "_" + read_filter_treshold + ".xlsx")

        ## add already existing metadata back to the df
        if len(TaXon_table_df_metadata.columns) != 1:
            TaXon_table_df_filtered = add_metadata(TaXon_table_df_filtered, TaXon_table_df_metadata)

        ## write file
        TaXon_table_df_filtered.to_excel(output_name, index=False)

        ## print results for the user
        n_old_OTUs = len(TaXon_table_df["ID"].values.tolist())
        n_remaining_OTUs = len(TaXon_table_df_filtered["ID"].values.tolist())
        diff_abs = n_old_OTUs - n_remaining_OTUs
        diff_rel = round(100 - n_remaining_OTUs / n_old_OTUs * 100, 2)

        ## finish script
        closing_text = "Removed " + str(diff_abs) + " OTUs (" + str(diff_rel) + "%).\n\n" + "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        log_text = str(read_filter_treshold) + " ; " + str(diff_abs) + " OTUs ; " + str(diff_rel) + "%"
        ttt_log("absolute read filter", "processing", TaXon_table_file.name, output_name.name, log_text, path_to_outdirs)

    elif read_filter_method == "relative_filtering":
        ## transform to percentage
        read_filter_rel = float(read_filter_treshold) / 100
        for sample in samples:
            ## transform to array
            a = np.array(TaXon_table_df[sample].values.tolist())
            ## calculate threshold for each sample
            sample_threshold = sum(a) * read_filter_rel
            ## apply filter to dataframe
            TaXon_table_df[sample] = np.where(a < int(sample_threshold), 0, a).tolist()

        ## remove OTUs that have 0 reads after filtering
        row_filter_list = []
        for row in TaXon_table_df.values.tolist():
            reads = set(row[10:])
            if reads != {0}:
                row_filter_list.append(row)
        TaXon_table_df_filtered = pd.DataFrame(row_filter_list)
        TaXon_table_df_filtered.columns = TaXon_table_df.columns.tolist()

        ## save filtered dataframe to file
        file_name = TaXon_table_file.stem
        output_name = Path(str(path_to_outdirs) + "/TaXon_tables/" + file_name + "_" + read_filter_treshold + ".xlsx")

        ## add already existing metadata back to the df
        if len(TaXon_table_df_metadata.columns) != 1:
            TaXon_table_df_filtered = add_metadata(TaXon_table_df_filtered, TaXon_table_df_metadata)

        ## write file
        TaXon_table_df_filtered.to_excel(output_name, index=False)

        ## print results for the user
        n_old_OTUs = len(TaXon_table_df["ID"].values.tolist())
        n_remaining_OTUs = len(TaXon_table_df_filtered["ID"].values.tolist())
        diff_abs = n_old_OTUs - n_remaining_OTUs
        diff_rel = round(100 - n_remaining_OTUs / n_old_OTUs * 100, 2)

        ## finish script
        closing_text = "Removed " + str(diff_abs) + " OTUs (" + str(diff_rel) + "%).\n\n" + "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        log_text = str(read_filter_treshold) + " ; " + str(diff_abs) + " OTUs ; " + str(diff_rel) + "%"
        ttt_log("relative read filter", "processing", TaXon_table_file.name, output_name.name, log_text, path_to_outdirs)

def subtract_NCs(TaXon_table_xlsx, path_to_outdirs, negative_controls):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)

    ##negative_controls = ["NC_3", "NC_1", "NC_2"]
    ## collect samples
    samples = [sample for sample in TaXon_table_df.columns.to_list()[10:] if sample not in negative_controls]
    ## calculate sum of NCs
    df_nc_sum = TaXon_table_df[negative_controls].sum(axis=1)
    ## create a new dataframe
    df_out = TaXon_table_df[TaXon_table_df.columns.tolist()[0:10]]

    # subtract the sum of reads found in the NCs from each OTU of the samples
    for sample in samples:
        df_out.insert(10, sample, (TaXon_table_df[sample] - df_nc_sum).values.tolist())

    ## replace negative values with 0
    num = df_out._get_numeric_data()
    num[num < 0] = 0

    ## remove empty OTUs
    out_list = [OTU for OTU in df_out.values.tolist() if sum(OTU[10:]) != 0]

    ## check if the still contains reads
    if df_out.empty:
        sg.PopupError('Filter theshold were to harsh: Nothing to print', title="Error", keep_on_top=True)

    else:
        output_xlsx = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem + "_NCsub.xlsx")
        df_out = pd.DataFrame(out_list, columns=df_out.columns.tolist()).replace("nan", "")

        ## add already existing metadata back to the df
        if len(TaXon_table_df_metadata.columns) != 1:
            df_out = add_metadata(df_out, TaXon_table_df_metadata)

        df_out.to_excel(output_xlsx, index=False)

        from taxontabletools.create_log import ttt_log
        ttt_log("nc subtract", "processing", TaXon_table_xlsx.name, output_xlsx.name, "nan", path_to_outdirs)

        ## finish script
        closing_text = str(len(TaXon_table_df) - len(df_out)) + " OTUs were removed. The Taxon table is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)



#
