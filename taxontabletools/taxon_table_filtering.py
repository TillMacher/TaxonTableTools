## TaXon table filtering
#1.a create mask for user input
def create_taxa_mask(TaXon_table_xlsx, mask):

    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    data = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_taxa = set(data[mask].values.tolist())
    available_taxa = [x for x in available_taxa if str(x) != 'nan']
    available_taxa = sorted(list(available_taxa))

    return available_taxa

# 1.b filter taxa from taxon table
def taxon_table_filter(TaXon_table_xlsx, exclude_taxa, mask, appendix_name, threshold, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    # convert taxa to exclude to a list if only one taxon is given (which is then string)
    if type(exclude_taxa) == str:
        exclude_taxa = [exclude_taxa]

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    # check for taxa to filter
    mask_position = list(df.columns).index(mask)
    df_columns = df.columns

    rows_to_keep = []

    df_rows = df.values.tolist()
    for row in df_rows:
        taxon_to_evaluate = row[mask_position]
        if taxon_to_evaluate not in exclude_taxa:
            if str(taxon_to_evaluate) != 'nan':
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
        print('\n' + 'Filter theshold were to harsh: Nothing to print', '\n')

    else:
        df_out.columns = df_columns

        # write output file
        file_name = TaXon_table_file.stem
        output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")
        threshold_output = "Similarity threshold = " + str(threshold)
        exclude_taxa.append(threshold_output)
        df_exclude_taxa = pd.DataFrame(exclude_taxa)
        df_exclude_taxa.columns = ['Filter criteria']
        writer = pd.ExcelWriter(output_name, engine = 'xlsxwriter')
        df_out.to_excel(writer, sheet_name = 'TaXon table', index=False)
        df_exclude_taxa.to_excel(writer, sheet_name = 'Filter criteria', index=False)
        writer.save()
        writer.close()

        closing_text = "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        ttt_log("taxon filter", "processing", TaXon_table_file.name, output_name.name, "nan", path_to_outdirs)

# 2.a create mask for user input
def create_sample_mask(TaXon_table_xlsx):

    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_samples = df.columns.tolist()[10:]

    return available_samples

# 2.b filter samples from taxon list
def exclude_samples(TaXon_table_xlsx, exluded_samples, appendix_name, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    TaXon_table_xlsx_path = TaXon_table_xlsx
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    if type(exluded_samples) == str:
        exluded_samples = [exluded_samples]

    for sample in exluded_samples:
        df = df.drop(sample, axis=1)

    header = df.columns.values.tolist()

    row_filter_list = []

    for row in df.values.tolist():
        reads = set(row[10:])
        if reads != {0}:
            row_filter_list.append(row)

    df = pd.DataFrame(row_filter_list)
    df.columns = header

    file_name = TaXon_table_file.stem
    output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")
    df.to_excel(output_name, sheet_name = 'TaXon table', index=False)

    closing_text = "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:])
    print(closing_text)
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("sample filter", "processing", TaXon_table_file.name, output_name.name, "nan", path_to_outdirs)
