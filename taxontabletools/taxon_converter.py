def taxon_converter(read_table_xlsx, taxonomy_results_xlsx, TaXon_table_name, sheet_name, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    taxonomy_results_xlsx =  Path(taxonomy_results_xlsx)
    read_table_xlsx = Path(read_table_xlsx)
    # create filename and path for output file
    Output_name = TaXon_table_name + ".xlsx"
    Output_file = path_to_outdirs / "TaXon_tables" / Output_name

    # store the file name for later use
    file_name = taxonomy_results_xlsx.name

    # create datafrmes for both files
    taxonomy_df = pd.read_excel(taxonomy_results_xlsx, sheet_name, header=0).set_index('IDs')
    if sheet_name == "BOLDigger hit":
        taxonomy_df = taxonomy_df.drop(columns=['Flags'])
    Read_df = pd.read_excel(read_table_xlsx, header=0).set_index('IDs')
    # extract OTU names
    OTU_names_taxonomy_df = taxonomy_df.index.values
    OTU_names_Read_df = Read_df.index.values
    # create a sorting column -> important to match the input files
    taxonomy_df_sorting_col = [OTU.replace("OTU_", "") for OTU in OTU_names_taxonomy_df]
    Read_df_sorting_col = [OTU.replace("OTU_", "") for OTU in OTU_names_Read_df]
    taxonomy_df["sort"] = taxonomy_df_sorting_col
    taxonomy_df = taxonomy_df.fillna("")
    Read_df["sort"] = Read_df_sorting_col
    # sort both dataframes by the sorting column
    taxonomy_df = taxonomy_df.sort_values(by=['sort'])
    Read_df = Read_df.sort_values(by=['sort'])

    # check if the sorting columns match each other -> if not break the and give an error message
    if taxonomy_df["sort"].values.tolist() != Read_df["sort"].values.tolist():
        error_message = "Fatal crash: Your files do not match!"
        sg.PopupError(error_message, title="Error")
        raise Exception()

    # create an empty list that will be used to create the output dataFrame
    dataframe_list = []
    # create a list of available OTUs (does not matter from which dataframe, because they have to be the same at this point)
    OTU_list = taxonomy_df.index.values

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(OTU_list) + 1
    ############################################################################

    ############################################################################
    event, values = window_progress_bar.read(timeout=10)
    if event == 'Cancel'  or event is None:
        window_progress_bar.Close()
        raise RuntimeError
    # update bar with loop value +1 so that bar eventually reaches the maximum
    progress_update += 0
    progress_bar.UpdateBar(progress_update)
    ############################################################################

    # iterate through the OTUs (=row of dataframe) and increase the counter with each step for visual updates
    for i in range(len(OTU_list)):
        # get the name of the OTU
        OTU_name = OTU_list[i]
        # get the taxonomy from the taxonomy file
        taxonomy = taxonomy_df.iloc[i].values.tolist()
        if len(taxonomy[5].split()) == 1:
            taxonomy[5] = taxonomy[4] + " " + taxonomy[5]
        # get all information from the Read table
        read_table = Read_df.iloc[i].values.tolist()
        # get the sequence information from the Read table
        sequence = read_table[-2]
        # get the read numbers from the Read table
        read_table = read_table[:-2]
        # get the sorting row
        sort = taxonomy_df["sort"][i]
        # append all collected information to a list from which the final dataframe will be created
        dataframe_list.append([OTU_name.replace(">", "")] + taxonomy[:-1] + sequence.split() + read_table + [int(sort)])

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    # create the column headers for the final file
    # get the names of the taxonomyigger file
    taxonomy_columns = list(taxonomy_df.columns)[:-1]
    # get the names of the Read table
    sample_names_columns = list(Read_df.columns)[:-2]
    # create the column headers
    header_row = ["IDs"] + taxonomy_columns + ["seq"] + sample_names_columns + ["sort"]

    # create a new dataframe with the just created headers
    df = pd.DataFrame(dataframe_list)
    df.columns = header_row
    # set the index to be the OTUs
    df.set_index('IDs')
    # sort the dataframe by the sorting column
    df = df.sort_values(by=["sort"])
    # remove the sorting column as it is no longer needed
    df = df.drop(["sort"], axis=1)
    # save the newly created Taxon table in TaXon format as excel file
    df.to_excel(Output_file, sheet_name='TaXon table', index=False)

    closing_text = "Taxon table is found under:\n" + '/'.join(str(Output_file).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    input = taxonomy_results_xlsx.name + " + " + read_table_xlsx.name
    ttt_log("taXon table converter", "processing", input, Output_file.name, "nan", path_to_outdirs)
