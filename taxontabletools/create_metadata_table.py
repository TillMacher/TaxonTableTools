def create_metadata_table(TaXon_table_xlsx, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    import sys, subprocess, os
    from pathlib import Path

    def open_table(table):
        if sys.platform == "win32":
            os.startfile(table)
        else:
            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
            subprocess.call([opener, table])

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_xslx_df = pd.read_excel(TaXon_table_xlsx)
    samples_list = TaXon_table_xslx_df.columns.tolist()[10:]
    samples_metadata_list = []

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_list) + 1
    ############################################################################

    for sample in samples_list:
        sample_metadata = []
        sample_metadata.append(sample)
        for part in sample.split("_"):
            sample_metadata.append(part)
        samples_metadata_list.append(sample_metadata)

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

    answer = "No"
    if Meta_data_table_xlsx.exists():
        answer = sg.PopupYesNo("Metadata tables already exists! Overwrite?")
        if answer == "Yes":
            metadata_df = pd.DataFrame(samples_metadata_list)
            metadata_df.columns  = ["Samples"] + ["col_" + str(column) for column in metadata_df.columns.tolist()[1:]]
            metadata_df.to_excel(Meta_data_table_xlsx, index=False)

            answer = sg.PopupYesNo("Open metadata table?", title="Finished", keep_on_top=True)
            if answer == "Yes":
                open_table(Meta_data_table_xlsx)

            from taxontabletools.create_log import ttt_log
            ttt_log("meta data table", "analysis", TaXon_table_xlsx.name, Meta_data_table_xlsx.name, "nan", path_to_outdirs)
    else:
        metadata_df = pd.DataFrame(samples_metadata_list)
        metadata_df.columns  = ["Samples"] + ["col_" + str(column) for column in metadata_df.columns.tolist()[1:]]
        metadata_df.to_excel(Meta_data_table_xlsx, index=False)

        answer = sg.PopupYesNo("Open metadata table?", title="Finished", keep_on_top=True)
        if answer == "Yes":
            open_table(Meta_data_table_xlsx)

        from taxontabletools.create_log import ttt_log
        ttt_log("meta data table", "analysis", TaXon_table_xlsx.name, Meta_data_table_xlsx.name, "nan", path_to_outdirs)

def modify_metadata_table(path_to_outdirs):
    import PySimpleGUI as sg
    from pathlib import Path
    import os, subprocess, sys

    metadata_tables_dir = Path(str(path_to_outdirs) + "/Meta_data_table")

    metadata_table = sg.PopupGetFile(initial_folder = metadata_tables_dir, message="Select a metadata table:")

    if metadata_table != None and metadata_table != '':
        if sys.platform == "win32":
            os.startfile(metadata_table)
        else:
            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
            subprocess.call([opener, metadata_table])
