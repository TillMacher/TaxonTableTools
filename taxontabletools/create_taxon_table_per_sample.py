import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path

def create_taxon_table_per_sample(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    TaXon_datasheet = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    TaXon_table = TaXon_datasheet.values.tolist()
    samples_to_process = TaXon_datasheet.columns[10:]
    first_ten_columns_header = TaXon_datasheet.columns[:10].values.tolist()
    first_ten_columns = TaXon_datasheet.iloc[:,[0,1,2,3,4,5,6,7,8,9]].values.tolist()
    OTU_list = TaXon_datasheet['ID'].values.tolist()

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_to_process) + 1
    ############################################################################

    for sample in samples_to_process:

        Output_name = Path(sample + ".xlsx")
        Output_file = path_to_outdirs / "TaXon_tables_per_sample" / Output_name

        read_numbers = TaXon_datasheet[sample].values.tolist()
        sample_rows_list = []

        for i, read_number in enumerate(read_numbers):
            if read_number > 0:
                sample_rows_list.append(first_ten_columns[i] + [read_number])

        headers_df = pd.DataFrame([first_ten_columns_header + [sample]])
        sample_df = pd.DataFrame(sample_rows_list)
        sample_df = headers_df.append(sample_df)

        sample_df.to_excel(Output_file, engine='xlsxwriter', sheet_name='TaXon table', index=False, header=False)

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

    closing_text = "\n" + "Taxon tables are found in: " + str(path_to_outdirs) + "/TaXon_tables_per_sample/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    placeholder = TaXon_table_file.name + " (multiple files)"
    ttt_log("taXon table per sample", "analysis", TaXon_table_file.name, placeholder, "nan", path_to_outdirs)
