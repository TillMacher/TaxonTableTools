def create_metadata_table(TaXon_table_xlsx, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

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
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    metadata_df = pd.DataFrame(samples_metadata_list)
    metadata_df = metadata_df.rename({0: 'Samples'}, axis='columns')
    metadata_df.to_excel(Meta_data_table_xlsx, index=False)

    closing_text = "Metadata table is found under:\n" + '/'.join(str(Meta_data_table_xlsx).split("/")[-4:])
    print("\n" + closing_text)
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
