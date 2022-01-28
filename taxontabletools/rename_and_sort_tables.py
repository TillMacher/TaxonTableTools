import pandas as pd
from pathlib import Path
import PySimpleGUI as sg
import numpy as np

def rename_samples(TaXon_table_xlsx, meta_data_to_test, path_to_outdirs):

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    new_names = Meta_data_table_df[meta_data_to_test].values.tolist()

    if len(set(new_names)) != len(new_names):
        sg.PopupError("Error: Found duplicates in the new sample names!", keep_on_top=True)
    else:
        answer = sg.PopupOKCancel("Warning: This will overwrite your Taxon table and rename all your samples!\n\nProceed anyways?")

        if answer == "OK":
            if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):
                if len(set(Meta_data_table_samples)) == len(Meta_data_table_samples):
                    ## collect new names and rename the dataframe
                    columns_list = TaXon_table_df.columns.tolist()
                    ## relabel samples according to their index
                    for sample, rename in Meta_data_table_df[["Samples", meta_data_to_test]].values.tolist():
                        i = columns_list.index(sample)
                        columns_list[i] = rename
                    # save the renamed taxon table
                    TaXon_table_df.columns = columns_list
                    Meta_data_table_df["Old names"] = Meta_data_table_samples
                    Meta_data_table_df["Samples"] = new_names
                    TaXon_table_df.to_excel(TaXon_table_xlsx, index=False, sheet_name="TaXon table")
                    Meta_data_table_df.to_excel(Meta_data_table_xlsx, index=False)
                    sg.Popup("Samples were successfully renamed according to the metadata table!")
                else:
                    sg.PopupError("Error: Found duplicates in the sample names!", keep_on_top=True)
            else:
                sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)

def sort_samples(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    answer = sg.PopupOKCancel("Warning: This will overwrite your Taxon table and rename all your samples!\n\nProceed anyways?")

    if answer == "OK":
        if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):
            if len(set(TaXon_table_samples)) == len(TaXon_table_samples):
                new_df = pd.DataFrame() #
                for sample in Meta_data_table_samples:
                    new_df[sample] = TaXon_table_df[sample].values.tolist()
                df1 = TaXon_table_df[TaXon_table_df.columns.tolist()[:10]]
                result = pd.concat([df1, new_df], axis=1, sort=False)
                result.to_excel(TaXon_table_xlsx, index=False)
                sg.Popup("Taxon table was successfully sorted according to the metadata table!")
            else:
                sg.PopupError("Error: Found duplicates in the sample names!", keep_on_top=True)
        else:
            sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)












#
