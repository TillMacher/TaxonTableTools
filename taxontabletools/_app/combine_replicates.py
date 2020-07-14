# remove PCR or Extraction Replicates
def combine_replicates(TaXon_table_xlsx, suffix_list, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    # create output file
    output_file = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + str(TaXon_table_file.stem) + "_derep.xlsx")
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    sample_names = df.columns[10:]
    unique_sample_names_list, samples_to_process_list = [], []

    for sample in sample_names:
        sample_name = sample.split("_")[0:-1]
        unique_sample_names_list.append("_".join(sample_name))

    unique_sample_names_set = sorted(set(unique_sample_names_list))

    replicates_dict ={}

    for sample in unique_sample_names_set:

        for i, suffix in enumerate(suffix_list):
            replicates_dict["rep_" + str(i)] = sample + "_" + suffix_list[i]

        combined = sample + "_comb"

        replicate_names_list = list(replicates_dict.values())

        try:
            df[combined] = df[replicate_names_list].sum(axis=1)
            df = df.drop(replicate_names_list, axis=1)
        except:
            print("Warning! No replicates found for:    " + sample)

    df.to_excel(output_file, index=False, sheet_name = 'TaXon table')
    closing_text = "Taxon table is found under:\n" + '/'.join(str(output_file).split("/")[-4:])
    print(closing_text)
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
