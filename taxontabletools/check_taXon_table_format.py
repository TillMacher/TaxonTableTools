import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path

# check the input format
def check_taXon_table_format(taXon_table):

    try:
        taXon_table_df = pd.read_excel(Path(taXon_table), "TaXon table")
        taXon_table_df = taXon_table_df.replace(np.nan, 'nan', regex=True)
    except:
        sg.PopupError("Could not find the TaXon table sheet", keep_on_top=True)
        raise


    ###################################
    # A) header prompt
    taXon_table_df_header = taXon_table_df.columns.tolist()[0:9]
    header_prompt = ["ID", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status"]

    if taXon_table_df_header != header_prompt:
        sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
        raise RuntimeError("Oops! Something is wrong with the header!")

    ###################################
    # B) OTUs prompt

    ## removed in 1.1.6

    ###################################
    # C) Species prompt
    # does not necessarily need to include both the Genus and epithet, but it is recommended

    species_list = taXon_table_df['Species'].values.tolist()
    epithet_only_list = []

    for species in species_list:
        if (species != 'nan' and len(species.split()) < 2):
            epithet_only_list.append(species)

    if epithet_only_list != []:
        WarningMessage = "Warning: There are " + str(len(epithet_only_list)) +  " species that do not fit the binomial nomenclature." + "\n" + "This error message will be ignored, but it's recommended to use binomial nomenclature!"
        sg.Popup(WarningMessage, title="Warning", keep_on_top=True)


    ###################################
    # D) Samples prompt

    samples_list = taXon_table_df.columns.tolist()[10:]

    for sample in samples_list:
        if " " in sample:
            ErrorMessage = "Please do not use spaces in the sample names:\n" + sample
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # E) Reads prompt

    OTU_reads_list = taXon_table_df.values.tolist()

    ask = False
    answer = False
    write = False

    for OTU_reads in OTU_reads_list:
        read_numbers = OTU_reads[10:]
        if sum(read_numbers) == 0:
            ask = True
        for read_number in read_numbers:
            try:
                read_number = int(read_number)
            except:
                OTU = OTU_reads[0]
                ErrorMessage = "Please check your read numbers in " + OTU + " -> " + str(read_number)
                sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
                raise RuntimeError(ErrorMessage)

    ###################################
    # F) Taxonomy consistency

    taxonomy_list = taXon_table_df[["ID", "Phylum", "Class", "Order", "Family", "Genus", "Species"]].values.tolist()
    for entry in taxonomy_list:
        if "nan" in entry:
            taxonomy = entry[1:]
            OTU = entry[0]
            n_nan = taxonomy.count("nan")
            if 6 - taxonomy.index("nan") != n_nan:
                answer = sg.PopupOKCancel("Internally missing taxonomy found!\nReplace with placeholder?\nWarning: This will overwrite the dataframe!")
                break
    if answer == 'OK':
        new_df_list = []
        for row in taXon_table_df.values.tolist():
            entry = row[0:7]
            if "nan" in entry:
                taxonomy = entry[1:]
                OTU = entry[0]
                n_nan = taxonomy.count("nan")
                if 6 - taxonomy.index("nan") != n_nan:
                    for item in [i for i,x in enumerate(entry) if x == "nan"]:
                        row[item] = ">FLAG<"
            row = ['' if x=='nan' else x for x in row]
            new_df_list.append(row)
        column_names = taXon_table_df.columns.tolist()
        taXon_table_df_new = pd.DataFrame(new_df_list, columns=column_names)
        write = True
        sg.Popup("Replaced missing taxonomy with a >FLAG< placeholder.\nPlease adjust the taxonomy in Excel.")

    if write == True:
        taXon_table_df_new.to_excel(taXon_table, sheet_name="TaXon table", index=False)


    sg.Popup("Your file looks great and is ready to use!", title="Taxonomy table check", keep_on_top=True)
