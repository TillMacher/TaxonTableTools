# check the input format
def check_taXon_table_format(taXon_table):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path

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
                ErrorMessage = "Please check the taxonomy of " + OTU + ".\nInternally missing taxonomy found!"
                sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
                raise RuntimeError(ErrorMessage)

    ###################################
    # write new dataframe if user chose to remove empty OTUs
    if ask == True:
        answer = sg.PopupOKCancel("OTUs with zero reads have been detected.\nRemove OTUs from dataframe?\nWarning: This will overwrite the dataframe!")
        if answer == 'OK':
            taXon_table_list = [OTU for OTU in taXon_table_df.values.tolist() if sum(OTU[10:]) != 0]
            taXon_table_df_new = pd.DataFrame(taXon_table_list, columns=taXon_table_df.columns.tolist()).replace("nan", "")
            taXon_table_df_new.to_excel(taXon_table, sheet_name="TaXon table", index=False)



    sg.Popup("Your file looks great and is ready to use!", title="Taxonomy table check", keep_on_top=True)
