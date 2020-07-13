# check the input format
def check_taXon_table_format(taXon_table):

    #taXon_table = "/Users/tillmacher/Desktop/Projects/Projects_Development/TaXon_Table_Tools/Projects/Default_project/TaXon_tables/Default_project_taxon_table.xlsx"

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
    header_prompt = ["IDs", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status"]

    if taXon_table_df_header != header_prompt:
        sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
        raise RuntimeError("Oops! Something is wrong with the header!")

    ###################################
    # B) OTUs prompt

    OTU_list = taXon_table_df['IDs'].values.tolist()

    # loop through all available OTUs
    for OTU in OTU_list:

        # define a variable Error message
        ErrorMessage = "Oops! Something is wrong with the OTU ID: " + OTU

        try:
            # split the OTU, which should be containing of two elements
            p1 = OTU.split("_")[0]
            p2 = OTU.split("_")[1]

        except:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

        # check if the second part is an integer
        # convert p2 to an integer
        try:
            p2 = int(p2)
        except :
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # C) Species prompt
    # does not necessarily need to include both the Genus and epithet, but it is recommended

    species_list = taXon_table_df['Species'].values.tolist()
    epithet_only_list = []

    for species in species_list:
        if (species != 'nan' and len(species.split()) < 2):
            epithet_only_list.append(species)

    if epithet_only_list != []:
        WarningMessage = "Warning: There are " + str(len(epithet_only_list)) +  " species that do not fit the binomial nomenclature." + "\n" "\n" + "This error message will be ignored, but it's recommended to use binomial nomenclature!"
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

    for OTU_reads in OTU_reads_list:
        read_numbers = OTU_reads[10:]
        for read_number in read_numbers:
            try:
                read_number = int(read_number)
            except:
                OTU = OTU_reads[0]
                ErrorMessage = "Please check your read numbers in " + OTU + " -> " + str(read_number)
                sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
                raise RuntimeError(ErrorMessage)


    ###################################
    # Wrap up

    sg.Popup("Your file looks great and is ready to use!", title="Taxonomy table check", keep_on_top=True)
