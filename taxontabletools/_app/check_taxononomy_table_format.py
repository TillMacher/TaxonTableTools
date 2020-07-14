# check the input format
def check_taxononomy_table_format(taxonomy_results_xlsx, sheet_name):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path

    try:
        taxonomy_table_df = pd.read_excel(Path(taxonomy_results_xlsx), sheet_name)
        taxonomy_table_df = taxonomy_table_df.replace(np.nan, 'nan', regex=True)
    except:
        sg.PopupError("Could not find the correct sheet: " + sheet_name, keep_on_top=True)
        raise RuntimeError("Could not find the correct sheet: " + sheet_name)

    ###################################
    # A) header prompt
    taxonomy_table_df_header = taxonomy_table_df.columns.tolist()
    header_prompt_1 = ["IDs", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status"]
    header_prompt_2 = ["IDs", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status", "Flags"]

    if sheet_name == "First hit":
        if taxonomy_table_df_header != header_prompt_1:
            sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
            raise RuntimeError("Oops! Something is wrong with the header!")

    if sheet_name == "JAMP hit":
        if taxonomy_table_df_header != header_prompt_1:
            sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
            raise RuntimeError("Oops! Something is wrong with the header!")

    if sheet_name == "BOLDigger hit":
        if taxonomy_table_df_header != header_prompt_2:
            sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
            raise RuntimeError("Oops! Something is wrong with the header!")

    ###################################
    # B) OTUs prompt

    OTU_list = taxonomy_table_df['IDs'].values.tolist()

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

        # check if the first part in named "OTU"
        if p1 != "OTU":
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

        # check if the second part is an integer
        # convert p2 to an integer
        p2 = int(p2)
        # now check if it is an integer
        if type(int(p2)) is not int:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # C) Species prompt
    # does not necessarily need to include both the Genus and epithet, but it is recommended

    species_list = taxonomy_table_df['Species'].values.tolist()
    epithet_only_list = []

    for species in species_list:
        if (species != 'nan' and len(species.split()) < 2):
            epithet_only_list.append(species)

    if epithet_only_list != []:
        WarningMessage = "Warning: There are species that do not fit the binomial nomenclature:" + "\n" + "\n" + str(len(epithet_only_list)) + " OTUs \n" + "\n" + "This error message will be ignored, but it's recommended to use binomial nomenclature!"
        sg.Popup(WarningMessage, title="Warning", keep_on_top=True)

    ###################################
    # Wrap up

    sg.Popup("Your file looks great and is ready to use!", title="Taxonomy table check", keep_on_top=True)
