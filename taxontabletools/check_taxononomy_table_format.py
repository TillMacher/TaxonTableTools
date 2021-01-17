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
    header_prompt_1 = ["ID", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status"]
    header_prompt_2 = ["ID", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status", "Flags"]

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

    ## removed in 1.1.6

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
