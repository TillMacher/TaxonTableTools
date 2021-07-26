def taxon_table_converter_ttt(read_table_xlsx, taxonomy_results_xlsx, TaXon_table_name, sheet_name, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    ## collect both input files
    taxonomy_results_xlsx =  Path(taxonomy_results_xlsx)
    read_table_xlsx = Path(read_table_xlsx)

    #ä create filename and path for output file
    Output_name = TaXon_table_name + ".xlsx"
    Output_file = path_to_outdirs / "TaXon_tables" / Output_name

    ## store the file name for later use
    file_name = taxonomy_results_xlsx.name

    ## create dataframes for both files
    taxonomy_df = pd.read_excel(taxonomy_results_xlsx, sheet_name, header=0)
    if sheet_name == "BOLDigger hit":
        taxonomy_df = taxonomy_df.drop(columns=['Flags'])
    read_table_df = pd.read_excel(read_table_xlsx, header=0)

    ## create a new dataframe
    TaXon_table_df = taxonomy_df

    # check if all OTU are correctly sorted and present in both files
    if taxonomy_df["ID"].to_list() == read_table_df["ID"].to_list():

        ## append the sequences to the TaXon stable
        TaXon_table_df["seq"] = read_table_df["Sequences"]

        ## remove the sequence column from the read table
        read_table_df.drop('Sequences', axis='columns', inplace=True)

        ## remove the ID column from the read table
        read_table_df.drop('ID', axis='columns', inplace=True)

        ## add samples to the dataframe
        TaXon_table_df = pd.concat([TaXon_table_df, read_table_df], axis=1)

        ## check if species are present as "Genus" + "Epithet"
        new_species_column = []
        for OTU in TaXon_table_df[["Genus", "Species"]].fillna("nan").values.tolist():
            if (OTU != ["nan", "nan"] and OTU[1] != 'nan'):
                if OTU[0] not in OTU[1]:
                    new_species_column.append(OTU[0] + " " + OTU[1])
                else:
                    new_species_column.append(OTU[1])
            else:
                new_species_column.append("")

        ## add new species column to the dataframe
        TaXon_table_df["Species"] = new_species_column

        ## save the newly created Taxon table in TaXon format as excel file
        TaXon_table_df.to_excel(Output_file, sheet_name='TaXon table', index=False)

        closing_text = "Taxon table is found under:\n" + '/'.join(str(Output_file).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        input = taxonomy_results_xlsx.name + " + " + read_table_xlsx.name
        ttt_log("taXon table converter", "processing", input, Output_file.name, "ttt", path_to_outdirs)

    else:
        sg.PopupError("Error: The IDs of the read table and taxonomy table do not match!")

def taxon_table_converter_qiime2(read_table_tsv, taxonomy_results_xlsx, TaXon_table_name, sheet_name, path_to_outdirs):

    # read_table_tsv = "/Users/tillmacher/Downloads/tutorial_read_table_qiime2.tsv"
    # taxonomy_results_xlsx = "/Users/tillmacher/Downloads/tutorial_taxonomy_table.xlsx"

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    taxonomy_results_xlsx =  Path(taxonomy_results_xlsx)
    read_table_tsv = Path(read_table_tsv)

    # create filename and path for output file
    Output_name = TaXon_table_name + ".xlsx"
    Output_file = path_to_outdirs / "TaXon_tables" / Output_name

    # store the file name for later use
    file_name = taxonomy_results_xlsx.name

    # create datafrmes for both files
    taxonomy_df = pd.read_excel(taxonomy_results_xlsx, sheet_name, header=0)
    if sheet_name == "BOLDigger hit":
        taxonomy_df = taxonomy_df.drop(columns=['Flags'])

    read_table_df = pd.read_csv(Path(read_table_tsv), sep="\t")

    # drop the first row
    read_table_df = read_table_df.iloc[1:]
    read_table_df = read_table_df.reset_index(drop=True)

    ## create a new dataframe
    TaXon_table_df = taxonomy_df

    # check if all OTU are correctly sorted and present in both files
    if taxonomy_df["ID"].to_list() == read_table_df["id"].to_list():

        ## append the sequences to the TaXon stable
        TaXon_table_df["seq"] = read_table_df["Sequence"].values.tolist()

        ## remove the sequence column from the read table
        read_table_df.drop('Sequence', axis='columns', inplace=True)

        ## remove the ID column from the read table
        read_table_df.drop('id', axis='columns', inplace=True)

        ## add samples to the dataframe
        TaXon_table_df = pd.concat([TaXon_table_df, read_table_df], axis=1)

        ## check if species are present as "Genus" + "Epithet"
        new_species_column = []
        for OTU in TaXon_table_df[["Genus", "Species"]].fillna("nan").values.tolist():
            if (OTU != ["nan", "nan"] and OTU[1] != 'nan'):
                if OTU[0] not in OTU[1]:
                    new_species_column.append(OTU[0] + " " + OTU[1])
                else:
                    new_species_column.append(OTU[1])
            else:
                new_species_column.append("")

        ## add new species column to the dataframe
        TaXon_table_df["Species"] = new_species_column

        ## save the newly created Taxon table in TaXon format as excel file
        TaXon_table_df.to_excel(Output_file, sheet_name='TaXon table', index=False)

        closing_text = "Taxon table is found under:\n" + '/'.join(str(Output_file).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        input = taxonomy_results_xlsx.name + " + " + read_table_tsv.name
        ttt_log("taXon table converter", "processing", input, Output_file.name, "qiime2", path_to_outdirs)

    else:
        sg.PopupError("Error: The IDs of the read table and taxonomy table do not match!")
