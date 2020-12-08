def convert_to_perlodes(TaXon_table_xlsx, operational_taxon_list, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    #get the taxonomy from the operational taxon list
    operational_taxon_list_df = pd.read_excel(Path(operational_taxon_list), header=2, sheet_name="Operationelle Taxaliste")
    taxonomy_list = operational_taxon_list_df["Taxonname\n(Perlodes-Datenbank)"].values.tolist()
    taxonomy_list = [x for x in taxonomy_list if str(x) != 'nan']

    # get the according IDs from the operational taxon list
    IDs_list = operational_taxon_list_df["ID_\nART"].values.tolist()
    IDs_list = [x for x in IDs_list if str(x) != 'nan']

    # create a dict to store both the ID and the taxonomy
    operational_taxon_list_dict = {}
    for ID, taxonomy in zip(IDs_list, taxonomy_list):
       operational_taxon_list_dict[taxonomy] = int(ID)

    # load the taxon table and create a list
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_taxonomy = TaXon_table_df.columns.tolist()[0:7]
    samples_list = TaXon_table_df.columns.tolist()[10:]

    # store hits and dropped OTUs
    hit_list, dropped_list, transversion_list = [], [], []

    # loop through the taxon table
    for taxonomy in TaXon_table_df[TaXon_table_taxonomy].drop_duplicates().values.tolist():

       # collect the OTU name, species, genus and family and convert to perlodes format
       OTU = taxonomy[0]
       Species = taxonomy[6]
       Species_group = str(taxonomy[6]) + "-Gr."
       Genus = str(taxonomy[5]) + " sp."
       Family = str(taxonomy[4]) + " Gen. sp."

       # test if the OTU has a hit at: Species level, Genus level or Family level
       if Species in operational_taxon_list_dict.keys():
           # add to hit list
           hit_list.append([OTU] + [str(operational_taxon_list_dict[Species])] + [Species])
           # add to perlodes log file
           transversion_list.append(taxonomy + [str(operational_taxon_list_dict[Species])] + [Species])

       elif Species_group in operational_taxon_list_dict.keys():
           hit_list.append([OTU] + [str(operational_taxon_list_dict[Species_group])] + [Species_group])
           transversion_list.append(taxonomy + [str(operational_taxon_list_dict[Species_group])] + [Species_group])

       elif Genus in operational_taxon_list_dict.keys():
           hit_list.append([OTU] + [str(operational_taxon_list_dict[Genus])] + [Genus])
           transversion_list.append(taxonomy + [str(operational_taxon_list_dict[Genus])] + [Genus])

       elif Family in operational_taxon_list_dict.keys():
           hit_list.append([OTU] + [str(operational_taxon_list_dict[Family])] + [Family])
           transversion_list.append(taxonomy + [str(operational_taxon_list_dict[Family])] + [Family])

       # otherwise store the hit with an "nan"
       else:
           hit_list.append([OTU] + ["nan"])
           dropped_list.append(OTU)
           transversion_list.append(taxonomy + ["", ""])

    # create an output list for perlodes
    # make read abundaces binary
    perlodes_input_list = []
    for hit, row in zip(hit_list, TaXon_table_df.values.tolist()):
       # skip OTUs that were not in the OPT
       if hit[1] != "nan":
           reads_list = []
           # loop through all the sample of the file
           for reads in row[10:]:
               # reads > 0 --> 1
               if reads > 0:
                   reads_list.append(1)
               # reads == 0 --> 0
               else:
                   reads_list.append(0)
           # now append the taxonomy and the presence/absence to the perlodes list
           perlodes_input_list.append([hit[1]] + [hit[2]] + reads_list)

    # print the number of dropped OTUs
    number_of_initial_OTUs = len(TaXon_table_df[TaXon_table_taxonomy].drop_duplicates().values.tolist())
    print("Warning: Dropped " + str(len(dropped_list)) + " of " + str(number_of_initial_OTUs) + " OTUs.\n")

    # write the perlodes output file
    perlodes_df = pd.DataFrame(perlodes_input_list)
    perlodes_df.columns = ["ID_ART", "TAXON_NAME"] + samples_list

    # Perlocdes sums up the counts of all df_duplicates! So these must be removed
    perlodes_df.index = perlodes_df["TAXON_NAME"]
    # create a set of all present taxa
    taxa_set =  set(perlodes_df["TAXON_NAME"].values.tolist())
    perlodes_filtered_list = []

    # loop through all target taxa
    for taxon in taxa_set:
       # collect the ID and Taxon name
       ID_taxon_name = perlodes_df.loc[taxon][["ID_ART", "TAXON_NAME"]].values.tolist()[0]
       # if there are duplicates:
       try:
           # calculate the sum for each sample
           sum_of_pa = list(perlodes_df.loc[taxon][samples_list].sum())
           # start a new list for each taxon
           pa_list = ID_taxon_name
           # replace the count by either 1 or 0
           for pa in sum_of_pa:
               if pa > 0:
                   pa_list.append(1)
               else:
                   pa_list.append(0)
           perlodes_filtered_list.append(pa_list)
       # if there are no duplicates
       except:
           # simply take the old presence/absence line
           perlodes_filtered_list.append(perlodes_df.loc[taxon].values.tolist())

    # write the filtered list to a dataframe
    perlodes_df = pd.DataFrame(perlodes_filtered_list)
    perlodes_df.columns = ["ID_ART", "TAXON_NAME"] + samples_list
    perlodes_directory = Path(str(path_to_outdirs) + "/" + "Perlodes" + "/" + TaXon_table_xlsx.stem)
    perlodes_xlsx = Path(str(perlodes_directory) + "_perlodes.xlsx")
    perlodes_df.to_excel(perlodes_xlsx, sheet_name='ImportList', index=False)

    # write the log file file to a different dataframe
    transversion_df = pd.DataFrame(transversion_list)
    transversion_df.columns = ["IDs", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ID_ART", "TAXON_NAME"]
    transversion_xlsx = Path(str(perlodes_directory) + "_perlodes_conversion_table.xlsx")
    transversion_df.to_excel(transversion_xlsx, index=False)

    closing_text = "Perlodes input file is found under:\n" + '/'.join(str(perlodes_xlsx).split("/")[-4:])
    sg.Popup(closing_text + "\n\nWarning: Please check the converted table for errors and compare the conversion results in the conversion table.\n\nThe conversion is still in beta and might contain errors!", title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("perlodes conversion", "processing", TaXon_table_xlsx.name, perlodes_xlsx.name, "nan", path_to_outdirs)
