def site_occupancy(TaXon_table_xlsx, meta_data_to_test, taxonomic_level, path_to_outdirs, x_site_occ, y_site_occ, font_site_occ):

    import os
    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path
    import matplotlib.pyplot as plt
    import PySimpleGUI as sg

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header = 0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]

    TaXon_table_n_samples = len(TaXon_table_samples)

    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    answer = "Ask"
    output_message = "No"

    try:
        if os.path.exists(Meta_data_table_xlsx):

            Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header = 0)
            Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

            n_sites = len(set(Meta_data_table_df[meta_data_to_test].tolist()))

            if (sorted(TaXon_table_samples) == sorted(Meta_data_table_samples) and TaXon_table_n_samples != n_sites):

                site_occupancy_dict = {}

                sites = set(Meta_data_table_df[meta_data_to_test].tolist())

                for site in sites:
                    # this can either be a species name or the above specified taxonomic level
                    present_OTU_list = []

                    # extract samples that belong to the site from the metadata file
                    included_samples_list = Meta_data_table_df[Meta_data_table_df.values  == site]['Samples'].values.tolist()

                    # count the number of samples per site to calculate the site occupancy
                    n_samples = len(included_samples_list)

                    # create a list of all species (or the specified taxonomic level)
                    if taxonomic_level == "OTUs":
                        taxonomic_level = "IDs"
                    overall_included_species_list = TaXon_table_df[taxonomic_level].values.tolist()
                    # make the list unique
                    overall_included_species_set = set(overall_included_species_list)
                    # remove potential 'nan's from the list
                    overall_included_species_set = [x for x in overall_included_species_set if str(x) != 'nan']

                    # create a set of species that is present at the sites
                    for sample in included_samples_list:

                        OTUs_per_species_list = []

                        # check the read abundaces for each sample
                        read_abundace_list = TaXon_table_df[sample].values.tolist()

                        # enumerate the read abundaces for each sample and collect all lines that have more than one read
                        for i, read_abundance in enumerate(read_abundace_list):
                            species = TaXon_table_df[taxonomic_level][i]
                            # if reads are present, collect the species name (or the specified taxonomic level) from the TaXon table
                            if read_abundance != 0:
                                OTUs_per_species_list.append(species)

                        # remove all nans
                        OTUs_per_species_list = [x for x in OTUs_per_species_list if str(x) != 'nan']
                        # make list unique
                        OTUs_per_species_list = list(set(OTUs_per_species_list))
                        # append to list of species for the current site
                        present_OTU_list.append(OTUs_per_species_list)

                    # flatten the list of present species per site
                    present_OTU_list_flattened = [val for sublist in present_OTU_list for val in sublist]

                    # store occupancy of each species in a dict, will be accessed by position in list
                    occupancy_dict = {}

                    # count the number of occurences for each species and calculate the occpancy based on the number of samples
                    for species in overall_included_species_set:
                        count = present_OTU_list_flattened.count(species)
                        occupancy = count / n_samples
                        occupancy_dict[species] = occupancy

                    occupancy_dict = {k: v for k, v in sorted(occupancy_dict.items(), key=lambda item: item[1])}
                    occupancy_list = list(occupancy_dict.values())
                    species_list = list(occupancy_dict.keys())

                    # create a simple bar plot for each site
                    plt.figure(figsize=(int(x_site_occ), int(y_site_occ)))
                    plt.bar(species_list, occupancy_list)
                    plt.title(site + " (" + taxonomic_level + ")")
                    plt.xticks(fontsize=font_site_occ, rotation=90)
                    plt.yticks(fontsize=font_site_occ)
                    plt.ylabel("occupancy", fontsize=font_site_occ)
                    plt.ylim(0,1)

                    occupancy_plot_directory = Path(str(path_to_outdirs) + "/" + "Site_occupancy_plots" + "/" + TaXon_table_xlsx.stem)
                    if not os.path.exists(occupancy_plot_directory):
                        os.mkdir(occupancy_plot_directory)

                    occupancy_plot_pdf = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".pdf")
                    occupancy_table = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".xlsx")

                    if answer == "Ask":
                        # show the first plot
                        plt.show(block=False)
                        # ask for answer
                        answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
                        # if answer is yes, save all upcoming plots
                        if answer == "Yes":
                            output_message = "Yes"
                            occupancy_df = pd.DataFrame(occupancy_list, species_list)
                            occupancy_df.columns = ["Occupancy"]
                            occupancy_df.index.name = "Taxon"
                            occupancy_df = occupancy_df.sort_values("Occupancy")
                            # sort the table numerical if OTUs were chosen
                            if taxonomic_level == "IDs":
                                sort_list = []
                                for OTU in occupancy_df.index.tolist():
                                    sort_list.append(int(OTU.split("_")[1]))
                                occupancy_df["sort"] = sort_list
                                occupancy_df = occupancy_df.sort_values("sort")
                                occupancy_df = occupancy_df.drop("sort", axis=1)
                            occupancy_df.to_excel(occupancy_table)
                            # save the plot
                            plt.savefig(occupancy_plot_pdf)
                            plt.close()
                        else:
                            # answer is "No"
                            plt.close()
                    elif answer == "Yes":
                        # write the other sites
                            occupancy_df = pd.DataFrame(occupancy_list, species_list)
                            occupancy_df.columns = ["Occupancy"]
                            occupancy_df.index.name = "Taxon"
                            occupancy_df = occupancy_df.sort_values("Occupancy")
                            # sort the table numerical if OTUs were chosen
                            if taxonomic_level == "IDs":
                                sort_list = []
                                for OTU in occupancy_df.index.tolist():
                                    sort_list.append(int(OTU.split("_")[1]))
                                occupancy_df["sort"] = sort_list
                                occupancy_df = occupancy_df.sort_values("sort")
                                occupancy_df = occupancy_df.drop("sort", axis=1)
                            occupancy_df.to_excel(occupancy_table)
                            # save the plot
                            plt.savefig(occupancy_plot_pdf)
                            plt.close()
                    else:
                            plt.close()

                if output_message == "Yes":
                    closing_text = "Site occupancy plots are found under:\n" + '/'.join(str(occupancy_plot_pdf).split("/")[-4:])
                    print(closing_text)
                    sg.Popup(closing_text, title="Finished", keep_on_top=True)

            else:
                sg.PopupError("Please check your Metadata file and Taxon table file: The samples do not match or the metadata is unique for all samples!", keep_on_top=True)
                print("Please check your metadata file and Taxon table file: The samples do not match or the metadata is unique for all samples!")

        else:
            sg.PopupError("Missing metadata file!", keep_on_top=True)
            print("Missing metadata file!")

    except:
        exception_text = "Something went wrong!" + "\n" + "Do not use numbers as metadata keys." + "\n" + "Please check your metadata file and read the manual."
        sg.PopupError(exception_text, title="Error", keep_on_top=True)
        print(exception_text)
