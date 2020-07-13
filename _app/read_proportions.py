def read_proportions(TaXon_table_xlsx, taxonomic_level, path_to_outdirs, x_read_props_scatter, y_read_props_scatter, font_read_props_scatter, no_hit_visualization, read_props_alternating_colors):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)
    samples_list = TaXon_table_df.columns.tolist()[10:]
    Species_read_proportion_dict = {}

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[samples_list].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        print("Please do not use presence absence data!")
        sg.Popup("Please do not use presence absence data!", title=("Error"))
        raise RuntimeError

    ##############################################################################
    # scatter plot

    read_proportion_plot_directory = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem)
    scatter_pdf = Path(str(read_proportion_plot_directory) + "_" + taxonomic_level + "_scatter.pdf")

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_list) + 1
    ############################################################################

    for sample in samples_list:

        # check if the option "best hit" is chosen
        if no_hit_visualization == "best_hit":
            TaXon_table_list = []

            taxon_levels_dict = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
            no_taxonomic_level = taxon_levels_dict[taxonomic_level]
            df = TaXon_table_df[['IDs', "Phylum", "Class", "Order", "Family", "Genus", "Species", sample]]
            df_2 = df[[sample]]/df[[sample]].sum()
            df = df.assign(perc=df_2.values)

            # loop through the OTU dataframe and create new dict to store the taxonomy and read counts for each OTU
            for OTU in df.values.tolist():
                OTU_dict = {}
                # if the taxonomic level that was chosen by the user is "nan" look for the best hit
                if OTU[no_taxonomic_level] == "nan":
                    test_taxonomy_list = []
                    # loop through all taxonomic levels and check were the first level was identified (== is not "nan")
                    # append the best hit to a list
                    for test in range(1,no_taxonomic_level-1):
                        if OTU[test] != "nan":
                            test_taxonomy_list.append([test, OTU[test]])
                    # collect the OTU ID, the taxon name and the read proportions
                    TaXon_table_list.append([OTU[0], OTU[test_taxonomy_list[-1][0]] + "*", OTU[-1]])
                # skip it the OTU was identified on the correct level
                else:
                    TaXon_table_list.append([OTU[0], OTU[no_taxonomic_level], OTU[-1]])

        else:
            df = TaXon_table_df[['IDs', taxonomic_level, sample]]
            df_2 = df[[sample]]/df[[sample]].sum()
            df = df.assign(perc=df_2.values)
            df = df.drop(columns=[sample])
            TaXon_table_list = df.values.tolist()

        for line in TaXon_table_list:

            species = line[1]

            reads = line[2] * 100
            dict_key = sample, species

            if dict_key not in Species_read_proportion_dict.keys():
                Species_read_proportion_dict[dict_key] = reads
            else:
                Species_read_proportion_dict[dict_key] = Species_read_proportion_dict[dict_key] + reads

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    sample_list = []
    species_list = []
    read_proportion_list = []

    for dict_key, value in list(Species_read_proportion_dict.items()):
        sample_list.append(dict_key[0])
        species_list.append(dict_key[1])
        read_proportion_list.append(value)

    # create data
    x = sample_list #['a', 'b', 'c']
    y = species_list #['d', 'f', 'g']-
    z = read_proportion_list #[199, 290, 290]

    # calculate the number of taxa
    n_taxa = len(set(species_list))
    # calculate the number of scatter entries for the plot
    n_species_list = len(species_list)

    # create an even color list if the number of taxa is even
    if (n_taxa % 2) == 0:
       clist = (["#52af2fff", "#4472c4ff"] * int(n_species_list / 2))
    # create an odd color list if the number of taxa is odd
    # this every sample will start with the same color
    else:
       clist = (["#52af2fff", "#4472c4ff"] * (int((n_taxa - 1) / 2)) + ["#52af2fff"]) * n_species_list
       clist = clist[:n_species_list]

    title = "Read proportions (" + taxonomic_level + ")"
    fig, ax = plt.subplots(figsize=(int(x_read_props_scatter), int(y_read_props_scatter)))
    if read_props_alternating_colors == True:
        scatter = ax.scatter(x, y, s=z, color=clist, alpha=0.9)
    else:
        scatter = ax.scatter(x, y, s=z, alpha=0.9)
    plt.xticks(fontsize=font_read_props_scatter, rotation=90)
    plt.yticks(fontsize=font_read_props_scatter)
    plt.title(title, fontsize=8)
    ax.set_axisbelow(True)
    ax.grid(color='gray', alpha=0.1)

    # produce a legend with a cross section of sizes from the scatter
    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.1, color="blue")
    legend2 = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), title="Proportion(%)")

    ########################

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        print("Saving...")
        plt.savefig(scatter_pdf)
        plt.close()
        closing_text = "Read proportion plot is found under:\n" + '/'.join(str(scatter_pdf).split("/")[-4:])
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()
