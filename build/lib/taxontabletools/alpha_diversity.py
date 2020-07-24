def alpha_diversity(TaXon_table_xlsx, meta_data_to_test, alpha_w, alpha_h, alpha_s, alpha_font, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path
    from matplotlib.pyplot import cm

    # Meta_data_table_xlsx = Path("/home/till/Desktop/Projects/Projects_Development/TaXon_Table_Tools/Projects/JDS/Meta_data_table/JDS_taxon_table_no_NC_metadata.xlsx")
    # TaXon_table_xlsx = Path("/home/till/Desktop/Projects/Projects_Development/TaXon_Table_Tools/Projects/JDS/TaXon_tables/JDS_taxon_table_no_NC.xlsx")
    # meta_data_to_test = "km"
    # alpha_diversity_estimate = "observed OTUs"
    # path_to_outdirs = "/home/till/Desktop/Projects/Projects_Development/TaXon_Table_Tools/Projects/Default_project"

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    Output_basic_alpha_plot_pdf = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + ".pdf")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        samples = Meta_data_table_samples
        OTU_abundances_dict = {}
        samples_metadata_list = []

        # remove samples that do not fit the format
        for i, sample in enumerate(samples):
            meta_data = str(Meta_data_table_df.loc[i][meta_data_to_test])
            samples_metadata_list.append([meta_data])

        #################################
        # Calculate Alpha diversity measurements (= observed_otus)

        observed_otus_dict = {}

        for sample in samples:
            observed_otus = len([taxon for taxon in TaXon_table_df[sample].values.tolist() if taxon != 0])
            observed_otus_dict[sample] = observed_otus

        # create list which contains the site for each sample
        categories, color_groups_dict = [], {}

        categories = Meta_data_table_df[meta_data_to_test].values.tolist()

        categories_set = set(categories)

        colors = [plt.cm.tab10(i/float(len(categories_set)-1)) for i in range(len(categories_set))]

        for i, category in enumerate(categories_set):
            color_groups_dict[category] = colors[i]

        ########################################
        # create the plot

        fig, ax = plt.subplots(figsize=(int(alpha_w), int(alpha_h)))
        for i, sample in enumerate(samples):
            group = samples_metadata_list[i][0]
            try:
                group_color = color_groups_dict[group]
            except:
                group_color = color_groups_dict[int(group)]
            plt.scatter(sample, observed_otus_dict[sample], s=int(alpha_s), label=group, c=[group_color])
        plt.xticks(fontsize=alpha_font, rotation=90)
        plt.yticks(fontsize=alpha_font)
        plt.title("Observed OTUs", fontsize=alpha_font)
        plt.ylabel("# OTUs", fontsize=alpha_font, rotation=90)
        ax.set_axisbelow(True)
        ax.grid(color='gray', alpha=0.1)

        # add a legend
        handles, labels = ax.get_legend_handles_labels()
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
        plt.legend(*zip(*unique), loc='center left', bbox_to_anchor=(1, 0.5))

        plt.draw()
        plt.pause(0.001)
        answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
        if answer == "Yes":
            plt.savefig(Output_basic_alpha_plot_pdf, bbox_inches='tight')
            plt.close()
            print("\n" + "Alpha diversity estimate plots are found in", path_to_outdirs, "/Alpha_diversity/")
            sg.Popup("Alpha diversity estimate are found in", path_to_outdirs, "/Alpha_diversity/", title="Finished", keep_on_top=True)
        else:
            plt.close()

    else:
        sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)
        print("Error: The samples between the taxon table and meta table do not match!")

    # except:
    #     sg.PopupError("Something went wrong! Make sure that the metadata is varying between samples!", keep_on_top=True)
    #     print("Something went wrong! Make sure that the metadata is varying between samples!")
