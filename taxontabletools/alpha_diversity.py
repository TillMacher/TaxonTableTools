def alpha_diversity_scatter_plot(TaXon_table_xlsx, meta_data_to_test, alpha_w, alpha_h, alpha_s, alpha_font, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path
    from matplotlib.pyplot import cm

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    Output_basic_alpha_plot_pdf = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_scatter_plot.pdf")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.PopupError("The meta data has to differ between samples!", title="Error")
        raise RuntimeError

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

def alpha_diversity_boxplot(TaXon_table_xlsx, meta_data_to_test, alpha_w, alpha_h, alpha_s, alpha_font, path_to_outdirs):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path
    from matplotlib.pyplot import cm

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    Output_basic_alpha_plot_pdf = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_boxplot.pdf")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.PopupError("The meta data has to differ between samples!", title="Error")
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        samples = Meta_data_table_samples
        OTU_abundances_dict = {}
        samples_metadata_list = []

        # remove samples that do not fit the format
        for i, sample in enumerate(samples):
            meta_data = str(Meta_data_table_df.loc[i][meta_data_to_test])
            samples_metadata_list.append(meta_data)

        #################################
        # Calculate Alpha diversity measurements (= observed_otus)

        observed_otus_dict = {}

        for i, sample in enumerate(samples):
            observed_otus = len([taxon for taxon in TaXon_table_df[sample].values.tolist() if taxon != 0])
            categorie = samples_metadata_list[i]
            if categorie not in observed_otus_dict.keys():
                observed_otus_dict[categorie] = [observed_otus]
            else:
                observed_otus_dict[categorie] = observed_otus_dict[categorie] + [observed_otus]

        # create boxplot
        labels, data = np.array(list(observed_otus_dict.keys()), dtype=object), np.array(list(observed_otus_dict.values()), dtype=object)
        plt.figure(figsize=(int(alpha_w), int(alpha_h)))
        plt.grid(color='gray', alpha=0.1)
        plt.boxplot(data,
            patch_artist=True,
            boxprops=dict(facecolor="#0000000a", color="black"),
            capprops=dict(color="black"),
            whiskerprops=dict(color="black"),
            flierprops=dict(color="black", markeredgecolor="black"),
            medianprops=dict(color="blue"),
            )
        plt.xticks(range(1, len(labels) + 1), labels, fontsize=alpha_font)
        plt.yticks(fontsize=alpha_font)
        plt.title("Observed OTUs", fontsize=alpha_font)
        plt.ylabel("# OTUs", fontsize=alpha_font, rotation=90)

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
