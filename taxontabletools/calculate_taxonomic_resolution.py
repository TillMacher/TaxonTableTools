def calculate_taxonomic_resolution(TaXon_table_xlsx, path_to_outdirs, x_tax_res, y_tax_res, font_tax_res, figure_type, ylim_tax_res):

    import glob
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    taxonomic_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    statistics_list, statistics_set, statistics_dict, highest_level_dict = [], [], {}, {}
    writer = pd.ExcelWriter(TaXon_table_xlsx, engine = 'xlsxwriter')

    for taxon_to_evaluate in taxonomic_levels:
        taxa_list = [x for x in df[taxon_to_evaluate].values.tolist() if str(x) != 'nan']
        statistics = taxon_to_evaluate, len(taxa_list)
        statistics_set.append(len(set(taxa_list)))
        statistics_list.append(list(statistics))
        statistics_dict[taxon_to_evaluate] = len(taxa_list)

    highest_level_dict["Phylum"] = statistics_dict["Phylum"] - statistics_dict["Class"]
    highest_level_dict["Class"] = statistics_dict["Class"] - statistics_dict["Order"]
    highest_level_dict["Order"] = statistics_dict["Order"] - statistics_dict["Family"]
    highest_level_dict["Family"] = statistics_dict["Family"] - statistics_dict["Genus"]
    highest_level_dict["Genus"] = statistics_dict["Genus"] - statistics_dict["Species"]
    highest_level_dict["Species"] = statistics_dict["Species"]

    taxon_levels = list(highest_level_dict.keys())
    highest_level_OTUs = list(highest_level_dict.values())
    total_OTUs = list(statistics_dict.values())

    # create plot
    # option A: scatter plot
    if figure_type == "a":

        output_file = Path(str(path_to_outdirs) + "/" + "Taxonomic_resolution_plots" + "/" + TaXon_table_file.stem + "_taxonomic_resolution_scatter.pdf")

        plt.figure(figsize=(int(x_tax_res), int(y_tax_res)))
        if ylim_tax_res != '':
            plt.ylim(0, int(ylim_tax_res))
        plt.bar(taxon_levels, highest_level_OTUs)
        plt.ylabel('# OTUs')
        plt.title('Taxonomic resolution (highest taxonomic level)')
        for i, v in enumerate(highest_level_OTUs):
            plt.text(i, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_res)

        plt.show(block=False)
        answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
        if answer == "Yes":
            plt.savefig(output_file, dpi=500)
            plt.close()
            closing_text = "\n" + "Taxonomic resolution plots are found in: " + str(path_to_outdirs) + "/Taxonomic_resolution_plots/"
            print(closing_text)
            sg.Popup(closing_text, title="Finished", keep_on_top=True)

            from taxontabletools.create_log import ttt_log
            ttt_log("taxonomic resolution", "analysis", TaXon_table_file.name, output_file.name, "plot a", path_to_outdirs)

        else:
            plt.close()

    # option B: bar plot
    else:

        output_file = Path(str(path_to_outdirs) + "/" + "Taxonomic_resolution_plots" + "/" + TaXon_table_file.stem + "_taxonomic_resolution_bar.pdf")

        plt.figure(figsize=(int(x_tax_res), int(y_tax_res)))
        if ylim_tax_res != '':
            plt.ylim(0, int(ylim_tax_res))
        plt.bar(taxon_levels, total_OTUs)
        plt.ylabel('# OTUs')
        plt.title('Taxonomic resolution (highest taxonomic level)')
        for i, v in enumerate(total_OTUs):
            plt.text(i, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_res)

        plt.show(block=False)
        answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
        if answer == "Yes":
            plt.savefig(output_file, dpi=500)
            plt.close()
            closing_text = "\n" + "Taxonomic resolution plots are found in: " + str(path_to_outdirs) + "/Taxonomic_resolution_plots/"
            print(closing_text)
            sg.Popup(closing_text, title="Finished", keep_on_top=True)

            from taxontabletools.create_log import ttt_log
            ttt_log("taxonomic resolution", "analysis", TaXon_table_file.name, output_file.name, "plot b", path_to_outdirs)

        else:
            plt.close()
