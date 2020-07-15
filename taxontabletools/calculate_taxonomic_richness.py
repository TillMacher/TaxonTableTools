def calculate_taxonomic_richness(TaXon_table_xlsx, path_to_outdirs, x_tax_rich, y_tax_rich, font_tax_rich, ylim_tax_rich):

    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    import matplotlib.pyplot as plt
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
    number_of_taxa_per_level = statistics_set

    # create plot
    output_file = Path(str(path_to_outdirs) + "/" + "Taxonomic_richness_plots" + "/" + TaXon_table_file.stem + "_taxonomic_richness_bar.pdf")

    plt.figure(figsize=(int(x_tax_rich), int(y_tax_rich)))
    if ylim_tax_rich != '':
        plt.ylim(0, int(ylim_tax_rich))
    plt.bar(taxon_levels, number_of_taxa_per_level)
    plt.ylabel('# taxa')
    plt.title('Number of taxa per taxonomic level')
    for i, v in enumerate(number_of_taxa_per_level):
        plt.text(i - 0.01, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_rich)

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()

        closing_text = "\n" + "Taxonomic richness plots are found in: " + str(path_to_outdirs) + "/taxonomic_richness_plots/"
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()
