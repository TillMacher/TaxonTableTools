def basic_stats(TaXon_table_xlsx, bstats_h, bstats_w, bstats_f, path_to_outdirs):

    import csv, glob, sys, os
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path
    import matplotlib.pyplot as plt

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)

    # number of samples
    n_samples = len(TaXon_table_df.columns[10:].tolist())

    # number of OTUs
    n_OTUs_total = len(TaXon_table_df['IDs'].tolist())

    # number of taxa per taxon level
    n_Phyla = list(set(TaXon_table_df['Phylum'].tolist()))
    if "nan" in n_Phyla:
        n_Phyla.remove("nan")
    n_Phyla = len(n_Phyla)

    n_Classes = list(set(TaXon_table_df['Class'].tolist()))
    if "nan" in n_Classes:
        n_Classes.remove("nan")
    n_Classes = len(n_Classes)

    n_Orders = list(set(TaXon_table_df['Order'].tolist()))
    if "nan" in n_Orders:
        n_Orders.remove("nan")
    n_Orders = len(n_Orders)

    n_Families = list(set(TaXon_table_df['Family'].tolist()))
    if "nan" in n_Families:
        n_Families.remove("nan")
    n_Families = len(n_Families)

    n_Genera = list(set(TaXon_table_df['Genus'].tolist()))
    if "nan" in n_Genera:
        n_Genera.remove("nan")
    n_Genera = len(n_Genera)

    n_Species = list(set(TaXon_table_df['Species'].tolist()))
    if "nan" in n_Species:
        n_Species.remove("nan")
    n_Species = len(n_Species)

    # number of respective status
    status_dict = {}
    status_entries_set = set(TaXon_table_df['Status'].tolist())
    for status in status_entries_set:
        count = TaXon_table_df['Status'].tolist().count(status)
        status_dict[status] = count

    # sequence lengths
    sequence_list = TaXon_table_df['seq'].tolist()
    sequence_len_list = []
    for sequence in sequence_list:
        sequence_len_list.append(len(sequence))
    sequence_len_set = set(sequence_len_list)
    min_len_seq = min(sequence_len_set)
    max_len_seq = max(sequence_len_set)
    avg_len_seq = round(sum(sequence_len_set) / len(sequence_len_set))

    # read stats per sample
    samples = TaXon_table_df.columns[10:].tolist()
    reads_dict = {}
    for sample in samples:
        # read stats
        reads_list = TaXon_table_df[sample].tolist()
        reads_sum = sum(reads_list)
        reads_avg = round(sum(reads_list) / len(reads_list))

        # OTU stats
        OTUs_list = []
        for OTU in reads_list:
            if OTU != 0:
                OTUs_list.append(OTU)
        n_OTUs = len(OTUs_list)

        # Species stats
        OTUs_species_list = TaXon_table_df[["Species", sample]].values.tolist()
        species_list = []
        for OTU in OTUs_species_list:
            if OTU[0] != 'nan' and OTU[1] != 0:
                species_list.append(OTU[0])
        n_species = len(set(species_list))

        # combine to dict
        reads_dict[sample] = [reads_sum, reads_avg, n_OTUs, n_species]

    # read stats total
    read_sum_total = 0
    for read_sum in reads_dict.values():
        read_sum_total += read_sum[0]

    #####################################################################################
    # Plot reads

    reads = [i[0] for i in reads_dict.values()]
    otus = [i[2] for i in reads_dict.values()]
    species =  [i[3] for i in reads_dict.values()]
    max_otus = max(otus) + 20

    bstats_w, bstats_h, bstats_f = int(bstats_w), int(bstats_h), int(bstats_f)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(bstats_w, bstats_h))
    ax1.bar(samples, reads)
    ax1.grid(color="grey", alpha=0.1)
    ax1.set_title('Reads', fontsize=bstats_f+2)
    ax1.tick_params(axis='y', which='major', labelsize=bstats_f)
    ax2.bar(samples, otus)
    ax2.grid(color="grey", alpha=0.1)
    ax2.set_title('OTUs', fontsize=bstats_f+2)
    ax2.tick_params(axis='y', which='major', labelsize=bstats_f)
    ax2.set_ylim(0, max_otus)
    ax3.bar(samples, species)
    ax3.grid(color="grey", alpha=0.1)
    ax3.set_title('OTUs on species level', fontsize=bstats_f+2)
    ax3.tick_params(axis='y', which='major', labelsize=bstats_f)
    ax3.set_ylim(0, max_otus)
    plt.xticks(fontsize=5, rotation=90)

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        basic_stats_directory = Path(str(path_to_outdirs) + "/" + "Basic_stats" + "/" + TaXon_table_xlsx.stem)
        basic_stats_plot = Path(str(basic_stats_directory) + "_basic_stats.pdf")
        fig.savefig(basic_stats_plot)
        plt.close()
        closing_text = "Basis stats plots are found in: " + str(path_to_outdirs) + "/Basic_stats/"
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()


    #####################################################################################
    output_list_1 = []
    output_list_2 = []
    #df = pd.DataFrame(simple_list,columns=['col1','col2'])

    output_list_1.append([' Samples',n_samples, ''])
    output_list_1.append([' OTUs',n_OTUs_total, ''])
    output_list_1.append(['Number of taxa per taxon level', '#', ''])
    output_list_1.append([' Phyla',n_Phyla, ''])
    output_list_1.append([' Classes',n_Classes, ''])
    output_list_1.append([' Orders',n_Orders, ''])
    output_list_1.append([' Families',n_Families, ''])
    output_list_1.append([' Genera',n_Genera, ''])
    output_list_1.append([' Species',n_Species, ''])
    output_list_1.append(['Database status','#', ''])
    for status, count in status_dict.items():
        output_list_1.append([" " + status,count, ''])
    output_list_1.append(['Sequence length','(bp)', ''])
    output_list_1.append([' Min', min_len_seq, ''])
    output_list_1.append([' Avg', avg_len_seq, ''])
    output_list_1.append([' Max', max_len_seq, ''])

    for sample, reads_stats in reads_dict.items():
        output_list_2.append([sample, reads_stats[1], reads_stats[0], reads_stats[2], reads_stats[3]])
    output_list_2.append(['Total reads', '', read_sum_total, '', ''])

    df_1 = pd.DataFrame(output_list_1,columns=['Category','#', ''])
    df_2 = pd.DataFrame(output_list_2,columns=["Sample", "avg reads", "total reads", "n OTUs", "n Species"])
    df_out = pd.concat([df_1, df_2], axis=1)
    df_out = df_out.replace(np.nan, '', regex=True)

    basic_stats_directory = Path(str(path_to_outdirs) + "/" + "Basic_stats" + "/" + TaXon_table_xlsx.stem)
    basic_stats_xlsx = Path(str(basic_stats_directory) + "_basic_stats.xlsx")
    df_out.to_excel(basic_stats_xlsx, index=False)

    table_1 = [['Category','#']] + df_1.values.tolist()
    table_2 = [["Sample", "avg reads", "total reads", "n OTUs", "n Species"]] + df_2.values.tolist()
    table_layout_1 = [[sg.Text('    '.join(list(map(str, row))), size=(70,1)) for col in range(1)] for row in table_1]
    layout = [ [sg.TabGroup([[sg.Tab('General information', table_layout_1),]])], [sg.Button("Close", key="Close")]]

    window_basic_stats = sg.Window('Basic stats', layout, keep_on_top=True)

    while True:
        event, values = window_basic_stats.Read()
        if event is None or event == 'Close':
            window_basic_stats.close()
            break
