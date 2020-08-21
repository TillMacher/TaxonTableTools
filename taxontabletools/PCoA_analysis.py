def PCoA_analysis(TaXon_table_xlsx, meta_data_to_test, pcoa_w, pcoa_h, pcoa_s, path_to_outdirs):

    import pandas as pd
    import numpy as np
    from skbio.diversity import beta_diversity
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import anosim
    import matplotlib.pyplot as plt
    from pathlib import Path
    import PySimpleGUI as sg

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    output_pdf = Path(str(path_to_outdirs) + "/" + "PCoA_plots" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + ".pdf")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "PCoA_plots" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + ".xlsx")

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test != {1,0}:
        print("Please use presence absence data!")
        sg.Popup("Please use presence absence data!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == len(Meta_data_table_df['Samples'].tolist()):
        print("The meta data is unique for all samples. Please adjust the meta data table!")
        sg.Popup("The meta data is unique for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        print("The meta data is similar for all samples. Please adjust the meta data table!")
        sg.Popup("The meta data is similar for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        samples = Meta_data_table_samples
        data = TaXon_table_df[TaXon_table_samples].transpose().values.tolist()
        jc_dm = beta_diversity("jaccard", data, samples)
        ordination_result = pcoa(jc_dm)
        metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

        anosim_results = anosim(jc_dm, metadata_list, permutations=999)
        anosim_r = round(anosim_results['test statistic'], 5)
        anosim_p = anosim_results['p-value']
        textbox = meta_data_to_test + "\nAnosim " + "R = " + str(anosim_r) + " " + "p = " + str(anosim_p)

        #######################################################################################
        # create window to ask for PCoA axis to test
        def slices(list, slice):
            for i in range(0, len(list), slice):
                yield list[i : i + slice]

        # collect the PCoA proportion explained values
        proportion_explained_list = []
        for i, pcoa_axis in enumerate(ordination_result.proportion_explained):
            proportion_explained_list.append("PC" + str(i+1) + " (" + str(round(pcoa_axis* 100, 2)) + " %)")

        pcoa_axis_checkboxes = list(slices([sg.Checkbox(name, key=name, size=(15,1)) for name in proportion_explained_list], 4))

        pcoa_window_layout = [
                    [sg.Text('Check two axis to be displayed')],
                    [sg.Frame(layout = pcoa_axis_checkboxes, title = '')],
                    [sg.Button('Plot', key='Plot'), sg.Checkbox("Invert axes", key="inverted_axis")],
                    [sg.Button('Back')],
                    ]

        pcoa_window = sg.Window('PCoA axis', pcoa_window_layout, keep_on_top=True)

        while True:
            event, values = pcoa_window.read()

            if event == 'Plot':
                # collect the pcoa axis values
                axis_to_plot = [key for key,value in values.items() if value == True and "PC" in key]
                # pass on only if two pcoa axes were checked
                if len(axis_to_plot) != 2:
                    sg.Popup("Please choose exactly two PCoA axes", title="Error", keep_on_top=True)
                else:
                    if values["inverted_axis"] == True:
                        cat1 = axis_to_plot[0].split()[0]
                        cat2 =  axis_to_plot[1].split()[0]
                    else:
                        cat2 = axis_to_plot[0].split()[0]
                        cat1 =  axis_to_plot[1].split()[0]
                    df_pcoa = ordination_result.samples[[cat1, cat2]]
                    df_pcoa.insert(2, "label", Meta_data_table_df[meta_data_to_test].values.tolist(), True)
                    plt.figure(figsize=(int(pcoa_w), int(pcoa_h)))
                    plt.grid(color='gray', alpha=0.1)
                    groups = df_pcoa.groupby("label")
                    for name, group in groups:
                        plt.plot(group[cat1], group[cat2], marker="o", linestyle="", label=name, ms=pcoa_s)
                    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    plt.xlabel(cat1 + " (" + str(round(ordination_result.proportion_explained[cat1]* 100, 2)) + " %)")
                    plt.ylabel(cat2 + " (" + str(round(ordination_result.proportion_explained[cat2]* 100, 2)) + " %)")
                    plt.title(textbox)

                    plt.show(block=False)
                    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
                    if answer == "Yes":
                        plt.savefig(output_pdf)
                        plt.close()
                        ordination_result.samples[[cat1, cat2]].to_excel(output_xlsx)
                        closing_text = "\n" + "PCoA plots are found in: " + str(path_to_outdirs) + "/PCoA_plots/"
                        print(closing_text)
                        sg.Popup(closing_text, title="Finished", keep_on_top=True)

                        from taxontabletools.create_log import ttt_log
                        ttt_log("pcoa analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)

                        break
                    else:
                        plt.close()

            if event == 'Back':
                break

        pcoa_window.close()
