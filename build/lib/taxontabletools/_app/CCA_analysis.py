def CCA_analysis(TaXon_table_xlsx, meta_data_to_test, cca_w, cca_h, cca_s, path_to_outdirs):

    import pandas as pd
    from skbio.stats.ordination import cca
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import PySimpleGUI as sg

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    IDs_list = TaXon_table_df["IDs"].values.tolist()
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    output_pdf = Path(str(path_to_outdirs) + "/" + "CCA_plots" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + ".pdf")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "CCA_plots" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + ".xlsx")

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test != {1,0}:
        print("Please use presence absence data!")
        sg.Popup("Please use presence absence data!", title=("Error"))
        raise RuntimeError

    # check for presence absence data
    # otherwise abort and print error message
    for i in Meta_data_table_df[meta_data_to_test]:
        if type(i) != int:
            print("Please use categorial numbers (ints) as meta data!")
            sg.Popup("Please use categorial numbers (ints) as meta data!", title=("Error"))
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

    # check for empty samples and later remove them from the dataframe
    # otherwise abort and print error message
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    empty_samples_list = []
    for sample in TaXon_table_samples:
        empty_sample_test = set(TaXon_table_df[sample].values.tolist())
        if empty_sample_test == {0}:
            empty_samples_list.append(sample)
    if empty_samples_list != []:
        print("Please remove empty samples first!")
        sg.Popup("Please remove empty samples first!", title=("Error"))
        raise RuntimeError

    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):
        # transpose the pa table
        df_features = TaXon_table_df[TaXon_table_samples].transpose()
        df_features.index = Meta_data_table_samples
        df_features = df_features.rename_axis("a")

        # create a constrains table from the metadata table
        constrains_list = []
        for col in Meta_data_table_df[["Samples", meta_data_to_test]].values.tolist():
            constrains_list.append(float(col[1]))
        df_constrains = pd.DataFrame(constrains_list, Meta_data_table_samples, [meta_data_to_test])
        df_constrains.index = Meta_data_table_samples
        df_constrains = df_constrains.rename_axis("a")

    ordination_result = cca(df_features, df_constrains)
    ordination_result.sample_constraints

    #######################################################################################
    # create window to ask for CCA axis to test
    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    # collect the cca proportion explained values
    proportion_explained_list = []
    for i, cca_axis in enumerate(ordination_result.proportion_explained):
        proportion_explained_list.append("CCA" + str(i+1) + " (" + str(round(cca_axis* 100, 2)) + " %)")

    cca_axis_checkboxes = list(slices([sg.Checkbox(name, key=name, size=(15,1)) for name in proportion_explained_list], 4))

    cca_window_layout = [
                [sg.Text('Check two axis to be displayed')],
                [sg.Frame(layout = cca_axis_checkboxes, title = '')],
                [sg.Button('Plot', key='Plot'), sg.Checkbox("Invert axes", key="inverted_axis")],
                [sg.Button('Back')],
                ]

    cca_window = sg.Window('CCA axis', cca_window_layout, keep_on_top=True)

    while True:
        event, values = cca_window.read()

        if event == 'Plot':
            # collect the cca axis values
            axis_to_plot = [key for key,value in values.items() if value == True and "CCA" in key]
            # pass on only if two cca axes were checked
            if len(axis_to_plot) != 2:
                sg.Popup("Please choose exactly two CCA axes", title="Error", keep_on_top=True)
            else:
                if values["inverted_axis"] == True:
                    cat2 = axis_to_plot[0].split()[0]
                    cat1 =  axis_to_plot[1].split()[0]
                else:
                    cat1 = axis_to_plot[0].split()[0]
                    cat2 =  axis_to_plot[1].split()[0]

                df_cca = ordination_result.samples[[cat1, cat2]]
                df_cca.insert(2, "label", df_constrains[meta_data_to_test].values.tolist(), True)
                plt.figure(figsize=(int(cca_w), int(cca_h)))
                plt.grid(color='gray', alpha=0.1)
                groups = df_cca.groupby("label")
                for name, group in groups:
                    plt.plot(group[cat1], group[cat2], marker="o", linestyle="", label=name, ms=cca_s)
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel(cat1 + " (" + str(round(ordination_result.proportion_explained[cat1]* 100, 2)) + " %)")
                plt.ylabel(cat2 + " (" + str(round(ordination_result.proportion_explained[cat2]* 100, 2)) + " %)")
                plt.title(meta_data_to_test)

                plt.show(block=False)
                answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
                if answer == "Yes":
                    plt.savefig(output_pdf)
                    plt.close()
                    ordination_result.samples[[cat1, cat2]].to_excel(output_xlsx)
                    closing_text = "Taxonomic resolution plots are found in: " + str(path_to_outdirs) + "/Taxonomic_resolution_plots/"
                    print(closing_text)
                    sg.Popup(closing_text, title="Finished", keep_on_top=True)
                    break
                else:
                    plt.close()

        if event == 'Back':
            break

    cca_window.close()















#
