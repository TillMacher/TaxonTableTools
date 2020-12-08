def replicate_analysis(TaXon_table_xlsx, height, width, suffix_list, path_to_outdirs, template, theme):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    from matplotlib_venn import venn3
    from matplotlib.pyplot import plot, ion, show
    import matplotlib.gridspec as gridspec
    import math, os
    import plotly.express as px

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    height = int(height)
    width = int(width)

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)

    sample_names = TaXon_table_df.columns[10:].tolist()
    OTUs = TaXon_table_df["IDs"].values.tolist()

    derep_sample_names_dict =  {}
    unique_sample_names_list = []
    replicates_dict = {}

    for sample in sample_names:
        sample_name = sample.split("_")[0:-1]
        unique_sample_names_list.append("_".join(sample_name))

    unique_sample_names_set = sorted(set(unique_sample_names_list))

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(unique_sample_names_set) + 1
    ############################################################################

    replicate_perc_shared_dict = {}

    replicate_analysis_name = Path(TaXon_table_xlsx).name.replace(".xlsx", "")
    dirName = Path(str(path_to_outdirs) + "/Replicate_analysis/" + replicate_analysis_name)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    for sample in unique_sample_names_set:

        for i, suffix in enumerate(suffix_list):
            replicates_dict["rep_" + str(i)] = sample + "_" + suffix_list[i]

        replicate_names_list = list(replicates_dict.values())

        try:
            # first collect the present OTUs from all replicates
            replicate_comparison_dict = {}
            for replicate in replicate_names_list:
                replicate_comparison_list = []
                for OTU in TaXon_table_df[["IDs", replicate]].values.tolist():
                    if OTU[1] != 0:
                        replicate_comparison_list.append(OTU[0])
                replicate_comparison_dict[replicate] = replicate_comparison_list

            # now compare them
            if len(replicate_comparison_dict.keys()) == 2:
                replicate_a = list(replicate_comparison_dict.keys())[0]
                replicate_b = list(replicate_comparison_dict.keys())[1]

                taxa_labels_a = sorted(replicate_comparison_dict[replicate_a])
                taxa_labels_b = sorted(replicate_comparison_dict[replicate_b])

                a_only = set(taxa_labels_a) - set(taxa_labels_b)
                len_a_only = len(a_only)
                b_only = set(taxa_labels_b) - set(taxa_labels_a)
                len_b_only = len(b_only)
                shared = set(taxa_labels_a) & set(taxa_labels_b)
                len_shared = len(shared)
                perc_shared = round(len_shared / (len_a_only + len_b_only + len_shared) * 100, 2)
                replicate_perc_shared_dict[sample] = perc_shared

                plt.figure(figsize=(10, 5))
                venn2(subsets = (len_a_only, len_b_only, len_shared), set_labels = (replicate_a, replicate_b))
                output_pdf = Path(str(dirName) + "/" + sample + ".pdf")
                plt.savefig(output_pdf, bbox_inches='tight')

            elif len(replicate_comparison_dict.keys()) == 3:

                replicate_a = list(replicate_comparison_dict.keys())[0]
                replicate_b = list(replicate_comparison_dict.keys())[1]
                replicate_c = list(replicate_comparison_dict.keys())[2]

                taxa_labels_a = sorted(replicate_comparison_dict[replicate_a])
                taxa_labels_b = sorted(replicate_comparison_dict[replicate_b])
                taxa_labels_c = sorted(replicate_comparison_dict[replicate_c])

                a_only = set(taxa_labels_a) - set(taxa_labels_b) - set(taxa_labels_c)
                len_a_only = len(a_only)
                b_only = set(taxa_labels_b) - set(taxa_labels_a) - set(taxa_labels_c)
                len_b_only = len(b_only)
                c_only = set(taxa_labels_c) - set(taxa_labels_a) - set(taxa_labels_b)
                len_c_only = len(c_only)

                shared_all = set(taxa_labels_a) & set(taxa_labels_b) & set(taxa_labels_c)
                len_shared_all = len(shared_all)
                shared_a_b = set(taxa_labels_a) & set(taxa_labels_b) - set(taxa_labels_c)
                len_shared_a_b = len(shared_a_b)
                shared_a_c = set(taxa_labels_a) & set(taxa_labels_c) - set(taxa_labels_b)
                len_shared_a_c = len(shared_a_c)
                shared_b_c = set(taxa_labels_b) & set(taxa_labels_c) - set(taxa_labels_a)
                len_shared_b_c = len(shared_b_c)

                perc_shared = round(len_shared / (len_a_only + len_b_only + len_c_only + len_shared_all + len_shared_a_b + len_shared_a_c + len_shared_b_c) * 100, 2)
                replicate_perc_shared_dict[sample] = perc_shared

                plt.figure(figsize=(10, 5))
                venn3(subsets = (len_a_only, len_b_only, len_shared_a_b, len_c_only, len_shared_a_c, len_shared_b_c, len_shared_all), set_labels = (replicate_a, replicate_b, replicate_c))
                output_pdf = Path(str(dirName) + "/" + sample + ".pdf")
                plt.savefig(output_pdf, bbox_inches='tight')
                plt.close()

            else:
                sg.Popup("Warning! More than 3 replicates are not supported for Venn diagrams:    " + sample)

            plt.close('all')

        except:
            sg.Popup("Warning! No replicates found for:    " + sample)

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    samples = list(replicate_perc_shared_dict.keys())
    shared_otus = list(replicate_perc_shared_dict.values())

    fig = px.bar(x=samples, y=shared_otus, labels={"y": "Shared OTUs", "x": "Sample", "text": "Shared OTUs"}, text=shared_otus)
    fig.update_yaxes(title='Shared OTUs (%)', range=[0, 100], dtick=10, autorange=False)
    fig.update_xaxes(title='')
    fig.update_layout(width=int(width), height=int(height), template=template)
    fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1.5, opacity=opacity_value)

    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        fig.show()
    bar_pdf = Path(str(dirName) + "/" + TaXon_table_xlsx.name.replace(".xlsx", "_replicate_analysis.pdf"))
    bar_html = Path(str(dirName) + "/" + TaXon_table_xlsx.name.replace(".xlsx", "_replicate_analysis.html"))
    fig.write_image(str(bar_pdf))
    fig.write_html(str(bar_html))
    closing_text = "Plots are found under:\n" + "Projects/Replicate_analysis/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
    from taxontabletools.create_log import ttt_log
    ttt_log("replicate analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, "nan", path_to_outdirs)
