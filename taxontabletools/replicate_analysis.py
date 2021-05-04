def replicate_analysis(TaXon_table_xlsx, height, width, suffix_list, path_to_outdirs, template, theme, font_size, custom_colors, clustering_unit):

    # TaXon_table_xlsx = "/Users/tillmacher/Desktop/Projects/TTT_Projects/Projects/Robo_Test_2/TaXon_tables/Robo_additional_data_samples.xlsx"
    # height = 1000
    # width = 1000
    # suffix_list = ["A", "B"]
    # path_to_outdirs = "/Users/tillmacher/Desktop/Projects/TTT_Projects/Projects/Robo_Test_2"
    # template = "simple_white"
    # theme = ["Blue", "Black", 1]
    # font_size = 16
    # import plotly.express as px
    # color_discrete_sequence = px.colors.qualitative.Plotly
    # add_annotations = True

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from statistics import mean
    from pathlib import Path
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    from matplotlib_venn import venn3
    from matplotlib.pyplot import plot, ion, show
    import matplotlib.gridspec as gridspec
    import math, os, webbrowser
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    from collections import OrderedDict

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    height = int(height)
    width = int(width)

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)

    sample_names = TaXon_table_df.columns[10:].tolist()
    OTUs = TaXon_table_df["ID"].values.tolist()

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
    fig_main_dict = {}
    reads_dict = {}

    ## create an output folder
    replicate_analysis_name = Path(TaXon_table_xlsx).name.replace(".xlsx", "")
    dirName = Path(str(path_to_outdirs) + "/Replicate_analysis/" + replicate_analysis_name)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    for sample in unique_sample_names_set:
        for i, suffix in enumerate(suffix_list):
            replicates_dict["rep_" + str(i)] = sample + "_" + suffix_list[i]

        replicate_names_list = list(replicates_dict.values())

        try:
            ## calculate the number of shared OTUs
            shared_OTUs_list = [row for row in TaXon_table_df[replicate_names_list].values.tolist() if 0 not in row]
            present_OTUs_list = [row for row in TaXon_table_df[replicate_names_list].values.tolist() if row != [0]*len(replicate_names_list)]
            perc_shared = round(len(shared_OTUs_list) / len(present_OTUs_list) * 100, 2)
            replicate_perc_shared_dict[sample] = perc_shared

            ## calculate the percentage of reads that is discarded and kept
            reads_total = sum([sum(row) for row in TaXon_table_df[replicate_names_list].values.tolist()])
            reads_kept_perc = round(sum([sum(row) for row in shared_OTUs_list]) / reads_total * 100, 2)
            reads_discarded_perc = round(100 - reads_kept_perc, 2)
            reads_dict[sample] = [reads_kept_perc, reads_discarded_perc]

            ## create left sided OTU plot
            fig_dict = {}
            for i, OTU in enumerate(present_OTUs_list):
                if 0 not in OTU:
                    fig_main_dict[i+1, sample, "Blue", "shared"] = [sum(OTU) / reads_total * 100]
                else:
                    fig_main_dict[i+1, sample, "Red", "non-shared"] = [sum(OTU) / reads_total * 100]
        except:
            print("Warning! No replicates found for:    " + sample)

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

    ########################################################################################################################
    ## figure 1 shared OTUs
    samples = list(replicate_perc_shared_dict.keys())
    shared_otus = list(replicate_perc_shared_dict.values())
    y_title = "shared " + clustering_unit
    fig = px.bar(x=samples, y=shared_otus, labels={"y": y_title, "x": "Sample", "text": y_title}, text=shared_otus)
    y_title = 'shared ' + clustering_unit + ' (%)'
    fig.update_yaxes(title=y_title, range=[0, 100], dtick=10, autorange=False)
    fig.update_xaxes(title='', tickmode='linear')
    fig.update_xaxes(tickangle=-90)
    fig.update_layout(width=int(width), height=int(height), template=template, font_size=font_size, title_font_size=font_size)
    fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1.5, opacity=opacity_value)

    ## write files
    output_pdf = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared_" + clustering_unit + ".pdf")
    output_html = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared_" + clustering_unit + ".html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ########################################################################################################################
    ## figure 2 kept/discarded reads
    samples = list(reads_dict.keys())
    #discarded_reads = [reads[1] for reads in list(reads_dict.values())]
    shared_reads = [reads[0] for reads in list(reads_dict.values())]
    fig = px.bar(x=samples, y=shared_reads, labels={"y": "shared reads (%)", "x": "Sample", "text": "shared reads (%)"}, text=shared_reads)
    fig.update_yaxes(title='shared reads (%)', range=[0, 100], dtick=10, autorange=False)
    fig.update_xaxes(title='', tickmode='linear')
    fig.update_xaxes(tickangle=-90)
    fig.update_layout(width=int(width), height=int(height), template=template, font_size=font_size, title_font_size=font_size)
    fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1.5, opacity=opacity_value)

    ## write files
    output_pdf2 = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared" + clustering_unit + "_reads.pdf")
    output_html2 = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared" + clustering_unit + "_reads.html")
    fig.write_image(str(output_pdf2))
    fig.write_html(str(output_html2))

    ########################################################################################################################
    ## figure 3 OTU left side plot
    ## sort the dict by abundance
    fig_main_dict_sorted = dict(sorted(fig_main_dict.items(), key=lambda item: item[1], reverse=True))
    ## collect y values >> read abundances
    ## collect x values >> rank
    # y1, x1, m1, t1 = [], [], [], []
    # y2, x2, m2, t2 = [], [], [], []
    # i = 0
    #
    # for key, value in fig_main_dict_sorted.items():
    #     if key[3] == 'shared':
    #         y1 = y1 + value
    #         x1.append(i)
    #         m1.append(key[2])
    #         t1.append(key[3])
    #         i += 1
    #     else:
    #         y2 = y2 + value
    #         x2.append(i)
    #         m2.append(key[2])
    #         t2.append(key[3])
    #         i += 1
    #
    # max_reads = math.ceil(max(y1 + y2)) +1
    # n_ranked_OTUs = i
    # n_shared = len(y1)
    # n_nonshared = len(y2)
    # name1 = "shared (n=" + str(n_shared) + ")"
    # name2 = "non shared (n=" + str(n_nonshared) + ")"
    # c1 = color_discrete_sequence[0]
    # c2 = color_discrete_sequence[1]
    #
    # fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    # fig.add_trace(go.Scatter(x=x1, y=y1, mode='markers', marker=dict(size=8, color=c1), name=name1),row=1, col=1)
    # fig.add_trace(go.Scatter(x=x2, y=y2, mode='markers', marker=dict(size=8, color=c2), name=name2),row=2, col=1)
    # fig.update_xaxes(title='', showticklabels=True, row=1, col=1)
    # fig.update_xaxes(title='ranked OTUs (by read abundance)', showticklabels=True, row=2, col=1)
    # if log_transform == True:
    #     fig.update_yaxes(title="reads (log)", range=[0,max_reads])
    #     out_sub = "_log_"
    # else:
    #     fig.update_yaxes(title="reads (%)", range=[0,105])
    #     out_sub = "_rel_"
    # fig.update_layout(width=int(width), height=int(height), template=template, font_size=font_size, title_font_size=font_size)
    #
    # ## add annotations to both plots
    # if add_annotations == True:
    #     for annotation in [1.0, 0.1, 0.01]:
    #         ## store annotation as text
    #         text = "<" + str(annotation) + "%"
    #         ## search for the x-axis rank of the annotation
    #         x_pos = []
    #         try:
    #             for rank, reads in zip(x1, y1):
    #                 if round(reads, 3) <= annotation:
    #                     x_pos = rank
    #                     break
    #         except:
    #             pass
    #         if x_pos != []:
    #             ## y_pos is annotation
    #             fig.add_annotation(x=x_pos, y=2.5,
    #                             text=text, showarrow=True, font=dict(size=font_size-2),
    #                             align="center", arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="Black", ax=20, ay=-30,
    #                             bordercolor=c1, borderwidth=1, borderpad=3, bgcolor=c1, opacity=0.9,
    #                             row=1, col=1)
    #
    #         ## search for the x-axis rank of the annotation
    #         x_pos = []
    #         try:
    #             for rank, reads in zip(x2, y2):
    #                 if round(reads, 3) <= annotation:
    #                     x_pos = rank
    #                     break
    #         except:
    #             pass
    #         if x_pos != []:
    #             ## y_pos is annotation
    #             fig.add_annotation(x=x_pos, y=2.5,
    #                             text=text, showarrow=True, font=dict(size=font_size-2),
    #                             align="center", arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="Black", ax=20, ay=-30,
    #                             bordercolor=c2, borderwidth=1, borderpad=3, bgcolor=c2, opacity=0.9,
    #                             row=2, col=1)

    ########################################################################################################################
    ## figure 4 OTU bar plot v2
    y1, y2 = [], []

    for key, value in fig_main_dict_sorted.items():
        if key[3] == 'shared':
            y1 = y1 + value
        else:
            y2 = y2 + value

    categories = [[100, 10], [10,1], [1,0.1], [0.1,0]]
    bar_plot_dict = {}
    n_OTUs_shared, n_OTUs_nonshared, names = [], [], []

    fig = go.Figure()

    for category in categories:
        upper = category[0]
        lower = category[1]
        shared = len([y for y in y1 if (y > lower and y < upper)])
        nonshared = len([y for y in y2 if (y > lower and y < upper)])
        n_OTUs = shared + nonshared
        shared_perc = shared / n_OTUs * 100
        nonshared_perc = nonshared / n_OTUs * 100
        n_OTUs_shared.append(shared_perc)
        n_OTUs_nonshared.append(nonshared_perc)
        if category != [0.1,0]:
            text = str(category[0]) + "%-" + str(category[1]) + "%"
            names.append(text)
            fig.add_annotation(x=text, y=100, text="n="+str(n_OTUs), font=dict(size=font_size-2), showarrow=False, yshift=10)
        else:
            text = "<0.1%"
            names.append(text)
            fig.add_annotation(x=text, y=100, text="n="+str(n_OTUs), font=dict(size=font_size-2), showarrow=False, yshift=10)


    fig.add_trace(go.Bar(x=names, y=n_OTUs_shared, name='shared', marker_color=custom_colors[0]))
    fig.add_trace(go.Bar(x=names, y=n_OTUs_nonshared, name='non-shared', marker_color=custom_colors[1]))

    fig.update_layout(width=int(width), height=int(height), template=template, font_size=font_size, title_font_size=font_size)
    y_title = 'rel. ' + clustering_unit + ' per bin'
    fig.update_yaxes(title=y_title)
    fig.update_xaxes(title='read abundance')

    ## write files
    output_pdf3 = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared_nonshared.pdf")
    output_html3 = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_shared_nonshared.html")
    fig.write_image(str(output_pdf3))
    fig.write_html(str(output_html3))

    ########################################################################################################################

    ## write statistics file
    output_txt = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_stats.txt")
    f = open(output_txt, "w")
    avg_shared_otus = round(mean(shared_otus),2 )
    avg_shared_reads = round(mean(shared_reads),2 )
    n_samples = len(samples)
    text = "Average shared " + clustering_unit + ": " + str(avg_shared_otus) + "%\n" + "Average shared reads: " + str(avg_shared_reads) + "%\n" + "Number of samples: " + str(n_samples)
    f.write(text)
    f.close()

    ## ask to show file
    answer = sg.PopupYesNo(text + '\n\nShow all three plots?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html3))
        webbrowser.open('file://' + str(output_html2))
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "The three plots are found under:\n" + "Projects/Replicate_analysis/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write log
    from taxontabletools.create_log import ttt_log
    ttt_log("replicate analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, "nan", path_to_outdirs)
