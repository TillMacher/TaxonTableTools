def PCoA_analysis(TaXon_table_xlsx, meta_data_to_test, width, height, pcoa_s, draw_mesh, path_to_outdirs, template):
    import pandas as pd
    import numpy as np
    from skbio.diversity import beta_diversity
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import anosim
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    from pathlib import Path
    import PySimpleGUI as sg
    import os, webbrowser

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test != {1,0}:
        sg.Popup("Please use presence absence data!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == len(Meta_data_table_df['Samples'].tolist()):
        sg.Popup("The meta data is unique for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
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
            if round(pcoa_axis* 100, 2) >= 1:
                proportion_explained_list.append("PC" + str(i+1) + " (" + str(round(pcoa_axis* 100, 2)) + " %)")

        pcoa_axis_checkboxes = list(slices([sg.Checkbox(name, key=name, size=(15,1)) for name in proportion_explained_list], 10))

        pcoa_window_layout = [
                    [sg.Text('Check two axis to be displayed')],
                    [sg.Frame(layout = pcoa_axis_checkboxes, title = '')],
                    [sg.Text('Only axes >= 1 % explained variance are shown')],
                    [sg.Text('')],
                    [sg.Button('Plot', key='Plot'), sg.Text(''), sg.Button('Plot matrix')],
                    [sg.Button('Back')],
                    ]

        pcoa_window = sg.Window('PCoA axis', pcoa_window_layout, keep_on_top=True)

        while True:
            event, values = pcoa_window.read()

            if event is None or event == 'Back':
                break

            if event == 'Plot':

                ## create a subfolder for better sorting and overview
                dirName = Path(str(path_to_outdirs) + "/" + "PCoA_plots" + "/" + TaXon_table_xlsx.stem + "/")
                if not os.path.exists(dirName):
                    os.mkdir(dirName)

                # collect the pcoa axis values
                axis_to_plot = [key for key,value in values.items() if value == True and "PC" in key]
                # pass on only if two pcoa axes were checked
                if len(axis_to_plot) == 2:
                    cat1 = axis_to_plot[1].split()[0]
                    cat2 = axis_to_plot[0].split()[0]

                    df_pcoa = ordination_result.samples[[cat1, cat2]]
                    df_pcoa.insert(2, "Metadata", Meta_data_table_df[meta_data_to_test].values.tolist(), True)
                    df_pcoa.insert(3, "Samples", Meta_data_table_df["Samples"].values.tolist(), True)

                    fig = px.scatter(df_pcoa, x=cat1, y=cat2, color="Metadata", text="Samples", title=textbox)
                    fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                    fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True)
                    fig.update_xaxes(title=axis_to_plot[1])
                    fig.update_yaxes(title=axis_to_plot[0])

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + ".pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + ".html")
                    output_xlsx = Path(str(dirName) + "/" + meta_data_to_test + ".xlsx")

                    ## write files
                    fig.write_image(str(output_pdf))
                    fig.write_html(str(output_html))
                    ordination_result.samples[[cat1, cat2]].to_excel(output_xlsx)

                    ## ask to show file
                    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
                    if answer == "Yes":
                        webbrowser.open('file://' + str(output_html))

                    ## print closing text
                    closing_text = "\n" + "PCoA plots are found in: " + str(path_to_outdirs) + "/PCoA_plots/"
                    sg.Popup(closing_text, title="Finished", keep_on_top=True)

                    ## write to log
                    from taxontabletools.create_log import ttt_log
                    ttt_log("pcoa analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)
                    break

                elif len(axis_to_plot) == 3:
                    cat1 = axis_to_plot[0].split()[0]
                    cat2 = axis_to_plot[1].split()[0]
                    cat3 = axis_to_plot[2].split()[0]

                    df_pcoa = ordination_result.samples[[cat1, cat2, cat3]]
                    df_pcoa.insert(3, "Metadata", Meta_data_table_df[meta_data_to_test].values.tolist(), True)
                    df_pcoa.insert(4, "Samples", Meta_data_table_df["Samples"].values.tolist(), True)

                    from itertools import combinations

                    ## check if lines are to be drawn between the dots
                    if draw_mesh == True:
                        combinations_list =[]
                        for metadata in df_pcoa["Metadata"]:
                            ## collect all entries for the respective metadata
                            arr = df_pcoa.loc[df_pcoa['Metadata'] == metadata][[cat1, cat2, cat3, "Metadata", "Samples"]].to_numpy()
                            ## create a df for all possible combinations using itertools combinations
                            for entry in list(combinations(arr, 2)):
                                combinations_list.append(list(entry[0]))
                                combinations_list.append(list(entry[1]))
                        ## create a dataframe to draw the plot from
                        df = pd.DataFrame(combinations_list)
                        df.columns = [cat1, cat2, cat3, "Metadata", "Samples"]
                        ## draw the plot
                        fig = px.scatter_3d(df, x=cat1, y=cat2, z=cat3, color="Metadata", text="Samples", title=textbox)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers+lines", line=dict(width=0.5))
                        fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True)
                        fig.update_layout(scene = dict(xaxis_title=axis_to_plot[0],yaxis_title=axis_to_plot[1],zaxis_title=axis_to_plot[2]))
                    else:
                        fig = px.scatter_3d(df_pcoa, x=cat1, y=cat2, z=cat3, color="Metadata", text="Samples")
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                        fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True)
                        fig.update_layout(scene = dict(xaxis_title=axis_to_plot[0],yaxis_title=axis_to_plot[1],zaxis_title=axis_to_plot[2]))

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_3d.pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + "_3d.html")
                    output_xlsx = Path(str(dirName) + "/" + meta_data_to_test + "_3d.xlsx")

                    ## write output files
                    fig.write_image(str(output_pdf))
                    fig.write_html(str(output_html))
                    ordination_result.samples[[cat1, cat2]].to_excel(output_xlsx)

                    ## ask to show file
                    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
                    if answer == "Yes":
                        webbrowser.open('file://' + str(output_html))

                    ## print closing text
                    closing_text = "PCoA plots are found in: " + str(path_to_outdirs) + "/PCoA_plots/"
                    sg.Popup(closing_text, title="Finished", keep_on_top=True)

                    ## write log file
                    from taxontabletools.create_log import ttt_log
                    ttt_log("pcoa analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)
                    break

                else:
                    sg.Popup("Please choose not more than 3 PCoA axes", title="Error", keep_on_top=True)

            if event == 'Plot matrix':
                if len(proportion_explained_list) >= 4:

                    ## create a subfolder for better sorting and overview
                    dirName = Path(str(path_to_outdirs) + "/" + "PCoA_plots" + "/" + TaXon_table_xlsx.stem + "/")
                    if not os.path.exists(dirName):
                        os.mkdir(dirName)

                    df_pcoa = ordination_result.samples[["PC1", "PC2", "PC3", "PC4"]]
                    df_pcoa.insert(4, "Metadata", Meta_data_table_df[meta_data_to_test].values.tolist(), True)
                    df_pcoa.insert(5, "Sample", Meta_data_table_df["Samples"].values.tolist(), True)

                    fig = make_subplots(rows=4, cols=4)
                    ########### 1 ###########
                    fig.add_trace(go.Scatter(),row=1, col=1)
                    fig.update_layout(template=template)
                    text = "PC1 (" + str(round(ordination_result.proportion_explained["PC1"]* 100, 2)) + " %)"
                    fig.add_annotation(text=text, showarrow=False)
                    fig.update_xaxes(showticklabels=False, showgrid=False)
                    fig.update_yaxes(showticklabels=False, showgrid=False)
                    ########### 2 ###########
                    df = df_pcoa[["PC1", "PC2", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        #fig = px.scatter(df_pcoa, x="PC1", y="PC2", , )
                        fig.add_trace(go.Scatter(   x=df_metadata["PC1"].values.tolist(),
                                                    y=df_metadata["PC2"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    text=df_metadata["Sample"].values.tolist()),row=1, col=2)
                    ########### 3 ###########
                    df = df_pcoa[["PC1", "PC3", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        #fig = px.scatter(df_pcoa, x="PC1", y="PC2", , )
                        fig.add_trace(go.Scatter(   x=df_metadata["PC1"].values.tolist(),
                                                    y=df_metadata["PC3"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    showlegend=False,
                                                    text=df_metadata["Sample"].values.tolist()),row=1, col=3)
                    ########### 4 ###########
                    df = df_pcoa[["PC1", "PC4", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        fig.add_trace(go.Scatter(   x=df_metadata["PC1"].values.tolist(),
                                                    y=df_metadata["PC4"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    showlegend=False,
                                                    text=df_metadata["Sample"].values.tolist()),row=1, col=4)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                        fig.update_xaxes(showgrid=False,row=1, col=4)
                        fig.update_yaxes(showgrid=False,row=1, col=4)
                    ########### 5 ###########
                    fig.add_trace(go.Scatter(),row=2, col=2)
                    fig.update_layout(template=template)
                    text = "PC2 (" + str(round(ordination_result.proportion_explained["PC2"]* 100, 2)) + " %)"
                    fig.add_annotation(text=text, showarrow=False, row=2, col=2)
                    ########### 6 ###########
                    df = df_pcoa[["PC2", "PC3", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        #fig = px.scatter(df_pcoa, x="PC1", y="PC2", , )
                        fig.add_trace(go.Scatter(   x=df_metadata["PC2"].values.tolist(),
                                                    y=df_metadata["PC3"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    showlegend=False,
                                                    text=df_metadata["Sample"].values.tolist()),row=2, col=3)
                    ########### 7 ###########
                    df = df_pcoa[["PC2", "PC4", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        fig.add_trace(go.Scatter(   x=df_metadata["PC2"].values.tolist(),
                                                    y=df_metadata["PC4"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    showlegend=False,
                                                    text=df_metadata["Sample"].values.tolist()),row=2, col=4)
                    ########### 8 ###########
                    fig.add_trace(go.Scatter(),row=3, col=3)
                    fig.update_layout(template=template)
                    text = "PC3 (" + str(round(ordination_result.proportion_explained["PC3"]* 100, 2)) + " %)"
                    fig.add_annotation(text=text, showarrow=False, row=3, col=3)
                    ########### 9 ###########
                    df = df_pcoa[["PC3", "PC4", "Metadata", "Sample"]]
                    for metadata in set(metadata_list):
                        df_metadata = df[df['Metadata'] == metadata]
                        #fig = px.scatter(df_pcoa, x="PC1", y="PC2", , )
                        fig.add_trace(go.Scatter(   x=df_metadata["PC3"].values.tolist(),
                                                    y=df_metadata["PC4"].values.tolist(),
                                                    mode='markers',
                                                    name=metadata,
                                                    showlegend=False,
                                                    text=df_metadata["Sample"].values.tolist()),row=3, col=4)
                    ########### 5 ###########
                    fig.add_trace(go.Scatter(),row=4, col=4)
                    fig.update_layout(template=template)
                    text = "PC4 (" + str(round(ordination_result.proportion_explained["PC4"]* 100, 2)) + " %)"
                    fig.add_annotation(text=text, showarrow=False, row=4, col=4)

                    ######################
                    fig.update_xaxes(showline=True, mirror=True, linewidth=1, linecolor='black')
                    fig.update_yaxes(showline=True, mirror=True, linewidth=1, linecolor='black')
                    fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                    # finish plot matrix
                    fig.update_layout(height=1000, width=1000, title_text=textbox)

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_matrix.pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + "_matrix.html")

                    ## write output files
                    fig.write_image(str(output_pdf))
                    fig.write_html(str(output_html))

                    ## ask to show file
                    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
                    if answer == "Yes":
                        webbrowser.open('file://' + str(output_html))

                    ## print closing text
                    closing_text = "\n" + "PCoA plots are found in: " + str(path_to_outdirs) + "/PCoA_plots/"
                    sg.Popup(closing_text, title="Finished", keep_on_top=True)

                    ## write to log file
                    from taxontabletools.create_log import ttt_log
                    ttt_log("pcoa analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)
                    break
                else:
                    sg.Popup("There must be at least 4 PCoA axis available to plot the matrix!")

        pcoa_window.close()

    else:
        sg.PopupError("The sample of both the TaXon table and the metadata table have to match!")
