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
from itertools import combinations
from taxontabletools.taxontable_manipulation import strip_metadata
from statistics import pvariance


def PCoA_analysis(TaXon_table_xlsx, meta_data_to_test, taxonomic_level, width, height, pcoa_s, path_to_outdirs, template, font_size, color_discrete_sequence, pcoa_dissimilarity):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('unidentified')
    TaXon_table_df = strip_metadata(TaXon_table_df)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]

    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()
    metadata_loc = Meta_data_table_df.columns.tolist().index(meta_data_to_test)

    ## drop samples with metadata called nan (= empty)
    drop_samples = [i[0] for i in Meta_data_table_df.values.tolist() if i[metadata_loc] == "nan"]

    if drop_samples != []:
        ## filter the TaXon table
        TaXon_table_df = TaXon_table_df.drop(drop_samples, axis=1)
        TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
        ## also remove empty OTUs
        row_filter_list = []
        for row in TaXon_table_df.values.tolist():
            reads = set(row[10:])
            if reads != {0}:
                row_filter_list.append(row)
        columns = TaXon_table_df.columns.tolist()
        TaXon_table_df = pd.DataFrame(row_filter_list, columns=columns)
        Meta_data_table_df = pd.DataFrame([i for i in Meta_data_table_df.values.tolist() if i[0] not in drop_samples], columns=Meta_data_table_df.columns.tolist())
        Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

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

        ## extract the relevant data
        TaXon_table_df = TaXon_table_df[[taxonomic_level] + samples]
        ## define an aggregation function to combine multiple hit of one taxonimic level
        aggregation_functions = {}
        ## define samples functions
        for sample in samples:
            ## 'sum' will calculate the sum of p/a data
            aggregation_functions[sample] = 'sum'
        ## define taxon level function
        aggregation_functions[taxonomic_level] = 'first'
        ## create condensed dataframe
        TaXon_table_df = TaXon_table_df.groupby(TaXon_table_df[taxonomic_level]).aggregate(aggregation_functions)
        if 'unidentified' in TaXon_table_df.index:
            TaXon_table_df = TaXon_table_df.drop('unidentified')

        # ## calculate distance matrix
        # if pcoa_dissimilarity == 'jaccard':
        #     distance_matrix = []
        #     for s1 in samples:
        #         distances = []
        #         array1 = [1 if i != 0 else 0 for i in TaXon_table_df[s1].values.T.tolist()]
        #         for s2 in samples:
        #             array2 = [1 if i != 0 else 0 for i in TaXon_table_df[s2].values.T.tolist()]
        #             distances.append(distance.jaccard(array1, array2))
        #         distance_matrix.append(distances)
        #
        # elif pcoa_dissimilarity == 'braycurtis':
        #     distance_matrix = []
        #     for s1 in samples:
        #         distances = []
        #         array1 = TaXon_table_df[s1].values.T.tolist()
        #         for s2 in samples:
        #             array2 = TaXon_table_df[s2].values.T.tolist()
        #             distances.append(distance.braycurtis(array1, array2))
        #         distance_matrix.append(distances)
        #
        # ## execute pcoa
        # mds = manifold.MDS(n_components=10, max_iter=30000, dissimilarity="precomputed", n_jobs=1)
        # positions = mds.fit(distance_matrix).embedding_
        # df = pd.DataFrame(positions, index=samples)
        #
        # ## calculate proportion explained
        # vars = [pvariance(df[i].values.tolist()) for i in df.columns]
        # sm = sum(vars)
        # proportion_explained = [round(i/sm*100,2) for i in vars]
        #

        jc_dm = beta_diversity(pcoa_dissimilarity, data, samples)
        ordination_result = pcoa(jc_dm)
        metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

        anosim_results = anosim(jc_dm, metadata_list, permutations=999)
        anosim_r = round(anosim_results['test statistic'], 5)
        anosim_p = anosim_results['p-value']
        textbox = meta_data_to_test + ", " + taxon_title + "<br>Anosim " + "R = " + str(anosim_r) + " " + "p = " + str(anosim_p)

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
                    [sg.Text('Check up to four axes to be displayed')],
                    [sg.Frame(layout = pcoa_axis_checkboxes, title = '')],
                    [sg.Text('Only axes >= 1 % explained variance are shown')],
                    [sg.CB("Connect categories", default=True, key="draw_mesh")],
                    [sg.Text('')],
                    [sg.Button('Plot', key='Plot')],
                    [sg.Button('Back')],
                    ]

        pcoa_window = sg.Window('PCoA axis', pcoa_window_layout, keep_on_top=True)

        while True:
            event, values = pcoa_window.read()

            draw_mesh = values["draw_mesh"]

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

                    if draw_mesh == True:
                        combinations_list =[]
                        for metadata in df_pcoa["Metadata"]:
                            ## collect all entries for the respective metadata
                            arr = df_pcoa.loc[df_pcoa['Metadata'] == metadata][[cat1, cat2, "Metadata", "Samples"]].to_numpy()
                            ## create a df for all possible combinations using itertools combinations
                            for entry in list(combinations(arr, 2)):
                                combinations_list.append(list(entry[0]))
                                combinations_list.append(list(entry[1]))
                        ## create a dataframe to draw the plot from
                        df = pd.DataFrame(combinations_list)
                        df.columns = [cat1, cat2, "Metadata", "Samples"]

                        fig = px.scatter(df, x=cat1, y=cat2, color="Metadata", text="Samples", title=textbox, color_discrete_sequence=color_discrete_sequence)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers+lines")
                        fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True, font_size=font_size, title_font_size=font_size)
                        fig.update_xaxes(title=axis_to_plot[1])
                        fig.update_yaxes(title=axis_to_plot[0])

                    else:
                        fig = px.scatter(df_pcoa, x=cat1, y=cat2, color="Metadata", text="Samples", title=textbox, color_discrete_sequence=color_discrete_sequence)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                        fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True, font_size=font_size, title_font_size=font_size)
                        fig.update_xaxes(title=axis_to_plot[1])
                        fig.update_yaxes(title=axis_to_plot[0])

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + ".pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + ".html")
                    output_xlsx = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + ".xlsx")

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
                        fig = px.scatter_3d(df, x=cat1, y=cat2, z=cat3, color="Metadata", text="Samples", title=textbox, color_discrete_sequence=color_discrete_sequence)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers+lines", line=dict(width=0.5))
                        fig.update_layout(height=int(height), width=int(width), template=template, title=textbox, showlegend=True, font_size=font_size, title_font_size=font_size)
                        fig.update_layout(scene = dict(xaxis_title=axis_to_plot[0],yaxis_title=axis_to_plot[1],zaxis_title=axis_to_plot[2]))
                    else:
                        fig = px.scatter_3d(df_pcoa, x=cat1, y=cat2, z=cat3, color="Metadata", text="Samples", color_discrete_sequence=color_discrete_sequence)
                        fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                        fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True, title=textbox, font_size=font_size, title_font_size=font_size)
                        fig.update_layout(scene = dict(xaxis_title=axis_to_plot[0],yaxis_title=axis_to_plot[1],zaxis_title=axis_to_plot[2]))

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_3d.pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_3d.html")
                    output_xlsx = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_3d.xlsx")

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
                    fig.update_layout(template=template, font_size=font_size, title_font_size=font_size)
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
                    fig.update_layout(template=template, font_size=font_size, title_font_size=font_size)
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
                    fig.update_layout(template=template, font_size=font_size, title_font_size=font_size)
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
                    fig.update_layout(template=template, font_size=font_size, title_font_size=font_size)
                    text = "PC4 (" + str(round(ordination_result.proportion_explained["PC4"]* 100, 2)) + " %)"
                    fig.add_annotation(text=text, showarrow=False, row=4, col=4)

                    ######################
                    fig.update_xaxes(showline=True, mirror=True, linewidth=1, linecolor='black')
                    fig.update_yaxes(showline=True, mirror=True, linewidth=1, linecolor='black')
                    fig.update_traces(marker_size=int(pcoa_s), mode="markers")
                    # finish plot matrix
                    fig.update_layout(height=1000, width=1000, title_text=textbox)

                    ## define output files
                    output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_matrix.pdf")
                    output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_matrix.html")

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
