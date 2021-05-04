def NMDS_analysis(TaXon_table_xlsx, meta_data_to_test, taxonomic_level, width, height, nmds_s, max_iter_val, n_init_val, path_to_outdirs, template, font_size, color_discrete_sequence, nmds_dissimilarity):

    import pandas as pd
    import numpy as np
    from skbio.diversity import beta_diversity
    from sklearn.manifold import MDS
    import plotly.graph_objects as go
    import plotly.express as px
    from pathlib import Path
    import PySimpleGUI as sg
    import os, webbrowser
    from itertools import combinations

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    ## create a subfolder for better sorting and overview
    dirName = Path(str(path_to_outdirs) + "/" + "NMDS_plots" + "/" + TaXon_table_xlsx.stem + "/")
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == len(Meta_data_table_df['Samples'].tolist()):
        sg.Popup("The meta data is unique for all samples. Please adjust the meta data table!", title=("Error"))
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
        df_new = TaXon_table_df.groupby(TaXon_table_df[taxonomic_level]).aggregate(aggregation_functions)
        if 'unidentified' in df_new.index:
            df_new = df_new.drop('unidentified')

        ## collect reads
        data = df_new[samples].transpose().values.tolist()
        ## calculate jaccard distances
        jaccard_dm = beta_diversity(nmds_dissimilarity, data, samples)

        ## NMDS function
        def nmds_function(matrix, dimensions):
            nmds = MDS(n_components=dimensions, metric=False, dissimilarity='precomputed', max_iter=int(max_iter_val), n_init=int(n_init_val))
            nmds_results = nmds.fit(jaccard_dm[:100])
            stress = round(nmds_results.stress_, 2)
            nmds_array = nmds_results.embedding_
            return({"stress":stress,"nmds_results":nmds_array})


        answer = sg.PopupOKCancel("The NMDS calculation may take a while. Continue?")

        if answer == "OK":
            ## test different dimensions
            nmds_results_dict = {}
            stress_dict = {}
            for i in range(1,11):
                nmds_results = nmds_function(jaccard_dm, i)
                nmds_results_dict[i] = nmds_results
                stress_dict[i] = nmds_results["stress"]

            ####################################################################################################
            win2_active = True

            layout2 = [[sg.Text("NMDS analysis options", size=(20,1))],
                    [sg.CB("Show stress plot", default=True, key="stress_plot")],
                    [sg.CB("Show NMDS 2D plot", default=True, key="2d_plot")],
                    [sg.CB("Show NMDS 3D plot", default=True, key="3d_plot")],
                    [sg.CB("Connect categories", default=True, key="draw_mesh")],
                    [sg.Text("")],
                    [sg.Button("Apply")]]

            win2 = sg.Window('NMDS analysis', layout2, keep_on_top=False)

            while True:
                event2, values2 = win2.Read()
                if event2 is None or event2 == 'Apply':
                    win2.close()
                    win2_active = False
                    break

            ####################################################################################################

            ## plot stress and dimensions
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=list(stress_dict.keys()), y=list(stress_dict.values()), mode='markers+lines', name=sample, marker=dict(color="Blue", size=int(10))))
            fig.update_layout(showlegend=False, xaxis_title="Dimensions", yaxis_title="Stress")
            fig.update_layout(height=int(600), width=int(800), template=template, showlegend=False, font_size=font_size, title_font_size=font_size)

            ## define output files
            output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_stress.pdf")
            output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_stress.html")
            ## write output files
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))
            ## ask to show file
            if values2['stress_plot'] == True:
                webbrowser.open('file://' + str(output_html))

            ####################################################################################################

            ## plot 2D
            stress = stress_dict[2]
            if values2["draw_mesh"] == True:
                ## create dataframe from NMDS results
                nmds_results_df = pd.DataFrame(nmds_results_dict[2]["nmds_results"], index=[samples])
                nmds_results_df.rename(columns={0:'NMDS1', 1:'NMDS2'}, inplace=True)
                nmds_results_df["Sample"] = samples
                nmds_results_df[meta_data_to_test] = Meta_data_table_df[meta_data_to_test].values.tolist()

                combinations_list =[]
                for metadata in nmds_results_df[meta_data_to_test]:
                    ## collect all entries for the respective metadata
                    arr = nmds_results_df.loc[nmds_results_df[meta_data_to_test] == metadata][['NMDS1', 'NMDS2', meta_data_to_test, "Sample"]].to_numpy()
                    ## create a df for all possible combinations using itertools combinations
                    for entry in list(combinations(arr, 2)):
                        combinations_list.append(list(entry[0]))
                        combinations_list.append(list(entry[1]))

                ## create a dataframe to draw the plot from
                df = pd.DataFrame(combinations_list)
                df.columns = ['NMDS1', 'NMDS2', meta_data_to_test, "Sample"]

                ## plot NMDS
                fig = go.Figure()
                fig = px.scatter(df, x="NMDS1", y="NMDS2", hover_data=['Sample'], color=meta_data_to_test, color_discrete_sequence=color_discrete_sequence)
                fig.update_traces(marker_size=int(nmds_s), mode="markers+lines", line=dict(width=0.5))
                fig.update_layout(title="Stress=" + str(stress), yaxis_title="NMDS1", xaxis_title="NMDS2")
                fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True, font_size=font_size, title_font_size=font_size)

            else:

                ## create dataframe from NMDS results
                nmds_results_df = pd.DataFrame(nmds_results_dict[2]["nmds_results"], index=[samples])
                nmds_results_df.rename(columns={0:'X', 1:'Y'}, inplace=True)
                nmds_results_df[meta_data_to_test] = Meta_data_table_df[meta_data_to_test].values.tolist()
                nmds_results_df["Sample"] = samples
                ## plot NMDS
                fig = go.Figure()
                fig = px.scatter(nmds_results_df, x="X", y="Y", hover_data=['Sample'], color=meta_data_to_test, color_discrete_sequence=color_discrete_sequence)
                fig.update_traces(marker_size=int(nmds_s), mode="markers")
                fig.update_layout(title="Stress=" + str(stress), yaxis_title="NMDS1", xaxis_title="NMDS2")
                fig.update_layout(height=int(height), width=int(width), template=template, showlegend=True, font_size=font_size, title_font_size=font_size)

            ## define output files
            output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_2d.pdf")
            output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_2d.html")
            ## write output files
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))
            ## ask to show file
            if values2['2d_plot'] == True:
                webbrowser.open('file://' + str(output_html))

            ####################################################################################################

            ## plot 3D
            stress = stress_dict[3]

            if values2["draw_mesh"] == True:
                ## create dataframe from NMDS results
                nmds_results_df = pd.DataFrame(nmds_results_dict[3]["nmds_results"], index=[samples])
                nmds_results_df["Sample"] = samples
                nmds_results_df[meta_data_to_test] = Meta_data_table_df[meta_data_to_test].values.tolist()
                nmds_results_df.rename(columns={0:'NMDS1', 1:'NMDS2', 2:'NMDS3'}, inplace=True)

                combinations_list =[]
                for metadata in nmds_results_df[meta_data_to_test]:
                    ## collect all entries for the respective metadata
                    arr = nmds_results_df.loc[nmds_results_df[meta_data_to_test] == metadata][['NMDS1', 'NMDS2', 'NMDS3', meta_data_to_test, "Sample"]].to_numpy()
                    ## create a df for all possible combinations using itertools combinations
                    for entry in list(combinations(arr, 2)):
                        combinations_list.append(list(entry[0]))
                        combinations_list.append(list(entry[1]))

                ## create a dataframe to draw the plot from
                df = pd.DataFrame(combinations_list)
                df.columns = ['NMDS1', 'NMDS2', 'NMDS3', meta_data_to_test, "Sample"]

                ## plot NMDS
                fig = go.Figure()
                ## draw the plot
                fig = px.scatter_3d(df, x="NMDS1", y="NMDS2", z="NMDS3", color=meta_data_to_test, text="Sample", title="textbox", color_discrete_sequence=color_discrete_sequence)
                fig.update_traces(marker_size=int(12), mode="markers+lines", line=dict(width=1))
                fig.update_layout(height=int(height), width=int(width), template=template, title="Stress=" + str(stress), showlegend=True, font_size=font_size, title_font_size=font_size)
                fig.update_layout(scene = dict(xaxis_title="NMDS1",yaxis_title="NMDS2",zaxis_title="NMDS3"))

            else:
                ## create dataframe from NMDS results
                nmds_results_df = pd.DataFrame(nmds_results_dict[3]["nmds_results"], index=[samples])
                nmds_results_df["Sample"] = samples
                nmds_results_df[meta_data_to_test] = Meta_data_table_df[meta_data_to_test].values.tolist()
                nmds_results_df.rename(columns={0:'NMDS1', 1:'NMDS2', 2:'NMDS3'}, inplace=True)
                ## plot NMDS
                fig = go.Figure()
                ## draw the plot
                fig = px.scatter_3d(nmds_results_df, x="NMDS1", y="NMDS2", z="NMDS3", color=meta_data_to_test, color_discrete_sequence=color_discrete_sequence)
                fig.update_traces(marker_size=int(12), mode="markers", line=dict(width=1))
                fig.update_layout(height=int(height), width=int(width), template=template, title="Stress=" + str(stress), showlegend=True, font_size=font_size, title_font_size=font_size)
                fig.update_layout(scene = dict(xaxis_title="NMDS1",yaxis_title="NMDS2",zaxis_title="NMDS3"))


            ## define output files
            output_pdf = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_3d.pdf")
            output_html = Path(str(dirName) + "/" + meta_data_to_test + "_" + taxon_title + "_3d.html")
            ## write output files
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))
            ## ask to show file
            if values2['3d_plot'] == True:
                webbrowser.open('file://' + str(output_html))

            ####################################################################################################

            ## print closing text
            closing_text = "NMDS plots are found in: " + str(path_to_outdirs) + "/NMDS_plots/"
            sg.Popup(closing_text, title="Finished", keep_on_top=True)
            ## write log file
            from taxontabletools.create_log import ttt_log
            ttt_log("nmds analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)
