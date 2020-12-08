def read_proportions_heatmap(TaXon_table_xlsx, taxonomic_level, path_to_outdirs, width_value, height_value, template):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import plotly.express as px
    from pathlib import Path
    import os

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)
    samples_list = TaXon_table_df.columns.tolist()[10:]
    Species_read_proportion_dict = {}

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[samples_list].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        sg.Popup("Please do not use presence absence data!", title=("Error"))
        raise RuntimeError

    if taxonomic_level != "OTUs":
        ## replace nan with the best hit
        taxon_levels_dict = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        value_taxonomic_level = taxon_levels_dict[taxonomic_level]
        best_hit_list = []
        for taxon in TaXon_table_df[list(taxon_levels_dict.keys())].values.tolist():
            ## human readable range => e.g. from 5 to 0 for species level
            for test in range(value_taxonomic_level-1,-1,-1):
                if taxon[test] != "nan":
                    best_hit_list.append(taxon[test])
                    break
        TaXon_table_df[taxonomic_level] = best_hit_list
    else:
        taxonomic_level = "IDs"


    ##############################################################################
    ## create a subfolder for better sorting and overview
    dirName = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem + "/")
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    output_pdf = Path(str(dirName) + "/" + taxonomic_level + "_heatmap.pdf")
    output_html = Path(str(dirName) + "/" + taxonomic_level + "_heatmap.html")
    output_xlsx = Path(str(dirName) + "/" + taxonomic_level + "_heatmap.xlsx")

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_list) + 1
    ############################################################################

    TaXon_table_df_2 = ""

    for sample in samples_list:
        df = TaXon_table_df[['IDs', "Phylum", "Class", "Order", "Family", "Genus", "Species", sample]]
        df_2 = df[[sample]]/df[[sample]].sum()
        df = df.assign(perc=df_2.values)
        df["perc"] = df.groupby([taxonomic_level])['perc'].transform('sum')
        df_3 = df.drop_duplicates(subset=[taxonomic_level, 'perc'])
        df_3 = df_3.drop([sample], axis=1)
        df_3 = df_3.rename(columns={"perc": sample})
        if TaXon_table_df_2 is "":
            TaXon_table_df_2 = df_3
        else:
            TaXon_table_df_2 = TaXon_table_df_2.join(df_3[[sample]])

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

    ## create plot
    ## ask if as subplot shall be generated
    plot_df = TaXon_table_df_2[samples_list]
    plot_df.index = TaXon_table_df_2[taxonomic_level]

    fig = px.imshow(plot_df)
    fig.update_traces(showlegend=False)
    fig.update_layout(width=int(width_value), height=int(height_value), template=template)

    # finish script
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        fig.show()
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    closing_text = "\n" + "Read proportion plots are found in: " + str(path_to_outdirs) + "/Read_proportion_plots/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("read proportions heatmap", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)

def read_proportions_bar(TaXon_table_xlsx, taxonomic_level, path_to_outdirs, width_value, height_value, template):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import plotly.express as px
    from pathlib import Path
    import os

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)
    samples_list = TaXon_table_df.columns.tolist()[10:]
    Species_read_proportion_dict = {}

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[samples_list].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        sg.Popup("Please do not use presence absence data!", title=("Error"))
        raise RuntimeError

    if taxonomic_level != "OTUs":
        ## replace nan with the best hit
        taxon_levels_dict = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        value_taxonomic_level = taxon_levels_dict[taxonomic_level]
        best_hit_list = []
        for taxon in TaXon_table_df[list(taxon_levels_dict.keys())].values.tolist():
            ## human readable range => e.g. from 5 to 0 for species level
            for test in range(value_taxonomic_level-1,-1,-1):
                if taxon[test] != "nan":
                    best_hit_list.append(taxon[test])
                    break
        TaXon_table_df[taxonomic_level] = best_hit_list
    else:
        taxonomic_level = "IDs"

    ##############################################################################
    ## create a subfolder for better sorting and overview
    dirName = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem + "/")
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    output_pdf = Path(str(dirName) + "/" + taxonomic_level + "_bar.pdf")
    output_html = Path(str(dirName) + "/" + taxonomic_level + "_bar.html")
    output_xlsx = Path(str(dirName) + "/" + taxonomic_level + "_bar.xlsx")

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_list) + 1
    ############################################################################

    TaXon_table_df_2 = ""

    for sample in samples_list:
        df = TaXon_table_df[['IDs', "Phylum", "Class", "Order", "Family", "Genus", "Species", sample]]
        df_2 = df[[sample]]/df[[sample]].sum()
        df = df.assign(perc=df_2.values*100)
        df["perc"] = df.groupby([taxonomic_level])['perc'].transform('sum')
        df_3 = df.drop_duplicates(subset=[taxonomic_level, 'perc'])
        df_3 = df_3.drop([sample], axis=1)
        df_3 = df_3.rename(columns={"perc": sample})
        if TaXon_table_df_2 is "":
            TaXon_table_df_2 = df_3
        else:
            TaXon_table_df_2 = TaXon_table_df_2.join(df_3[[sample]])

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    ## create dataframe for plotly express
    plot_df = ''
    for sample in samples_list:
        if plot_df is '':
            plot_df = pd.DataFrame([[sample] + entry for entry in TaXon_table_df_2[[taxonomic_level, sample]].values.tolist()], columns=["Sample", "Taxon", "Reads"])
        else:
            df = pd.DataFrame([[sample] + entry for entry in TaXon_table_df_2[[taxonomic_level, sample]].values.tolist()], columns=["Sample", "Taxon", "Reads"])
            plot_df = plot_df.append(df)

    n_taxa = len(TaXon_table_df_2[taxonomic_level].values.tolist())
    plot_df["Color"] = list(np.linspace(0,100,n_taxa)) * len(samples_list)

    fig = px.bar(plot_df, x="Sample", y="Reads", color="Taxon", color_continuous_scale='Cividis_r', labels={"Color": "Taxon"})
    fig.update_layout(barmode='stack', width=int(width_value), height=int(height_value), template=template)
    fig.update_yaxes(title_text="Reads (%)")
    fig.update_xaxes(title_text="")

    answer = sg.PopupYesNo('Open html?', keep_on_top=True)
    if answer == "Yes":
        fig.show()
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))
    closing_text = "Read proportion plot is found under:\n" + '/'.join(str(output_pdf).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
    from taxontabletools.create_log import ttt_log
    ttt_log("read proportions bar plot", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)

def read_proportions_pie(TaXon_table_xlsx, taxonomic_level, path_to_outdirs, width_value, height_value, template):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from pathlib import Path
    import os

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)
    samples_list = TaXon_table_df.columns.tolist()[10:]
    Species_read_proportion_dict = {}

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[samples_list].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        sg.Popup("Please do not use presence absence data!", title=("Error"))
        raise RuntimeError

    if taxonomic_level != "OTUs":
        ## replace nan with the best hit
        taxon_levels_dict = {"Phylum": 1, "Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
        value_taxonomic_level = taxon_levels_dict[taxonomic_level]
        best_hit_list = []
        for taxon in TaXon_table_df[list(taxon_levels_dict.keys())].values.tolist():
            ## human readable range => e.g. from 5 to 0 for species level
            for test in range(value_taxonomic_level-1,-1,-1):
                if taxon[test] != "nan":
                    best_hit_list.append(taxon[test])
                    break
        TaXon_table_df[taxonomic_level] = best_hit_list
    else:
        taxonomic_level = "IDs"

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples_list) + 1
    ############################################################################

    TaXon_table_df_2 = ""

    for sample in samples_list:
        df = TaXon_table_df[['IDs', "Phylum", "Class", "Order", "Family", "Genus", "Species", sample]]
        df_2 = df[[sample]]/df[[sample]].sum()
        df = df.assign(perc=df_2.values)
        df["perc"] = df.groupby([taxonomic_level])['perc'].transform('sum')
        df_3 = df.drop_duplicates(subset=[taxonomic_level, 'perc'])
        df_3 = df_3.drop([sample], axis=1)
        df_3 = df_3.rename(columns={"perc": sample})
        if TaXon_table_df_2 is "":
            TaXon_table_df_2 = df_3
        else:
            TaXon_table_df_2 = TaXon_table_df_2.join(df_3[[sample]])

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    ## create dataframe for plot
    plot_df = TaXon_table_df_2[samples_list]
    plot_df.index = TaXon_table_df_2[taxonomic_level]

    ##############################################################################
    ## create a subfolder for better sorting and overview
    dirName = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem + "/")
    dirName_samples = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem + "/samples")
    if not os.path.exists(dirName):
        os.mkdir(dirName)
    if not os.path.exists(dirName_samples):
        os.mkdir(dirName_samples)

    ## read abundance pie chart per sample
    for sample in samples_list:
        sample_df = plot_df.loc[plot_df[sample] > 0.0, [sample]]
        labels = sample_df.index.tolist()
        values = sample_df[sample].values.tolist()
        fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.3)])
        fig.update_layout(template=template, title=sample, annotations=[dict(text=taxonomic_level, x=0.5, y=0.5, showarrow=False)])
        fig.update_traces(textposition='inside')
        fig.update_layout(width=int(width_value), height=int(height_value), template=template)

        output_pdf = Path(str(dirName_samples) + "/" + sample + "_" + taxonomic_level + "_pie.pdf")
        output_html = Path(str(dirName_samples) + "/" + sample + "_" + taxonomic_level + "_pie.html")

        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

    ## main read abundance pie chart
    main_df = pd.DataFrame(TaXon_table_df[taxonomic_level].values.tolist(), list(TaXon_table_df[samples_list].sum(axis=1)), columns=["Taxon"])
    main_df["Reads"] = main_df.index
    df_2 = main_df["Reads"]/main_df["Reads"].sum()
    main_df = main_df.assign(perc=df_2.values*100)

    fig = go.Figure(data=[go.Pie(labels=main_df["Taxon"], values=main_df["perc"], hole=.3)])
    fig.update_traces(textposition='inside')
    fig.update_layout(annotations=[dict(text=taxonomic_level, x=0.5, y=0.5, showarrow=False)])
    fig.update_layout(width=int(width_value), height=int(height_value), template=template)

    output_pdf = Path(str(dirName) + "/" + taxonomic_level + "_pie.pdf")
    output_html = Path(str(dirName) + "/" + taxonomic_level + "_pie.html")
    output_xlsx = Path(str(dirName) + "/" + taxonomic_level + "_pie.xlsx")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## finish script
    answer = sg.PopupYesNo('Open html?', keep_on_top=True)
    if answer == "Yes":
        fig.show()
    closing_text = "Read proportion plot is found under:\n" + '/'.join(str(output_pdf).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
    from taxontabletools.create_log import ttt_log
    ttt_log("read proportions pie chart", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)
