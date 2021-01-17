def beta_diversity(TaXon_table_xlsx, beta_w, beta_h, beta_cmap, meta_data_to_test, taxonomic_level, path_to_outdirs, template, font_size):

    import pandas as pd
    import numpy as np
    from skbio.diversity import beta_diversity
    from skbio.stats.distance import anosim
    import plotly.express as px
    from pathlib import Path
    import PySimpleGUI as sg
    import webbrowser

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == len(Meta_data_table_df['Samples'].tolist()):
        sg.Popup("The meta data is unique for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.Popup("The meta data is similar for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        ## collect samples for plot
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
        jaccard_dm = beta_diversity("jaccard", data, samples)

        anosim_results = anosim(jaccard_dm, metadata_list, permutations=999)
        anosim_r = round(anosim_results['test statistic'], 5)
        anosim_p = anosim_results['p-value']
        textbox = "Anosim (" + meta_data_to_test + ", " + taxonomic_level + ")<br>" + "R = " + str(anosim_r) + "<br>" + "p = " + str(anosim_p)

        matrix = jaccard_dm.data
        matrix_df = pd.DataFrame(matrix)
        matrix_df.columns = samples
        matrix_df.index = samples

        # create plot
        fig = px.imshow(matrix, x=samples,y=samples, color_continuous_scale=beta_cmap, labels=dict(color="Jaccard distance"))
        fig.update_layout(height=int(beta_h), width=int(beta_w), template=template, showlegend=True, title=textbox, font_size=font_size, title_font_size=font_size)

        # finish script
        output_pdf = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxonomic_level + "_jc.pdf")
        output_html = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxonomic_level + "_jc.html")
        output_xlsx = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxonomic_level + "_jc.xlsx")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))
        matrix_df.to_excel(output_xlsx)

        ## ask to show plot
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## write to log file
        sg.Popup("Beta diversity estimate are found in", path_to_outdirs, "/Beta_diversity/", title="Finished", keep_on_top=True)
        from taxontabletools.create_log import ttt_log
        ttt_log("beta diversity", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)

    else:
        sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)

def jaccard_clustering(TaXon_table_xlsx, height, width, threshold, taxonomic_level, path_to_outdirs, template, font_size):
    from scipy.cluster.hierarchy import dendrogram, linkage
    import plotly.figure_factory as ff
    import numpy as np
    import pandas as pd
    from skbio.diversity import beta_diversity
    from pathlib import Path
    import PySimpleGUI as sg
    import webbrowser

    ## import table
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")

    ## collect samples for plot
    samples = TaXon_table_df.columns.tolist()[10:]

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
    jaccard_dm = beta_diversity("jaccard", data, samples)
    ## convert to distance matrix
    X1 = jaccard_dm.data
    matrix_df = pd.DataFrame(X1)
    matrix_df.columns = samples
    matrix_df.index = samples
    ## convert to 2D array
    X2 = jaccard_dm.condensed_form()
    ## cluster dendrogram
    fig = ff.create_dendrogram(X1, labels=samples, color_threshold=float(threshold), orientation="left", linkagefun=lambda x: linkage(X2, 'complete', metric='jaccard'))
    fig.update_yaxes(ticks="")
    fig.update_xaxes(title="A")
    fig.update_layout(xaxis_title="Jaccard distance", height=int(height), width=int(width), template=template, font_size=font_size, title_font_size=font_size)

    # finish script
    output_pdf = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxonomic_level + "_dendrogram.pdf")
    output_html = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxonomic_level + "_dendrogram.html")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxonomic_level + "_dendrogram.xlsx")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))
    matrix_df.to_excel(output_xlsx)

    ## ask to show plot
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## write to log file
    sg.Popup("Jaccard clustering dendrograms are found in", path_to_outdirs, "/Beta_diversity/", title="Finished", keep_on_top=True)
    from taxontabletools.create_log import ttt_log
    ttt_log("jaccard clustering", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)
