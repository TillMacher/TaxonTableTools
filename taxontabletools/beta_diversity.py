import plotly.figure_factory as ff
import numpy as np
import pandas as pd
from skbio.diversity import beta_diversity
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram, linkage
import PySimpleGUI as sg
import webbrowser
from skbio.stats.distance import anosim
import plotly.express as px
from taxontabletools.taxontable_manipulation import strip_metadata


def beta_diversity_heatmap(TaXon_table_xlsx, width, heigth, cmap, meta_data_to_test, taxonomic_level, path_to_outdirs, template, font_size, diss_metric):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
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

    metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level

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
        ## calculate dissimilarity distances
        dissimilarity_dm = beta_diversity(diss_metric, data, samples)

        anosim_results = anosim(dissimilarity_dm, metadata_list, permutations=999)
        anosim_r = round(anosim_results['test statistic'], 5)
        anosim_p = anosim_results['p-value']
        textbox = "Anosim (" + meta_data_to_test + ", " + taxon_title + ")<br>" + "R = " + str(anosim_r) + "<br>" + "p = " + str(anosim_p)

        matrix = dissimilarity_dm.data
        matrix_df = pd.DataFrame(matrix)
        matrix_df.columns = samples
        matrix_df.index = samples

        # create plot
        color_label = diss_metric + " distance"
        fig = px.imshow(matrix, x=samples,y=samples, color_continuous_scale=cmap, labels=dict(color=color_label))
        fig.update_layout(height=int(heigth), width=int(width), template=template, showlegend=True, title=textbox, font_size=font_size, title_font_size=font_size)

        # finish script
        output_pdf = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_" + diss_metric + ".pdf")
        output_html = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_" + diss_metric + ".html")
        output_xlsx = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_" + diss_metric + ".xlsx")
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

def betadiv_clustering(TaXon_table_xlsx, height, width, threshold, betadiv_linkage, taxonomic_level, path_to_outdirs, template, font_size, diss_metric):

    ## import table
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_df = strip_metadata(TaXon_table_df)

    ## create a y axis title text
    taxon_title = taxonomic_level.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

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
    dissimilarity_dm = beta_diversity(diss_metric, data, samples)
    ## convert to distance matrix
    X1 = dissimilarity_dm.data
    matrix_df = pd.DataFrame(X1)
    matrix_df.columns = samples
    matrix_df.index = samples
    ## convert to 2D array
    X2 = dissimilarity_dm.condensed_form()
    ## cluster dendrogram
    fig = ff.create_dendrogram(X1, labels=samples, color_threshold=float(threshold), orientation="left", linkagefun=lambda x: linkage(X2, betadiv_linkage, metric=diss_metric))
    fig.update_yaxes(ticks="")
    fig.update_xaxes(title="A")
    title = str(diss_metric) + " distance"
    fig.update_layout(xaxis_title=title, height=int(height), width=int(width), template=template, font_size=font_size, title_font_size=font_size)

    # finish script
    output_pdf = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxon_title + "_dendrogram_" + diss_metric + ".pdf")
    output_html = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxon_title + "_dendrogram_" + diss_metric + ".html")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + taxon_title + "_dendrogram_" + diss_metric + ".xlsx")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))
    matrix_df.to_excel(output_xlsx)

    ## ask to show plot
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## write to log file
    sg.Popup(diss_metric + " clustering dendrograms are found in", path_to_outdirs, "/Beta_diversity/", title="Finished", keep_on_top=True)
    from taxontabletools.create_log import ttt_log
    ttt_log(diss_metric + " clustering", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)
