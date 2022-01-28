import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path
import webbrowser
import plotly.graph_objects as go
import PySimpleGUI as sg
from scipy import stats
import itertools

def alpha_diversity_scatter_plot(TaXon_table_xlsx, meta_data_to_test, width, heigth, scatter_size, taxonomic_level, path_to_outdirs, template, theme, font_size, color_discrete_sequence):

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

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

    categories = Meta_data_table_df[meta_data_to_test].tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.PopupError("The meta data has to differ between samples!", title="Error")
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        ## define variables
        samples = Meta_data_table_samples
        OTU_abundances_dict = {}
        samples_metadata_list = []

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

        # remove samples that do not fit the format
        for i, sample in enumerate(samples):
            meta_data = str(Meta_data_table_df.loc[i][meta_data_to_test])
            samples_metadata_list.append(meta_data)

        #################################
        # Calculate Alpha diversity measurements (= observed_otus)

        observed_otus_dict = {}
        samples_dict = {}

        for i, sample in enumerate(samples):
            observed_otus = len([taxon for taxon in TaXon_table_df[sample].values.tolist() if taxon != 0])
            category = samples_metadata_list[i]
            if category not in observed_otus_dict.keys():
                observed_otus_dict[category] = [observed_otus]
                samples_dict[category] = [sample]
            else:
                observed_otus_dict[category] = observed_otus_dict[category] + [observed_otus]
                samples_dict[category] = samples_dict[category] + [sample]

        ########################################
        # create the plot
        title = "# " + taxon_title

        ## create a list of colors that is longer than the list of categories
        ## plotly express does this automatically, plotly go not!
        while True:
            if len(color_discrete_sequence) < len(categories):
                color_discrete_sequence  = color_discrete_sequence + color_discrete_sequence
            else:
                break

        fig = go.Figure()
        for category, color in zip(sorted(set(categories)), color_discrete_sequence):
            fig.add_trace(go.Scatter(x=samples_dict[category], y=observed_otus_dict[category], mode='markers', name=category, marker=dict(color=color, size=int(scatter_size))))
        fig.update_layout(height=int(heigth), width=int(width), template=template, yaxis_title=title, showlegend=True, font_size=font_size, title_font_size=font_size)

        ## finish script
        output_pdf = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_scatter_plot.pdf")
        output_html = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_scatter_plot.html")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        sg.Popup("Alpha diversity estimate are found in", path_to_outdirs, "/Alpha_diversity/", title="Finished", keep_on_top=True)
        from taxontabletools.create_log import ttt_log
        ttt_log("alpha diversity scatter", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)

    else:
        sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)

def alpha_diversity_boxplot(TaXon_table_xlsx, meta_data_to_test, width, heigth, taxonomic_level, path_to_outdirs, template, theme, font_size, color_discrete_sequence):

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

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

    categories = Meta_data_table_df[meta_data_to_test].tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.PopupError("The meta data has to differ between samples!", title="Error")
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        ## define variables
        samples = Meta_data_table_samples
        OTU_abundances_dict = {}
        samples_metadata_list = []

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

        # remove samples that do not fit the format
        for i, sample in enumerate(samples):
            meta_data = str(Meta_data_table_df.loc[i][meta_data_to_test])
            samples_metadata_list.append(meta_data)

        #################################
        # Calculate Alpha diversity measurements (= observed_otus)

        observed_otus_dict = {}

        for i, sample in enumerate(samples):
            observed_otus = len([taxon for taxon in TaXon_table_df[sample].values.tolist() if taxon != 0])
            category = samples_metadata_list[i]
            if category not in observed_otus_dict.keys():
                observed_otus_dict[category] = [observed_otus]
            else:
                observed_otus_dict[category] = observed_otus_dict[category] + [observed_otus]

        ########################################
        ## kurskal test
        ## for samples of an independant dataset
        kruskal_list = []
        for a, b in itertools.combinations(observed_otus_dict.keys(), 2):
            res = stats.kruskal(observed_otus_dict[a], observed_otus_dict[b])
            len_a = len(observed_otus_dict[a])
            len_b = len(observed_otus_dict[b])
            if res.pvalue <= 0.05:
                if len_a >= 5 and len_b >= 5:
                    kruskal_list = kruskal_list + [[a, len_a, b, len_b, round(res.statistic, 2), round(res.pvalue, 4), "*"]]
                else:
                    kruskal_list = kruskal_list + [[a, len_a, b, len_b, round(res.statistic, 2), round(res.pvalue, 4), "*, n > 5"]]
            else:
                kruskal_list = kruskal_list + [[a, len_a, b, len_b, round(res.statistic, 2), round(res.pvalue, 4), ""]]

        kruskal_df = pd.DataFrame(kruskal_list)
        kruskal_df.columns = ["ID 1", "len 1", "ID 2", "len 2", "H", "p", "#"]


        ########################################
        # create the plot
        title = "# " + taxon_title

        ## create a list of colors that is longer than the list of categories
        ## plotly express does this automatically, plotly go not!
        while True:
            if len(color_discrete_sequence) < len(categories):
                color_discrete_sequence  = color_discrete_sequence + color_discrete_sequence
            else:
                break

        fig = go.Figure()
        for category, color in zip(sorted(set(categories)), color_discrete_sequence):
            fig.add_trace(go.Box(y=observed_otus_dict[category], name=category, marker_color=color, marker_line_color="Black", marker_line_width=0.2, opacity=opacity_value))
        fig.update_layout(height=int(heigth), width=int(width), template=template, yaxis_title=title, showlegend=False, font_size=font_size, title_font_size=font_size)

        ## finish script
        output_pdf = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_boxplot.pdf")
        output_html = Path(str(path_to_outdirs) + "/" + "Alpha_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_" + taxon_title + "_boxplot.html")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        sg.Popup("Alpha diversity estimate are found in", path_to_outdirs, "/Alpha_diversity/", title="Finished", keep_on_top=True)
        from taxontabletools.create_log import ttt_log
        ttt_log("alpha diversity boxplot", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)

    else:
        sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)
