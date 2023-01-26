from pathlib import Path
import PySimpleGUI as sg
import statsmodels.api as sm
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import scipy.stats, math
import numpy as np
from scipy import stats
from statistics import mean
from scipy.spatial import distance
import plotly.figure_factory as ff
from statistics import mean
import webbrowser, os
from taxontabletools.taxontable_manipulation import strip_metadata

def tc_alpha_diversity(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, font_size, color_discrete_sequence):

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = "# " + taxonomic_level
        taxonomic_level = "ID"
    if taxonomic_level == "Phylum":
        taxon_title = "# phyla"
    if taxonomic_level == "Class":
        taxon_title = "# classes"
    if taxonomic_level == "Order":
        taxon_title = "# orders"
    if taxonomic_level == "Family":
        taxon_title = "# families"
    if taxonomic_level == "Genus":
        taxon_title = "# genera"
    if taxonomic_level == "Species":
        taxon_title = "# species"

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## count overall number of taxa
    df_1_taxa = len(list(set([i for i in df1[taxonomic_level].values.tolist() if i != "nan"])))
    df_2_taxa = len(list(set([i for i in df2[taxonomic_level].values.tolist() if i != "nan"])))

    ## count the number of taxa per sample for both tables
    df1_dict, df2_dict = {}, {}

    for sample in df1_samples:
        n_taxa = set([i[0] for i in df1[[taxonomic_level, sample]].values.tolist() if i[0] != "nan" and i[1] != 0])
        df1_dict[sample] = len(list(n_taxa))

    for sample in df2_samples:
        n_taxa = set([i[0] for i in df2[[taxonomic_level, sample]].values.tolist() if i[0] != "nan" and i[1] != 0])
        df2_dict[sample] = len(list(n_taxa))

    ## collect values
    y1 = list(df1_dict.values())
    y2 = list(df2_dict.values())
    text_values1 = list(df1_dict.keys())
    text_values2 = list(df2_dict.keys())

    ## create a boxlplot
    fig = go.Figure()
    fig.add_trace(go.Box(y=y1, name=name_1, text=text_values1, marker_color=color_discrete_sequence[0]))
    fig.add_trace(go.Box(y=y2, name=name_2, text=text_values2, marker_color=color_discrete_sequence[1]))
    fig.update_yaxes(title=taxon_title)
    fig.update_traces(boxpoints='all', jitter=0.5)
    fig.update_yaxes(rangemode="tozero")
    fig.update_layout(width=int(width_value), height=int(height_value), template=template, showlegend=False, font_size=font_size)

    ## create a folder if neccessary
    dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    ## write ouput files
    output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_alpha_diversity.pdf")
    output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_alpha_diversity.html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write to log
    from taxontabletools.create_log import ttt_log
    out = name_1 + ";" + name_2
    ttt_log("Table comparison alpha diversity", "analysis", out, output_pdf.name, "", path_to_outdirs)

def tc_overall_taxon_richness(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, taxonomic_level_2, font_size, color_discrete_sequence):

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = "# " + taxonomic_level
        taxonomic_level = "ID"
    if taxonomic_level == "Phylum":
        taxon_title = "# phyla"
    if taxonomic_level == "Class":
        taxon_title = "# classes"
    if taxonomic_level == "Order":
        taxon_title = "# orders"
    if taxonomic_level == "Family":
        taxon_title = "# families"
    if taxonomic_level == "Genus":
        taxon_title = "# genera"
    if taxonomic_level == "Species":
        taxon_title = "# species"

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## collect all taxa for the second (i.e. higher taxonomic level)
    taxa2_list = sorted(list(set([i for i in df1[taxonomic_level_2].values.tolist() + df2[taxonomic_level_2].values.tolist() if i != "nan"])))

    ## store the number of taxa per taxon (e.g. species per phylum)
    taxa_dict_1, taxa_dict_2 = {}, {}

    ## loop through taxa
    for taxon in taxa2_list:
        taxa_dict_1[taxon] = len(set([i[0] for i in df1[[taxonomic_level, taxonomic_level_2]].values.tolist() if i[0] != "nan" and i[1] == taxon]))
        taxa_dict_2[taxon] = len(set([i[0] for i in df2[[taxonomic_level, taxonomic_level_2]].values.tolist() if i[0] != "nan" and i[1] == taxon]))

    ## collect values
    x = list(taxa_dict_1.values())
    y = list(taxa_dict_2.values())
    text_values = list(taxa_dict_2.keys())

    ## calcuate spearman correlation
    spearman = scipy.stats.spearmanr(x, y)
    spearman_p = spearman[1]

    if spearman_p <= 0.05:
        spearman_rho = str(round(spearman[0], 3)) + "*"
    else:
        spearman_rho = str(round(spearman[0], 3))

    ## plot data
    fig = go.Figure()
    ## add scatter plot
    fig.add_trace(go.Scatter(name="", x=x, y=y, mode='markers', text=text_values, marker=dict(color=color_discrete_sequence[0], size=15)))
    ## add best fit
    bestfit = sm.OLS(y,sm.add_constant(x)).fit().fittedvalues
    fig.add_trace(go.Scatter(name="rho=" + spearman_rho, x=x, y=bestfit, mode='lines', opacity=0.2, marker=dict(color="black")))
    fig.update_yaxes(title=taxon_title + " (" + name_2 + ")")
    fig.update_xaxes(title=taxon_title + " (" + name_1 + ")")
    title_text = taxonomic_level + " per " + taxonomic_level_2.lower() + " (rho=" + spearman_rho + ")"
    fig.update_layout(title=title_text, width=int(width_value), height=int(height_value), template=template, showlegend=False, font_size=font_size)

    ## create a folder if neccessary
    dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    ## write ouput files
    output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "per" + taxonomic_level_2 + "_taxonomic_richness_scatter.pdf")
    output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "per" + taxonomic_level_2 + "_taxonomic_richness_scatter.html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write to log
    from taxontabletools.create_log import ttt_log
    out = name_1 + ";" + name_2
    ttt_log("Table comparison taxonomic richness scatter", "analysis", out, output_pdf.name, "", path_to_outdirs)

def tc_seperate_taxon_richness(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, taxonomic_level_2, font_size, color_discrete_sequence, min_show):

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## collect all taxa for the second (i.e. higher taxonomic level)
    taxa2_list = sorted(list(set([i for i in df1[taxonomic_level_2].values.tolist() + df2[taxonomic_level_2].values.tolist() if i != "nan"])))

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = "# " + taxonomic_level
        taxonomic_level = "ID"
    if taxonomic_level == "Phylum":
        taxon_title = "# phyla"
    if taxonomic_level == "Class":
        taxon_title = "# classes"
    if taxonomic_level == "Order":
        taxon_title = "# orders"
    if taxonomic_level == "Family":
        taxon_title = "# families"
    if taxonomic_level == "Genus":
        taxon_title = "# genera"
    if taxonomic_level == "Species":
        taxon_title = "# species"

    non_shown = []
    y1_dict = {}
    y2_dict = {}


    ## loop through taxa
    for taxon in taxa2_list:

        y1, y2, text_values1, text_values2 = [], [], [], []
        min_show = int(min_show)

        for sample in df1_samples:
            n1 = len(set([i[0] for i in df1[[taxonomic_level, taxonomic_level_2, sample]].values.tolist() if i[0] != "nan" and i[1] == taxon and i[2] != 0]))
            y1.append(n1)
            text_values1.append(sample)

        for sample in df2_samples:
            n2 = len(set([i[0] for i in df2[[taxonomic_level, taxonomic_level_2, sample]].values.tolist() if i[0] != "nan" and i[1] == taxon and i[2] != 0]))
            y2.append(n2)
            text_values2.append(sample)

        if max(y1+y2) >= min_show:
            y1_dict[taxon] = [y1, text_values1]
            y2_dict[taxon] = [y2, text_values2]
        else:
            non_shown.append(taxon)

    ## plot
    ## create a boxlplot
    n_cols = 2
    n_rows = math.ceil(len(y1_dict.keys()) / n_cols)
    fig = make_subplots(rows=n_rows, cols=n_cols, vertical_spacing = 0.1, horizontal_spacing = 0.1, subplot_titles=list(y1_dict.keys()))
    col_counter = 1
    row_counter = 1

    for taxon in y1_dict.keys():
        fig.add_trace(go.Box(y=y1_dict[taxon][0], name=name_1, text=y1_dict[taxon][1], marker_color=color_discrete_sequence[0], marker_line_color="Black", marker_line_width=0.2), col=col_counter, row=row_counter)
        fig.add_trace(go.Box(y=y2_dict[taxon][0], name=name_2, text=y2_dict[taxon][1], marker_color=color_discrete_sequence[1], marker_line_color="Black", marker_line_width=0.2), col=col_counter, row=row_counter)
        fig.update_yaxes(rangemode="tozero", col=col_counter, row=row_counter)
        fig.update_traces(boxpoints='all', jitter=0.4, col=col_counter, row=row_counter)

        if col_counter == 1:
            fig.update_yaxes(title=taxon_title, col=col_counter, row=row_counter)

        ## iterate through rows and columns
        col_counter += 1
        if col_counter == n_cols+1:
            col_counter = 1
            row_counter += 1


    fig.update_layout(showlegend=False, width=int(width_value), height=int(height_value), template=template, font_size=font_size)

    ## create a folder if neccessary
    dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    ## write ouput files
    output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_per_" + taxonomic_level_2 + "_taxonomic_richness_boxplots.pdf")
    output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_per_" + taxonomic_level_2 + "_taxonomic_richness_boxplots.html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = str(len(non_shown)) + " taxa fell below the treshold and are not shown.\n\nPlots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write to log
    from taxontabletools.create_log import ttt_log
    out = name_1 + ";" + name_2
    ttt_log("Table comparison taxonomic richness boxplot", "analysis", out, output_pdf.name, "", path_to_outdirs)

def tc_pairwise_sample_comparison(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, font_size, color_discrete_sequence):

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"
    if taxonomic_level == "Phylum":
        taxon_title = "phyla"
    if taxonomic_level == "Class":
        taxon_title = "classes"
    if taxonomic_level == "Order":
        taxon_title = "orders"
    if taxonomic_level == "Family":
        taxon_title = "families"
    if taxonomic_level == "Genus":
        taxon_title = "genera"
    if taxonomic_level == "Species":
        taxon_title = "species"

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## collect all taxa for the second (i.e. higher taxonomic level)
    taxa_blacklist = sorted(list(set([i for i in df1[taxonomic_level].values.tolist() + df2[taxonomic_level].values.tolist() if i != "nan"])))

    ## only proceed if both files have the same samples!
    if sorted(df1_samples) != sorted(df2_samples):
        sg.PopupError("Samples do not match between files!")
    else:

        jaccard_dict = {}
        jaccard_values = []

        for sample in df1_samples:

            n1 = set([i[0] for i in df1[[taxonomic_level, sample]].values.tolist() if i[0] != 'nan' and i[1] != 0])
            n2 = set([i[0] for i in df2[[taxonomic_level, sample]].values.tolist() if i[0] != 'nan' and i[1] != 0])

            col1 = []
            col2 = []
            for taxon in list(set(list(n1) + list(n2))):
                if taxon in n1:
                    col1.append(1)
                else:
                    col1.append(0)
                if taxon in n2:
                    col2.append(1)
                else:
                    col2.append(0)

            dist = round(distance.jaccard(col1, col2), 2)

            jaccard_dict[sample] = str(dist)
            jaccard_values.append(dist)

        avg_jaccard = round(mean(jaccard_values),2)
        min_jaccard = min(jaccard_values)
        max_jaccard = max(jaccard_values)

        jaccard_dict = dict(sorted(jaccard_dict.items(), key=lambda item: item[1]))

        fig = make_subplots(rows=1, cols=3, column_widths=[0.8, 0.2, 0.2], vertical_spacing = 0.05, horizontal_spacing = 0.05)

        y_table_1_list = []
        y_shared_list = []
        y_table_2_list = []
        x_samples = []
        y_n_taxa = []

        for sample in jaccard_dict.keys():

            table_1_list = []
            shared_list = []
            table_2_list = []

            n1 = set([i[0] for i in df1[[taxonomic_level, sample]].values.tolist() if i[1] != 0 and i[0] != "nan"])
            n2 = set([i[0] for i in df2[[taxonomic_level, sample]].values.tolist() if i[1] != 0 and i[0] != "nan"])

            n_taxa = len(set(list(n1) + list(n2)))

            if n_taxa != 0:
                if len(n1) != 0 and len(n2) != 0:
                    table_1_list.append(len(n1 - n2))
                    shared_list.append(len(n1 & n2))
                    table_2_list.append(len(n2 - n1))
                elif len(n1) == 0:
                    table_1_list.append(0)
                    shared_list.append(0)
                    table_2_list.append(len(n2))
                elif len(n2) == 0:
                    table_1_list.append(len(n1))
                    shared_list.append(0)
                    table_2_list.append(0)

            table_1_only = sum(table_1_list) / n_taxa * 100
            shared = sum(shared_list) / n_taxa * 100
            table_2_only = sum(table_2_list) / n_taxa * 100

            y_table_1_list.append(table_1_only)
            y_shared_list.append(shared)
            y_table_2_list.append(table_2_only)
            x_samples.append(sample)
            y_n_taxa.append(n_taxa)


        ## add the values used for the plot to a dataframe
        df_out = pd.DataFrame()
        df_out["Sample"] = x_samples
        df_out["Only "+ name_1] = y_table_1_list
        df_out["Shared"] = y_shared_list
        df_out["Only " + name_2] = y_table_2_list

        ## Shared / Only barchart
        fig.add_trace(go.Bar(name="Shared", orientation='h', y=x_samples, x=y_shared_list, marker_color="rgb(102,194,164)"), row=1, col=1)
        fig.add_trace(go.Bar(name=name_1, orientation='h', y=x_samples, x=y_table_1_list, marker_color="rgb(141,160,203)"), row=1, col=1)
        fig.add_trace(go.Bar(name=name_2, orientation='h', y=x_samples, x=y_table_2_list, marker_color="rgb(252,141,98)"), row=1, col=1)
        fig.update_layout(barmode='stack', showlegend=False, width=int(width_value), height=int(height_value), template=template, title="", font_size=font_size)
        fig.update_yaxes(tickmode = 'linear', showgrid=False, row=1, col=1)
        fig.update_xaxes(title=taxon_title + " (%)", showgrid=True, row=1, col=1)

        ## Jaccard plot
        y = list(jaccard_dict.keys())
        x = [float(i) for i in list(jaccard_dict.values())]
        df_out["Jaccard dissimilarity"] = x
        fig.add_trace(go.Bar(y=y, x=x, name="Jaccard", orientation='h', marker_color="lightgrey"), row=1, col=2)
        fig.update_yaxes(tickmode = 'linear', showticklabels=False, showgrid=False, row=1, col=2)
        fig.update_xaxes(title="jaccard diss.", showgrid=True, range=[0, 1], autorange=False, tick0=0, dtick=0.5, row=1, col=2)

        y = x_samples
        x = y_n_taxa
        df_out["# taxa"] = x
        fig.add_trace(go.Bar(y=y, x=x, name="taxa", orientation='h', marker_color="lightgrey"), row=1, col=3)
        fig.update_yaxes(tickmode = 'linear', showticklabels=False, showgrid=False, row=1, col=3)
        fig.update_xaxes(title="# " + taxon_title, showgrid=True, row=1, col=3)

        min_shared = round(min(y_shared_list),2)
        avg_shared = round(mean(y_shared_list),2)
        max_shared = round(max(y_shared_list),2)

        df_out = df_out.sort_values("Jaccard dissimilarity", ascending=False)

        ## create a folder if neccessary
        dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

        ## write ouput files
        output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_overview_chart.pdf")
        output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_overview_chart.html")
        output_xlsx = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_overview_chart.xlsx")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))
        df_out.to_excel(output_xlsx, index=False)

        ## ask to show file
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## print closing text
        closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        out = name_1 + ";" + name_2
        ttt_log("Table comparison overview", "analysis", out, output_pdf.name, "", path_to_outdirs)

def tc_pairwise_taxon_comparison(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, font_size, color_discrete_sequence, method):

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## collect all taxa for the second (i.e. higher taxonomic level)
    taxa_blacklist = sorted(list(set([i for i in df1[taxonomic_level].values.tolist() + df2[taxonomic_level].values.tolist() if i != "nan"])))

    levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    ## only proceed if both files have the same samples!
    if sorted(df1_samples) != sorted(df2_samples):
        sg.PopupError("Samples do not match between files!")
    else:
        ## count taxa
        taxon_dict_x, taxon_dict_y = {}, {}

        ## unique taxon list for each sample
        sample_taxa_1, sample_taxa_2 = {}, {}
        n_taxa_1, n_taxa_2 = {}, {}

        ## store information for each sample in the two dicts
        for sample in df1_samples:
            ## collect taxa for each sample
            lst1 = [i[:-1] for i in df1[levels + [sample]].values.tolist() if i[-1] != 0]
            lst2 = [i[:-1] for i in df2[levels + [sample]].values.tolist() if i[-1] != 0]

            ## remove duplicates
            lst1 = [list(x) for x in set(tuple(x) for x in lst1)]
            lst2 = [list(x) for x in set(tuple(x) for x in lst2)]

            ## count overall number of taxa for rel. abundance
            n_taxa_1[sample] = len(lst1)
            n_taxa_2[sample] = len(lst2)

            ## flatten list
            sample_taxa_1[sample] = [item for sublist in lst1 for item in sublist]
            sample_taxa_2[sample] = [item for sublist in lst2 for item in sublist]

        n_significant, values_dict = [], {}

        for taxon in taxa_blacklist:
            ## store x and y values for every taxon
            x_values = []
            y_values = []

            if method == "Relative":
                ## collect information for each sample
                for sample in df1_samples:
                    x_values.append(sample_taxa_1[sample].count(taxon) / n_taxa_1[sample] * 100)
                    y_values.append(sample_taxa_2[sample].count(taxon) / n_taxa_2[sample] * 100)
            else:
                ## collect information for each sample
                for sample in df1_samples:
                    x_values.append(sample_taxa_1[sample].count(taxon))
                    y_values.append(sample_taxa_2[sample].count(taxon))

            ## test for spearman correlation
            spearman = scipy.stats.spearmanr(x_values, y_values)
            if spearman[1] <= 0.05:
                n_significant.append(1)
                R = " (rho=" + str(round(spearman[0], 2)) + "*)"
                values_dict[taxon] = [taxon + R], x_values, y_values, df1_samples
            elif spearman[1] > 0.05:
                n_significant.append(1)
                R = " (rho=" + str(round(spearman[0], 2)) + ")"
                values_dict[taxon] = [taxon + R], x_values, y_values, df1_samples

        ## create the plot
        columns = 3
        rows = math.ceil(len(n_significant) / columns)

        ## create a new figure
        titles = [i[0][0] for i in values_dict.values()]
        fig = make_subplots(rows=rows, cols=columns, subplot_titles=titles)

        col_count = 1
        row_count = 1
        for taxon in values_dict.keys():
            ## collect x values (aka df1 values)
            x = values_dict[taxon][1]
            ## collect x values (aka df2 values)
            y = values_dict[taxon][2]
            ## collect sample names
            text = values_dict[taxon][3]

            ## add scatter plot
            fig.add_trace(go.Scatter(name=taxon, x=x, y=y, mode='markers', text=text, marker_color="Blue"), row=row_count, col=col_count)

            ## add best fit
            bestfit = sm.OLS(y,sm.add_constant(x)).fit().fittedvalues
            fig.add_trace(go.Scatter(name="", x=x, y=bestfit, mode='lines', opacity=0.2, marker=dict(color="black")), row=row_count, col=col_count)

            ## add y axis titles
            if col_count == 1:
                if method == "Relative":
                    y_title = "taxa % (" + name_2 + ")"
                else:
                    y_title = "# taxa (" + name_2 + ")"
                fig.update_yaxes(title=y_title, row=row_count, col=col_count)

            ## add x axis titles
            if row_count == rows:
                if method == "Relative":
                    x_title = "taxa % (" + name_1 + ")"
                else:
                    x_title = "# taxa (" + name_1 + ")"
                fig.update_xaxes(title=x_title, row=row_count, col=col_count)

            ## increase column count
            col_count += 1
            ## increase row count inf max is reached and reset column count
            if col_count > columns:
                col_count = 1
                row_count += 1

        fig.update_layout(barmode='stack', showlegend=False, width=int(width_value), height=int(height_value), template=template, title="", font_size=font_size)

        ## create a folder if neccessary
        dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

        ## write ouput files
        output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_overview_chart.pdf")
        output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "_overview_chart.html")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

        ## ask to show file
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## print closing text
        closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        out = name_1 + ";" + name_2
        ttt_log("Table comparison overview", "analysis", out, output_pdf.name, "", path_to_outdirs)

def tc_taxon_occurrence(TaXon_table_xlsx_1, TaXon_table_xlsx_2, path_to_outdirs, name_1, name_2, width_value, height_value, template, taxonomic_level, taxonomic_level_2, font_size, color_discrete_sequence):

    ## load dataframes and sample names
    df1 = pd.read_excel(TaXon_table_xlsx_1).fillna('nan')
    df1 = strip_metadata(df1)
    df1_samples = df1.columns.tolist()[10:]
    df2 = pd.read_excel(TaXon_table_xlsx_2).fillna('nan')
    df2 = strip_metadata(df2)
    df2_samples = df2.columns.tolist()[10:]

    ## collect all taxa for the second (i.e. higher taxonomic level)
    taxa2_list = sorted(list(set([i for i in df1[taxonomic_level_2].values.tolist() + df2[taxonomic_level_2].values.tolist() if i != "nan"])))

    ## number of unique species
    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = "# " + taxonomic_level
        taxonomic_level = "ID"
    if taxonomic_level == "Phylum":
        taxon_title = "# unique phyla"
    if taxonomic_level == "Class":
        taxon_title = "# unique classes"
    if taxonomic_level == "Order":
        taxon_title = "# unique orders"
    if taxonomic_level == "Family":
        taxon_title = "# unique families"
    if taxonomic_level == "Genus":
        taxon_title = "# unique genera"
    if taxonomic_level == "Species":
        taxon_title = "# unique species"

    taxa_dict = {}
    ## loop through taxa
    for taxon in taxa2_list:
        n1 = set([i[0] for i in df1[[taxonomic_level, taxonomic_level_2]].values.tolist() if i[0] != "nan" and i[1] == taxon])
        n2 = set([i[0] for i in df2[[taxonomic_level, taxonomic_level_2]].values.tolist() if i[0] != "nan" and i[1] == taxon])

        only_n1 = len(n1 - n2)
        only_n2 = len(n2 - n1)

        taxa_dict[taxon] = [only_n1, 0 - only_n2]

    taxa_dict = {k: v for k, v in sorted(taxa_dict.items(), key=lambda item: item[1], reverse=True)}

    fig = go.Figure()
    y1 = [i[0] for i in taxa_dict.values()]
    y2 = [i[1] for i in taxa_dict.values()]
    x = list(taxa_dict.keys())
    fig.add_trace(go.Bar(x=x, y=y2, base=0, name=name_2, marker_color=color_discrete_sequence[1]))
    fig.add_trace(go.Bar(x=x, y=y1, base=0, name=name_1, marker_color=color_discrete_sequence[0]))
    fig.update_layout(showlegend=True, width=int(width_value), height=int(height_value), template=template, title="", font_size=font_size, barmode="stack")
    fig.update_xaxes(tickmode='linear', tickangle=-90)
    fig.update_yaxes(title=taxon_title)

    ## create a folder if neccessary
    dirName = Path(str(path_to_outdirs) + "/Table_comparison/" + name_1 + "_" + name_2 + "/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

    ## write ouput files
    output_pdf = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "per" + taxonomic_level_2 + "_unique_taxa.pdf")
    output_html = Path(str(dirName) + "/" + name_1 + "_" + name_2 + "_" + taxonomic_level + "per" + taxonomic_level_2 + "_unique_taxa.html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Table_comparison/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write to log
    from taxontabletools.create_log import ttt_log
    out = name_1 + ";" + name_2
    ttt_log("Table comparison unique taxa", "analysis", out, output_pdf.name, "", path_to_outdirs)
















#
