def parcat_analysis(TaXon_table_xlsx, path_to_outdirs, template, height, width, meta_data_to_test, plotly_colors, available_taxonomic_levels_list, taxonomic_level, theme, font_size, color_discrete_sequence):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from pathlib import Path
    import webbrowser, random

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

    ## create a y axis title text
    taxon_title = taxonomic_level

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    if len(set(metadata_list)) == 1:
        sg.PopupError("Please choose more than one meta data category.")
    else:

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

            ## create a list of samples for each category
            category_dict = {}
            for sample, category in zip(TaXon_table_samples, metadata_list):
                if category not in category_dict.keys():
                    category_dict[category] = [sample]
                else:
                    category_dict[category] = category_dict[category] + [sample]

            ## create a color dict
            color_dict = {}
            color_discrete_sequence = color_discrete_sequence * len(TaXon_table_samples)
            for i, category in enumerate(sorted(set(metadata_list))):
                color_dict[category] = color_discrete_sequence[i]

            ## sum up the reads of all samples of a respective category
            ## store them in a dict
            sum_dict = {}
            for key, value in category_dict.items():
                sum_dict[key] = TaXon_table_df[value].sum(axis=1).values.tolist()

            ## collect all available taxa
            taxa = TaXon_table_df[taxonomic_level].values.tolist()

            ## store taxa and labels
            taxon_list = []
            label_list = []
            color_list = []

            ## loop through all categories
            for category, reads in sum_dict.items():
                ## loop through all taxa (in terms of read numbers)
                for i, entry  in enumerate(reads):
                    if entry != 0:
                        taxon_list.append(taxa[i])
                        label_list.append(category)
                        color_list.append(color_dict[category])

            if (taxonomic_level == "Species" or taxonomic_level == "Genus"):
                taxon_list = ["<i>" + taxon + "</i>" for taxon in taxon_list]

            ## create the parcats
            fig = go.Figure()
            fig = go.Figure(go.Parcats(dimensions=[{'label': taxon_title, 'values': taxon_list},{'label': meta_data_to_test, 'values': label_list}],line={'color': color_list, 'shape': 'hspline'}))

            fig.update_layout(height=int(height), width=int(width), template=template, yaxis_title=meta_data_to_test, showlegend=False, font_size=font_size, title_font_size=font_size)
            fig.update_layout(margin=dict(l=200, r=200, t=50, b=50))

            # finish script
            output_pdf = Path(str(path_to_outdirs) + "/" + "ParCat_plots" + "/" + TaXon_table_xlsx.stem + "_" + taxon_title + "_" + meta_data_to_test + "_parcat.pdf")
            output_html = Path(str(path_to_outdirs) + "/" + "ParCat_plots" + "/" + TaXon_table_xlsx.stem + "_" + taxon_title + "_" + meta_data_to_test + "_parcat.html")
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))

            ## ask to show plot
            answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
            if answer == "Yes":
                webbrowser.open('file://' + str(output_html))

            ## write to log file
            sg.Popup("Parallel category plots are found in", path_to_outdirs, "/ParCat_plots/", title="Finished", keep_on_top=True)
            from taxontabletools.create_log import ttt_log
            ttt_log("parallel category analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)

        else:
            sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)
