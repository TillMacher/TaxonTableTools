def site_occupancy_barchart(TaXon_table_xlsx, meta_data_to_test, taxonomic_level, path_to_outdirs, x_site_occ, y_site_occ, template, theme, font_size):

    import os, webbrowser
    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path
    import plotly.graph_objects as go
    import PySimpleGUI as sg

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    ## adjust taxonomic level if neccessary
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header = 0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header = 0).fillna("nan")
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

    TaXon_table_n_samples = len(TaXon_table_samples)
    n_sites = len(set(Meta_data_table_df[meta_data_to_test].tolist()))

    answer = "Ask"
    output_message = "No"

    if (sorted(TaXon_table_samples) == sorted(Meta_data_table_samples) and TaXon_table_n_samples != n_sites):

        site_occupancy_dict = {}

        sites = set(Meta_data_table_df[meta_data_to_test].tolist())

        for site in sites:
            # this can either be a species name or the above specified taxonomic level
            present_OTU_list = []

            # extract samples that belong to the site from the metadata file
            included_samples_list = Meta_data_table_df[Meta_data_table_df.values  == site]['Samples'].values.tolist()

            # count the number of samples per site to calculate the site occupancy
            n_samples = len(included_samples_list)

            # create a list of all species (or the specified taxonomic level)
            if taxonomic_level == "OTUs":
                taxonomic_level = "ID"
            overall_included_species_list = TaXon_table_df[taxonomic_level].values.tolist()
            # make the list unique
            overall_included_species_set = set(overall_included_species_list)
            # remove potential 'nan's from the list
            overall_included_species_set = [x for x in overall_included_species_set if str(x) != 'nan']

            # create a set of species that is present at the sites
            for sample in included_samples_list:

                OTUs_per_species_list = []

                # check the read abundaces for each sample
                read_abundace_list = TaXon_table_df[sample].values.tolist()

                # enumerate the read abundaces for each sample and collect all lines that have more than one read
                for i, read_abundance in enumerate(read_abundace_list):
                    species = TaXon_table_df[taxonomic_level][i]
                    # if reads are present, collect the species name (or the specified taxonomic level) from the TaXon table
                    if read_abundance != 0:
                        OTUs_per_species_list.append(species)

                # remove all nans
                OTUs_per_species_list = [x for x in OTUs_per_species_list if str(x) != 'nan']
                # make list unique
                OTUs_per_species_list = list(set(OTUs_per_species_list))
                # append to list of species for the current site
                present_OTU_list.append(OTUs_per_species_list)

            # flatten the list of present species per site
            present_OTU_list_flattened = [val for sublist in present_OTU_list for val in sublist]

            # store occupancy of each species in a dict, will be accessed by position in list
            occupancy_dict = {}

            # count the number of occurences for each species and calculate the occpancy based on the number of samples
            for species in overall_included_species_set:
                count = present_OTU_list_flattened.count(species)
                occupancy = count / n_samples * 100
                occupancy_dict[species] = occupancy

            occupancy_dict = {k: v for k, v in sorted(occupancy_dict.items(), key=lambda item: item[1])}
            occupancy_list = list(occupancy_dict.values())
            species_list = list(occupancy_dict.keys())

            if (taxonomic_level == "Species" or taxonomic_level == "Genus"):
                x_values = ["<i>" + taxon + "</i>" for taxon in species_list]
            else:
                x_values = species_list

            occupancy_plot_directory = Path(str(path_to_outdirs) + "/" + "Site_occupancy_plots" + "/" + TaXon_table_xlsx.stem)
            if not os.path.exists(occupancy_plot_directory):
                os.mkdir(occupancy_plot_directory)

            fig = go.Figure(data=[go.Bar(x=x_values, y=occupancy_list)])
            fig.update_traces(marker_color=color1, marker_line_color=color2,marker_line_width=0.6, opacity=opacity_value)
            fig.update_layout(title_text=site + " (" + taxonomic_level + ")", yaxis_title="occupancy (%)")
            fig.update_layout(height=int(y_site_occ), width=int(x_site_occ), template=template, font_size=font_size, title_font_size=font_size)
            fig.update_yaxes(range=[0,100])
            fig.update_xaxes(tickmode='linear')
            fig.update_xaxes(tickangle=-90)


            output_pdf = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".pdf")
            output_html = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".html")
            occupancy_table = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".xlsx")
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))
            occupancy_df = pd.DataFrame(occupancy_list, species_list)
            occupancy_df.columns = ["Occupancy"]
            occupancy_df.index.name = "Taxon"
            occupancy_df = occupancy_df.sort_values("Occupancy")
            # sort the table numerical if OTUs were chosen
            if taxonomic_level == "ID":
                sort_list = []
                for OTU in occupancy_df.index.tolist():
                    sort_list.append(int(OTU.split("_")[1]))
                occupancy_df["sort"] = sort_list
                occupancy_df = occupancy_df.sort_values("sort")
                occupancy_df = occupancy_df.drop("sort", axis=1)
            occupancy_df.to_excel(occupancy_table)

        ## ask to show file
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## print closing text
        closing_text = "Site occupancy plots are found under:\n" + '/'.join(str(output_pdf).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        placeholder = TaXon_table_xlsx.name + " (multiple site occupancy plots)"
        ttt_log("site occupancy", "analysis", TaXon_table_xlsx.name, placeholder, meta_data_to_test, path_to_outdirs)

    else:
        sg.PopupError("Please check your Metadata file and Taxon table file: The samples do not match or the metadata is unique for all samples!", keep_on_top=True)

def site_occupancy_heatmap(TaXon_table_xlsx, path_to_outdirs, template, height, width, meta_data_to_test, taxonomic_level, font_size, color_discrete_sequence, add_categories_sum):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from pathlib import Path
    import webbrowser, os

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0).fillna("unidentified")
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    ## drop samples with metadata called nan (= empty)
    drop_samples = [i[0] for i in Meta_data_table_df.values.tolist() if i[1] == "nan"]

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

    if len(set(metadata_list)) == 1:
        sg.PopupError("Please choose more than one meta data category.")
    else:

        if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

            ## define variables
            samples = TaXon_table_samples
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
            for sample, category in zip(Meta_data_table_samples, metadata_list):
                if category not in category_dict.keys():
                    category_dict[category] = [sample]
                else:
                    category_dict[category] = category_dict[category] + [sample]

            ## collect all available taxa
            taxa = TaXon_table_df[taxonomic_level].values.tolist()

            ## check if the respective species are present in the collections
            taxon_presence_dict = {}
            n_rows, row_heights = [], []

            color_discrete_sequence = color_discrete_sequence * len(category_dict.keys())

            if (taxonomic_level == "Species" or taxonomic_level == "Genus"):
                x_values = ["<i>" + taxon + "</i>" for taxon in taxa]
            else:
                x_values = taxa

            if add_categories_sum == True:
                for samples in category_dict.values():
                    row_heights.append(len(samples))
                row_heights.append(len(set(metadata_list)))
                fig = make_subplots(rows=len(set(metadata_list)) + 1, cols=1, shared_xaxes=True, vertical_spacing=0.05, row_heights=row_heights)
            else:
                for samples in category_dict.values():
                    row_heights.append(len(samples))
                fig = make_subplots(rows=len(set(metadata_list)), cols=1, shared_xaxes=True, vertical_spacing=0.05, row_heights=row_heights)

            row = 1
            for metadata, samples in category_dict.items():
                if type(samples) == "str":
                    samples = [samples]
                z_values = []
                for sample in samples:
                    reads = TaXon_table_df[sample].values.tolist()
                    z_values = z_values + [[1 if x > 0 else 0 for x in reads]]
                y_values = samples
                fig.add_trace(go.Heatmap(z=z_values, x=x_values, y=y_values, showscale=False, xgap=1, ygap=1, hoverongaps = False, colorscale=[[0, "White"], [1, color_discrete_sequence[row-1]]]), row=row, col=1)
                row += 1

            if add_categories_sum == True:
                z_values, y_values = [], []
                for metadata, samples in category_dict.items():
                    reads = [sum(reads) for reads in TaXon_table_df[samples].values.tolist()]
                    z_values = z_values + [[1 if x > 0 else 0 for x in reads]]
                    y_values.append(metadata)
                fig.add_trace(go.Heatmap(z=z_values[::-1], x=x_values, y=y_values[::-1], showscale=False, xgap=1, ygap=1, hoverongaps = False, colorscale=[[0, "White"], [1, "Grey"]]), row=row, col=1)
                row += 1

            fig.update_layout(width=int(width), height=int(height), template="seaborn", font_size=font_size, yaxis_nticks=5, title_font_size=font_size)
            fig.update_xaxes(tickmode='linear')
            fig.update_yaxes(tickmode='linear')
            fig.update_xaxes(tickangle=-90)

            occupancy_plot_directory = Path(str(path_to_outdirs) + "/" + "Site_occupancy_plots" + "/" + TaXon_table_xlsx.stem)
            if not os.path.exists(occupancy_plot_directory):
                os.mkdir(occupancy_plot_directory)

            ## define output files
            output_pdf = Path(str(occupancy_plot_directory) + "/" + taxonomic_level + "_" + meta_data_to_test + "_heatmap.pdf")
            output_html = Path(str(occupancy_plot_directory) + "/" + taxonomic_level + "_" + meta_data_to_test + "_heatmap.html")

            ## write output files
            fig.write_image(str(output_pdf))
            fig.write_html(str(output_html))

            ## ask to show file
            answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
            if answer == "Yes":
                webbrowser.open('file://' + str(output_html))

            ## print closing text
            closing_text = "Site occupancy heatmaps are found under:\n" + '/'.join(str(output_pdf).split("/")[-4:])
            sg.Popup(closing_text, title="Finished", keep_on_top=True)

            ## write to log
            from taxontabletools.create_log import ttt_log
            placeholder = TaXon_table_xlsx.name + " (multiple site occupancy plots)"
            ttt_log("site occupancy", "analysis", TaXon_table_xlsx.name, "", meta_data_to_test, path_to_outdirs)


        else:
            sg.Popup("The metdata table and taXon table are not matching!")
