import random
import PySimpleGUI as sg
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import math, webbrowser
import plotly.express as px
from statistics import mean
from taxontabletools.taxontable_manipulation import strip_metadata
from plotly.subplots import make_subplots

def rarefaction_curve_legacy(TaXon_table_xlsx, repetitions, path_to_outdirs, template, theme, font_size, taxonomic_level_1):

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    ## create a y axis title text
    taxon_title = taxonomic_level_1.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level_1 in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level_1
        taxonomic_level_1 = "ID"

    ## load TaxonTable
    TaXon_table_file = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx).fillna('nan')
    df = strip_metadata(df)

    available_samples = df.columns.tolist()[10:]
    sample_dict_clean = {}

    # iterate through all available samples
    for sample in available_samples:
        # create a dict for the read numbers of the respective sample for each species
        sample_OTU_list = df[[sample, taxonomic_level_1]].values.tolist()
        # select only the present Species
        sample_species_list = list(set([OTU[1] for OTU in sample_OTU_list if (OTU[0] != 0 and OTU[1] != "nan")]))
        # store the species in a dictionary
        sample_dict_clean[sample] = sample_species_list

    # draw once for each sample
    number_of_draws = len(sample_dict_clean.keys())

    # dictionary to store the drawing results
    draw_dictionary = {}

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / repetitions
    ############################################################################

    for n_reps in range(0, repetitions):
        # store the original dictionary to start over again
        # a copy of the original dictionary is required, because the samples will be removed with each draw
        # thus for each replicate a new dictionary to draw from has to be created
        sample_dict_to_draw = dict(sample_dict_clean)

        species_list = []
        species_set = []

        for i in range(0, number_of_draws):
            # choose a random sample from the dictionary
            random_choice = random.choice(list(sample_dict_to_draw.keys()))
            # extract the OTU IDs from the chosen sample and add them to the already existing OTU IDs
            species_list = species_list + sample_dict_clean[random_choice]
            # create a unique set
            species_set = set(species_list)
            # number of OTUs
            n_species = len(species_set)
            # now add the unique OTU list to the output dictionary
            # if the key is not in the dict, create a new entry (= OTU ID plus number of OTUs)
            if i not in draw_dictionary.keys():
                draw_dictionary[i] = [n_species]
            # if the key already exists, calculate the sum of the already existing number of OTUs and the new number of OTUs
            else:
                # create a new list to store the current number of OTUs
                add_species_list = draw_dictionary[i]
                add_species_list.append(n_species)
                draw_dictionary[i] = add_species_list

            # remove the sample to draw only once
            sample_dict_to_draw.pop(random_choice)

        ############################################################################
        event, values = window_progress_bar.read(timeout=1)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    # create a dict to store the average number of OTUs per draw
    rarefaction_dict_average, rarefaction_dict_stdef = {}, {}

    def average(lst):
        return sum(lst) / len(lst)

    # iterate through the draw_dictionary and calculate the average number of OTUs
    for key, value in draw_dictionary.items():
        average_species = average(draw_dictionary[key])
        stdef_species = np.std(draw_dictionary[key], dtype=np.float64)
        rarefaction_dict_average[key] = average_species
        rarefaction_dict_stdef[key] = stdef_species

    # draw the plot
    draws = [i+1 for i in rarefaction_dict_average.keys()]
    n_species = list(rarefaction_dict_average.values())
    error_bar = list(rarefaction_dict_stdef.values())
    y_axis_title = "# " + taxon_title
    fig = go.Figure(data=[go.Scatter(x=draws, y=n_species, error_y=dict(type='data', array=error_bar, thickness=0.5, width=3, visible=True))])
    fig.update_layout(title_text="repetitions = " + str(n_reps+1), yaxis_title=y_axis_title, xaxis_title="# samples")
    fig.update_traces(marker_color=color1, marker_line_color=color2, opacity=opacity_value)
    fig.update_layout(height=800, width=1200, template=template, showlegend=False, font_size=font_size, title_font_size=font_size)

    ## write files
    output_pdf = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_" + taxon_title + ".pdf")
    output_html = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_taxa" + taxon_title + ".html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Rarefaction curves are found in: " + str(path_to_outdirs) + "/rarefaction_curves/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write log
    from taxontabletools.create_log import ttt_log
    ttt_log("rarefaction curve all-in-one", "analysis", TaXon_table_file.name, output_pdf.name, "nan", path_to_outdirs)

def rarefaction_curve_taxa(TaXon_table_xlsx, repetitions, path_to_outdirs, template, font_size, taxonomic_level_1, taxonomic_level_2, color_discrete_sequence):

    ## load TaxonTable
    TaXon_table_file = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx).fillna('unidentified')
    df = strip_metadata(df)

    ## create a y axis title text
    taxon_title = taxonomic_level_1.lower()

    ## adjust taxonomic level if neccessary
    if taxonomic_level_1 in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level_1
        taxonomic_level_1 = "ID"

    ## collect available samples
    available_samples = df.columns.tolist()[10:]
    sample_dict_clean = {}

    ## collect available taxa
    available_taxa = [taxon for taxon in set(df[taxonomic_level_2].values.tolist()) if taxon != "nan"]

    ## create a color dict
    ## extend to the color dict
    color_discrete_sequence = color_discrete_sequence * len(available_taxa)
    color_dict = {}
    for i, taxon in enumerate(available_taxa):
        color_dict[taxon] = color_discrete_sequence[i]

    ## collect the increase for each taxon
    increase_dict = {}

    ## create an empty figure
    fig = go.Figure()

    for taxon in sorted(available_taxa):
        df_filtered = df.loc[df[taxonomic_level_2] == taxon]

        # iterate through all available samples
        for sample in available_samples:
            # create a dict for the read numbers of the respective sample for each species
            sample_OTU_list = df_filtered[[sample, taxonomic_level_1]].values.tolist()
            # select only the present Species
            sample_species_list = list(set([OTU[1] for OTU in sample_OTU_list if (OTU[0] != 0 and OTU[1] != "nan")]))
            # store the species in a dictionary
            sample_dict_clean[sample] = sample_species_list

        # draw once for each sample
        number_of_draws = len(sample_dict_clean.keys())

        # dictionary to store the drawing results
        draw_dictionary = {}

        for n_reps in range(0, repetitions):
            # store the original dictionary to start over again
            # a copy of the original dictionary is required, because the samples will be removed with each draw
            # thus for each replicate a new dictionary to draw from has to be created
            sample_dict_to_draw = dict(sample_dict_clean)

            species_list = []
            species_set = []

            for i in range(0, number_of_draws):
                # choose a random sample from the dictionary
                random_choice = random.choice(list(sample_dict_to_draw.keys()))
                # extract the OTU IDs from the chosen sample and add them to the already existing OTU IDs
                species_list = species_list + sample_dict_clean[random_choice]
                # create a unique set
                species_set = set(species_list)
                # number of OTUs
                n_species = len(species_set)
                # now add the unique OTU list to the output dictionary
                # if the key is not in the dict, create a new entry (= OTU ID plus number of OTUs)
                if i not in draw_dictionary.keys():
                    draw_dictionary[i] = [n_species]
                # if the key already exists, calculate the sum of the already existing number of OTUs and the new number of OTUs
                else:
                    # create a new list to store the current number of OTUs
                    add_species_list = draw_dictionary[i]
                    add_species_list.append(n_species)
                    draw_dictionary[i] = add_species_list

                # remove the sample to draw only once
                sample_dict_to_draw.pop(random_choice)

        # create a dict to store the average number of OTUs per draw
        rarefaction_dict_average, rarefaction_dict_stdef = {}, {}

        def average(lst):
            return sum(lst) / len(lst)

        # iterate through the draw_dictionary and calculate the average number of OTUs
        for key, value in draw_dictionary.items():
            average_species = average(draw_dictionary[key])
            stdef_species = np.std(draw_dictionary[key], dtype=np.float64)
            rarefaction_dict_average[key] = average_species
            rarefaction_dict_stdef[key] = stdef_species

        ## add to plot
        draws = [i+1 for i in rarefaction_dict_average.keys()]
        n_species = list(rarefaction_dict_average.values())
        increase_dict[taxon] = n_species
        error_bar = list(rarefaction_dict_stdef.values())
        fig.add_trace(go.Scatter(x=draws, y=n_species, name=taxon, marker_color=color_dict[taxon], error_y=dict(type='data', array=error_bar, thickness=0.5, width=3, visible=True)))

    ## update figure
    y_axis_title = "# " + taxon_title
    fig.update_layout(title_text="repetitions = " + str(n_reps+1), yaxis_title=y_axis_title, xaxis_title="# samples")
    fig.update_layout(height=700, width=1200, template="simple_white", showlegend=True, font_size=font_size, title_font_size=font_size)
    fig.update_xaxes(rangemode="tozero")
    fig.update_yaxes(rangemode="tozero")
    fig.update_layout(height=800, width=1200, template=template, showlegend=True, font_size=font_size, title_font_size=font_size)

    ## write files
    out_name = taxonomic_level_1.lower() + "_" + taxonomic_level_2.lower()
    output_pdf = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_" + out_name + ".pdf")
    output_html = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_" + out_name + ".html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Rarefaction curves are found in: " + str(path_to_outdirs) + "/rarefaction_curves/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write log
    from taxontabletools.create_log import ttt_log
    ttt_log("rarefaction curve per taxon", "analysis", TaXon_table_file.name, output_pdf.name, "nan", path_to_outdirs)

def rarefaction_curve_reads(TaXon_table_xlsx, repetitions, width, height, path_to_outdirs, template, theme, font_size):

    ## load TaxonTable
    TaXon_table_file = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df = strip_metadata(TaXon_table_df)
    samples = TaXon_table_df.columns.tolist()[10:]
    scatter_size = 5

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    height = int(height)

    ## count rows and columns to create subplots
    n_rows = math.ceil(len(samples) / 4)
    n_columns = 5
    column_count = 1
    row_count = 1
    fig = make_subplots(rows=n_rows, cols=4, subplot_titles=samples, shared_yaxes=True)

    ## calculate maximum number of OTUs
    max_OTUs = []
    for sample in samples:
        max_OTUs.append(len([OTU for OTU in TaXon_table_df[sample] if OTU != 0]))
    y_limit = max(max_OTUs) + 20

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples)
    ############################################################################

    ############################################################################
    event, values = window_progress_bar.read(timeout=1)
    if event == 'Cancel'  or event is None:
        window_progress_bar.Close()
        raise RuntimeError
    # update bar with loop value +1 so that bar eventually reaches the maximum
    progress_update += 0
    progress_bar.UpdateBar(progress_update)
    ############################################################################

    for sample in samples:

        ## filter sample from data
        read_df = TaXon_table_df[[sample, "ID"]]
        ## drop empty OTUs
        read_df = read_df[read_df[sample] != 0]
        ## create read list to draw the subsamples from
        read_list = pd.Series(np.repeat(read_df['ID'].to_list(), read_df[sample].to_list()))

        output = []

        ## draw random sample
        for perc in np.arange(0.00, 1.05, 0.05):
            ## calculate sample size
            sub_sample_size = int(len(read_list) * perc)

            ## draw X subsamples of that size
            mean_species = mean([read_list.sample(n = sub_sample_size).nunique() for i in range(repetitions)])

            output.append(mean_species)

        output = pd.DataFrame({'percentage': np.arange(0.00, 1.05, 0.05), 'mean_OTUs': output})

        ## write plot
        fig.add_trace(go.Scatter(x=output["percentage"], y=output["mean_OTUs"], name=sample, mode='markers+lines', marker=dict(size=int(scatter_size))), row=row_count, col=column_count)
        fig.update_traces(marker_color=color1, marker_line_color=color2, opacity=opacity_value, row=row_count, col=column_count)
        fig.update_yaxes(range=[0, y_limit], row=row_count, col=column_count)

        ## add a y axis title to all left bound plots
        if column_count == 1:
            fig.update_yaxes(title_text="# OTUs", row=row_count, col=column_count)

        ## add x axis title to all plots in the last row
        # if row_count == n_rows:
        #     fig.update_xaxes(title_text="subsample (%)", row=row_count, col=column_count)

        column_count += 1
        if column_count == n_columns:
            column_count = 1
            row_count += 1
            height += 100

        ############################################################################
        event, values = window_progress_bar.read(timeout=1)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    fig.update_layout(height=int(height), width=int(width), template=template, font_size=font_size, title_font_size=font_size, showlegend=False)
    fig.update_yaxes(rangemode="tozero")
    fig.update_xaxes(rangemode="tozero")

    ## write files
    output_pdf = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_reads.pdf")
    output_html = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_reads.html")
    fig.write_image(str(output_pdf))
    fig.write_html(str(output_html))

    ## ask to show file
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(output_html))

    ## print closing text
    closing_text = "Rarefaction curves are found in: " + str(path_to_outdirs) + "/rarefaction_curves/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    ## write log
    from taxontabletools.create_log import ttt_log
    ttt_log("rarefaction curve reads", "analysis", TaXon_table_file.name, output_pdf.name, repetitions, path_to_outdirs)
