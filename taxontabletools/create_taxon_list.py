def create_taxon_list(TaXon_table_xlsx, taxon_list_output_file_name, create_gbif_link, calc_dist, calc_occupancy, taxon_tools_version, path_to_outdirs):

    import requests_html, json
    import re
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    from datetime import date
    import numpy as np
    from pathlib import Path
    from Bio import pairwise2
    from Bio.Seq import Seq

    # load the TaXon table
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    TaXon_table_df = TaXon_table_df.fillna("nan")
    samples = TaXon_table_df.columns.tolist()[10:]
    n_samples = len(samples)
    total_reads = sum([sum(reads) for reads in TaXon_table_df[samples].values.tolist()])

    # create the output files
    if taxon_list_output_file_name == '':
        taxon_list_output_file_name = "my_taxon_list"
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".xlsx")

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Read counts and OTU numbers')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar1')],
              [sg.Text('Intraspecific distances')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar2')],
              [sg.Text('Occupancy')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar3')],
              [sg.Text('GBIF link')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar4')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar1 = window_progress_bar['progressbar1']
    progress_bar2 = window_progress_bar['progressbar2']
    progress_bar3 = window_progress_bar['progressbar3']
    progress_bar4 = window_progress_bar['progressbar4']
    ############################################################################

    # progress bar 1
    progress_update = 0
    progress_increase = 1000 / len(TaXon_table_df.values.tolist()) + 1

    # create a dict with read numbers from the TaXon table
    overall_included_OTU_set = []
    taxa_list_dict = {}
    for OTU in TaXon_table_df.values.tolist():
        species = ' '.join(OTU[1:7])
        taxonomy = OTU[1:7]
        overall_included_OTU_set.append(' '.join(taxonomy))
        n_reads = sum(OTU[10:])
        n_reads_rel = n_reads / total_reads * 100
        n_OTUs = 1
        if species not in taxa_list_dict.keys():
            taxa_list_dict[species] = taxonomy + [n_reads, n_reads_rel, n_OTUs]
        else:
            n_reads = taxa_list_dict[species][6] + n_reads
            n_reads_rel = taxa_list_dict[species][7] + n_reads_rel
            n_OTUs = taxa_list_dict[species][8] + n_OTUs
            taxa_list_dict[species] = taxonomy + [n_reads, n_reads_rel, n_OTUs]

        ############################################################################
        event, values = window_progress_bar.read(timeout=5)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar1.UpdateBar(progress_update)
        ############################################################################

    ##############################################################################
    # calculate intraspecific distances for Species with multiple OTUs

    def pw_distance(seq1, seq2):
        len1 = len(seq1)
        len2 = len(seq2)
        length = max([len1, len2])
        alignments = pairwise2.align.globalxx(seq1, seq2)
        max_score = max([alignment.score for alignment in alignments])
        diff = 100 - (max_score / length * 100)
        diff = round(diff, 2)
        return(diff)

    OTU_sequences_dict = {}
    distances_dict = {}

    # check if distance are to be calulated
    if calc_dist == True:
        for OTU in TaXon_table_df.values.tolist():
            species = OTU[6]
            if species != "nan":
                sequence = OTU[9]
                if species not in OTU_sequences_dict.keys():
                    OTU_sequences_dict[species] = [sequence]
                else:
                    OTU_sequences_dict[species] = OTU_sequences_dict[species] + [sequence]

        # progress bar 2
        progress_update = 0
        progress_increase = 1000 / len(OTU_sequences_dict.items()) + 1

        for species, sequences in OTU_sequences_dict.items():
            distances = []
            if len(sequences) > 1:
                for main_sequence in sequences:
                    for sub_sequence in sequences:
                        distances.append(pw_distance(main_sequence, sub_sequence))
                distances = list(filter((0.0).__ne__, distances))
                max_dist = max(distances)
                min_dist = min(distances)
                avg_dist = sum(distances) / len(distances)
                avg_dist = round(avg_dist, 2)
                distances_dict[species] = [max_dist, min_dist, avg_dist]

            ############################################################################
            event, values = window_progress_bar.read(timeout=5)
            if event == 'Cancel'  or event is None:
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar2.UpdateBar(progress_update)
            ############################################################################

    ##############################################################################
    # calculate sample occupancy

    if calc_occupancy == True:

        # progress bar 1
        progress_update = 0
        progress_increase = 1000 / len(samples) + 1

        present_OTU_list = []
        for sample in samples:
            OTUs_per_species_list = []
            # check the read abundaces for each sample
            read_abundace_list = TaXon_table_df[sample].values.tolist()
            # enumerate the read abundaces for each sample and collect all lines that have more than one read
            for i, read_abundance in enumerate(read_abundace_list):
                OTU = TaXon_table_df["IDs"][i]
                taxonomy = ' '.join(TaXon_table_df[TaXon_table_df['IDs'].str.contains(OTU + "$", regex=True)].values.tolist()[0][1:7])
                # if reads are present, collect the species name (or the specified taxonomic level) from the TaXon table
                if read_abundance != 0:
                    OTUs_per_species_list.append(taxonomy)
            # remove all nans
            OTUs_per_species_list = [x for x in OTUs_per_species_list if str(x) != 'nan']
            # make list unique
            OTUs_per_species_list = list(set(OTUs_per_species_list))
            # append to list of species for the current site
            present_OTU_list.append(OTUs_per_species_list)

            ############################################################################
            event, values = window_progress_bar.read(timeout=5)
            if event == 'Cancel'  or event is None:
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar3.UpdateBar(progress_update)
            ############################################################################

        # store occupancy of each species in a dict, will be accessed by position in list
        occupancy_dict = {}
        # flatten the list of present species per site
        present_OTU_list_flattened = [val for sublist in present_OTU_list for val in sublist]
        # count the number of occurences for each species and calculate the occpancy based on the number of samples
        for taxon in overall_included_OTU_set:
            count = present_OTU_list_flattened.count(taxon)
            occupancy = round(count / n_samples * 100, 2)
            if len(taxa_list_dict[taxon]) == 9:
                taxa_list_dict[taxon] = taxa_list_dict[taxon] + [occupancy]
    else:
        for taxon in overall_included_OTU_set:
            if len(taxa_list_dict[taxon]) == 9:
                taxa_list_dict[taxon] = taxa_list_dict[taxon] + [""]

    ##############################################################################

    def gbif_requester(species_name):

        ## create an html session
        with requests_html.HTMLSession() as session:
            ## generate html request name
            request_name = '%20'.join(species_name.split(' '))

            ## request that name
            r = session.get('https://api.gbif.org/v1/species?name=%s&limit=1' % request_name)

            ## parse json
            res = json.loads(r.text)

            ## return the link
            if res['results']:
                return 'https://www.gbif.org/species/%s' % res['results'][0]['key']
            else:
                return False

    gbif_df_list = []

    # progress bar 3
    progress_update = 0
    progress_increase = 1000 / len(taxa_list_dict.items()) + 1

    if create_gbif_link == True:

        for taxon, values in taxa_list_dict.items():
            species = values[-5]
            if create_gbif_link == True:
                if species != 'nan':
                    gbif_link = gbif_requester(species)
                    if gbif_link != False:
                        species_entry = values + [gbif_link]
                        gbif_df_list.append(species_entry)
                    else:
                        species_entry = values + [""]
                        gbif_df_list.append(species_entry)
                else:
                    species_entry = values + [""]
                    gbif_df_list.append(species_entry)
            else:
                species_entry = values + [""]
                gbif_df_list.append(species_entry)

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar4.UpdateBar(progress_update)
            ############################################################################

    else:
        for taxon, values in taxa_list_dict.items():
            species_entry = values + [""]
            gbif_df_list.append(species_entry)

    window_progress_bar.Close()

    # now append the distances for each species that has multiple OTUs
    # create the final list for the output dataframes
    df_out_list = []
    for entry in gbif_df_list:
        species = entry[5]
        if species in distances_dict.keys():
            df_out_list.append(entry + distances_dict[species])
        else:
            df_out_list.append(entry + ["", "", ""])

    df_out  = pd.DataFrame(df_out_list, columns=TaXon_table_df.columns.tolist()[1:7] + ["reads", "reads (%)", "OTUs", "occupancy (%)", "gbif", "dist max (%)", "dist min (%)", "dist avg (%)"])
    df_out = df_out[TaXon_table_df.columns.tolist()[1:7] + ["reads", "reads (%)", "OTUs", "occupancy (%)", "dist max (%)", "dist min (%)", "dist avg (%)", "gbif"]]
    df_out.to_excel(output_xlsx, sheet_name = 'Taxa', index=False)

    closing_text = "Taxon list is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("taxon list", "analysis", TaXon_table_xlsx.name, output_xlsx.name, "", path_to_outdirs)
