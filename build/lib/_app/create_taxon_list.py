def create_taxon_list(TaXon_table_xlsx, taxon_list_output_file_name, language, values, create_gbif_link, taxon_tools_version, path_to_outdirs):

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
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    TaXon_table_df = TaXon_table_df.fillna("nan")
    samples = TaXon_table_df.columns.tolist()[10:]
    n_samples = len(samples)

    # create the output files
    if taxon_list_output_file_name == '':
        taxon_list_output_file_name = "my_taxon_list"
    output_txt = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".txt")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".xlsx")

    # write the text file
    # open document
    document = open(output_txt, "w")
    # add explanation to beginning of the file
    if language == "German":

        # if german: add explanation to taxon list
        document.write('Taxonliste ' + "(" + str(date.today()) + ")" + "\n")
        document.write('\n')
        document.write('Datenblatt:' + "\n")
        document.write('Beschreibung:: ' + values["TL_description"] + "\n")
        document.write('Author(en): ' + values["TL_authors"] + "\n")
        document.write('Laborprotokoll: ' + values["TL_lab_protocol"] + "\n")
        document.write('Anz. Replikate: ' + values["TL_replicates"] + "\n")
        document.write('Anz. Negativ-Kontrollen: ' + values["TL_negative_controls"] + "\n")
        document.write('Primer: ' + values["TL_primers"] + "\n")
        document.write('Sequenzierlauf: ' + values["TL_sequencing_run"] + "\n")
        document.write('Bionf. Pipeline: ' + values["TL_pipeline"] + "\n")
        document.write('Taxon Table Tools: ' + taxon_tools_version + "\n")

        document.write('\n')
        document.write('Taxonliste: ')
        document.write('Folgende Taxonliste enthält die durch DNA-Metabarcoding ermittelten Taxa. Der GBIF Link befindet sich in der letzen Spalte.' + "\n")
        document.write('\n')
        document.write('Excel-Tabelle: ' + str(output_xlsx) + "\n")

    if language == "English":

        document.write('Taxon list ' + "(" + str(date.today()) + ")" + "\n")
        document.write('\n')
        document.write('Data sheet:' + "\n")
        document.write('Description: ' + values["TL_description"] + "\n")
        document.write('Author(s): ' + values["TL_authors"] + "\n")
        document.write('Lab protocol: ' + values["TL_lab_protocol"] + "\n")
        document.write('No. of replicates: ' + values["TL_replicates"] + "\n")
        document.write('No. of negative controls: ' + values["TL_negative_controls"] + "\n")
        document.write('Primers: ' + values["TL_primers"] + "\n")
        document.write('Sequencing run: ' + values["TL_sequencing_run"] + "\n")
        document.write('Bionf. pipeline: ' + values["TL_pipeline"] + "\n")
        document.write('Taxon Table Tools: ' + taxon_tools_version + "\n")
        document.write('\n')
        document.write('Taxon list: ')
        document.write('The following taxon list contains the via DNA metabarcoding identified taxa. The GBIF link is given in the last column.' + "\n")
        document.write('\n')
        document.write('Excel-Sheet: ' + str(output_xlsx))

    # write the taxa list
    ##############################################################################
    # create a dict with read numbers from the TaXon table
    overall_included_OTU_set = []
    taxa_list_dict = {}
    for OTU in TaXon_table_df.values.tolist():
        species = ' '.join(OTU[1:7])
        taxonomy = OTU[1:7]
        overall_included_OTU_set.append(' '.join(taxonomy))
        n_reads = sum(OTU[10:])
        n_OTUs = 1
        if species not in taxa_list_dict.keys():
            taxa_list_dict[species] = taxonomy + [n_reads, n_OTUs]
        else:
            n_reads = taxa_list_dict[species][6] + n_reads
            n_OTUs = taxa_list_dict[species][7] + n_OTUs
            taxa_list_dict[species] = taxonomy + [n_reads, n_OTUs]
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

    for OTU in TaXon_table_df.values.tolist():
        species = OTU[6]
        if species != "nan":
            sequence = OTU[9]
            if species not in OTU_sequences_dict.keys():
                OTU_sequences_dict[species] = [sequence]
            else:
                OTU_sequences_dict[species] = OTU_sequences_dict[species] + [sequence]

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

    ##############################################################################
    # calculate sample occupancy
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
    # store occupancy of each species in a dict, will be accessed by position in list
    occupancy_dict = {}
    # flatten the list of present species per site
    present_OTU_list_flattened = [val for sublist in present_OTU_list for val in sublist]
    # count the number of occurences for each species and calculate the occpancy based on the number of samples
    for taxon in overall_included_OTU_set:
        count = present_OTU_list_flattened.count(taxon)
        occupancy = round(count / n_samples, 2)
        if len(taxa_list_dict[taxon]) == 8:
            taxa_list_dict[taxon] = taxa_list_dict[taxon] + [occupancy]
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

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(taxa_list_dict.keys()) + 1
    ############################################################################

    gbif_df_list = []

    for taxon, values in taxa_list_dict.items():
        species = values[-4]
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
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

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

    df_out  = pd.DataFrame(df_out_list, columns=TaXon_table_df.columns.tolist()[1:7] + ["reads", "OTUs", "occupancy", "gbif", "dist max", "dist min", "dist avg"])
    df_out = df_out[TaXon_table_df.columns.tolist()[1:7] + ["reads", "OTUs", "occupancy", "dist max", "dist min", "dist avg", "gbif"]]
    df_out.to_excel(output_xlsx, sheet_name = 'Taxa', index=False)
    document.close()

    closing_text = "Taxon list is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
    print(closing_text)
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
