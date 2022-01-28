import requests_html, json
import re, sys, subprocess, os
import PySimpleGUI as sg
import pandas as pd
from pandas import DataFrame
from datetime import date
import numpy as np
from pathlib import Path
from Bio import pairwise2
from Bio.Seq import Seq

def create_taxon_list(TaXon_table_xlsx, taxon_list_output_file_name, create_gbif_link, calc_dist, use_metadata, taxon_tools_version, path_to_outdirs, clustering_unit):

    ##### functions
    def open_table(table):
        if sys.platform == "win32":
            os.startfile(table)
        else:
            opener = "open" if sys.platform == 'darwin' else 'xdg-open'
            subprocess.call([opener, table])

    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    def collect_metadata(Meta_data_table_df):
        available_metadata = Meta_data_table_df.columns.tolist()[1:]
        win2_active = True
        meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(available_metadata)], 5))
        layout2 = [[sg.Text("Site occupancy", size=(20,1))],
        [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
        [sg.Text('',size=(1,1))],
        [sg.Button('Calculate')],
        [sg.Button('Skip')]]


        win2 = sg.Window('Taxon list', layout2, keep_on_top=True)

        while True:
            event2, values2 = win2.Read()
            if event2 == 'Calculate':
                for input_value, test in values2.items():
                    if (test == True):
                        meta_data_to_test = input_value
                win2.Close()
                win2_active = False
                return meta_data_to_test

            if event2 is None or event2 == "Skip":
                win2.Close()
                win2_active = False
                return False

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
                species_id = res['results'][0]['key']
                link = ['https://www.gbif.org/species/%s' % res['results'][0]['key']]
            else:
                species_id = False
                link = [""]

            if species_id != False:
                ## collect vernacular name
                ## create an html session
                with requests_html.HTMLSession() as session:
                    ## request that name
                    r = session.get('https://api.gbif.org/v1/species/%s/vernacularNames' % species_id)

                    try:
                        ## parse json
                        res = json.loads(r.text)
                        ## search for english common name
                        name = list(set([name["vernacularName"] for name in res['results'] if name["language"] == "eng"]))
                        if len(name) > 1:
                            ## remove duplicates from the vernacular name list
                            name = list(set([i.lower() for i in name]))
                            name = [', '.join(name)]
                        if name == []:
                            name = [""]
                        return name + link
                    except:
                        return [""] + link
            else:
                return ["", ""]

    def pw_distance(seq1, seq2):
        len1 = len(seq1)
        len2 = len(seq2)
        length = max([len1, len2])
        alignments = pairwise2.align.globalxx(seq1, seq2)
        max_score = max([alignment.score for alignment in alignments])
        diff = 100 - (max_score / length * 100)
        diff = round(diff, 2)
        return(diff)

    # load the TaXon table
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    TaXon_table_df = TaXon_table_df.fillna("nan")
    samples = TaXon_table_df.columns.tolist()[10:]
    n_samples = len(samples)
    total_reads = sum([sum(reads) for reads in TaXon_table_df[samples].values.tolist()])

    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    if use_metadata == True:
        if Meta_data_table_xlsx.exists():
            ## load metadata list
            Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0).fillna("nan")
            Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
            meta_data_to_test = collect_metadata(Meta_data_table_df)

            if meta_data_to_test != False:

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

                metadata_samples_dict = {}

                for i in Meta_data_table_df[["Samples", meta_data_to_test]].values.tolist():
                    metadata = i[1]
                    sample = i[0]
                    if metadata in metadata_samples_dict.keys():
                        metadata_samples_dict[metadata] = metadata_samples_dict[metadata] + [sample]
                    else:
                        metadata_samples_dict[metadata] = [sample]
                for metadata, samples in metadata_samples_dict.items():
                    TaXon_table_df[metadata] = [sum(i) for i in TaXon_table_df[samples].values.tolist()]
                    TaXon_table_df = TaXon_table_df.drop(samples, axis=1)
                samples = TaXon_table_df.columns.tolist()[10:]
                n_samples = len(samples)
                total_reads = sum([sum(reads) for reads in TaXon_table_df[samples].values.tolist()])
        else:
            sg.PopupOK("Warning: Metadata table does not exists.\n\nPlease first create a metadata table!\n\nSkipping metadata occupancy.", title="Warning")

    taxonomic_levels = TaXon_table_df.columns.values.tolist()[1:7]
    TaXon_table_df["FullName"] = ['<>'.join(t) for t in TaXon_table_df[taxonomic_levels].values.tolist()]

    ## create an empty data frame
    taxon_list_df = pd.DataFrame()
    ## add unique species IDs
    unique_species = list(set(['<>'.join(t) for t in TaXon_table_df[taxonomic_levels].values.tolist()]))
    taxon_list_df["FullName"] = unique_species

    # create the output files
    if taxon_list_output_file_name == '':
        taxon_list_output_file_name = "my_taxon_list"
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".xlsx")

    answer = "Yes"
    if output_xlsx.exists():
        answer = sg.PopupYesNo("Taxon list already exists! Overwrite?")

    if answer == "Yes":
        ############################################################################
        ## create the progress bar window
        layout = [[sg.Text('Calculating statistics for taxon list.')],
                  [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
                  [sg.Cancel()]]
        window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
        progress_bar = window_progress_bar['progressbar']
        ############################################################################

        # progress bar 1
        progress_update = 0
        progress_increase = 1000 / len(taxon_list_df.values.tolist())

        ############################################################################
        # 1

        n_OTUs_list = []
        n_reads_list = []
        n_reads_rel_list = []
        max_dist_list = []
        min_dist_list = []
        avg_dist_list = []
        rel_occurrence_list = []
        abs_occurrences_list = []
        vernacular_name_list = []
        gbif_link_list = []

        for species in taxon_list_df["FullName"]:
            n_OTUs = len([i for i in TaXon_table_df.values.tolist() if species in i])
            n_reads = sum([sum(i[10:-1]) for i in TaXon_table_df.values.tolist() if species in i])
            n_reads_rel = round(n_reads / total_reads * 100, 3)

            if use_metadata == True:
                df = pd.DataFrame([i[10:-1] for i in TaXon_table_df.values.tolist() if species in i], columns=samples)
                sample_occurrences = []
                for sample in samples:
                    if sum(df[sample]) != 0:
                        sample_occurrences.append(sample)
                rel_occurences = round(len(sample_occurrences) / n_samples * 100, 3)

                abs_occurrences = []
                for sample in samples:
                    if sample in sample_occurrences:
                        abs_occurrences.append('x')
                    else:
                        abs_occurrences.append('')
            else:
                df = pd.DataFrame([i[10:-1] for i in TaXon_table_df.values.tolist() if species in i], columns=samples)
                sample_occurrences = []
                for sample in samples:
                    if sum(df[sample]) != 0:
                        sample_occurrences.append(sample)
                rel_occurences = round(len(sample_occurrences) / n_samples * 100, 3)

                abs_occurrences = []
                for sample in samples:
                    if sample in sample_occurrences:
                        abs_occurrences.append('x')
                    else:
                        abs_occurrences.append('')

            n_OTUs_list.append(n_OTUs)
            n_reads_list.append(n_reads)
            n_reads_rel_list.append(n_reads_rel)
            rel_occurrence_list.append(rel_occurences)
            abs_occurrences_list.append(abs_occurrences)

            if calc_dist == True:
                species_name = species.split("<>")[-1]
                if n_OTUs != 1 and species_name != "nan":
                    distances = []
                    sequences = [i[9] for i in TaXon_table_df.values.tolist() if species in i]
                    for main_sequence in sequences:
                        for sub_sequence in sequences:
                            distances.append(pw_distance(main_sequence, sub_sequence))
                    distances = list(filter((0.0).__ne__, distances))
                    max_dist = max(distances)
                    min_dist = min(distances)
                    avg_dist = sum(distances) / len(distances)
                    avg_dist = round(avg_dist, 2)
                    max_dist_list.append(max_dist)
                    min_dist_list.append(min_dist)
                    avg_dist_list.append(avg_dist)
                else:
                    max_dist_list.append('')
                    min_dist_list.append('')
                    avg_dist_list.append('')

            if create_gbif_link == True:
                species_name = species.split("<>")[-1]
                if species_name != "nan":
                    res = gbif_requester(species_name)
                    vernacular_name_list.append(res[0])
                    gbif_link_list.append(res[1])
                else:
                    vernacular_name_list.append('')
                    gbif_link_list.append('')

            ############################################################################
            event, values = window_progress_bar.read(timeout=5)
            if event == 'Cancel'  or event is None:
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        taxon_list_df["Phylum"] = [i.split("<>")[0] if i.split("<>")[0] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df["Class"] = [i.split("<>")[1] if i.split("<>")[1] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df["Order"] = [i.split("<>")[2] if i.split("<>")[2] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df["Family"] = [i.split("<>")[3] if i.split("<>")[3] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df["Genus"] = [i.split("<>")[4] if i.split("<>")[4] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df["Species"] = [i.split("<>")[5] if i.split("<>")[5] != "nan" else "" for i in taxon_list_df["FullName"].values.tolist()]
        taxon_list_df.drop('FullName', axis='columns', inplace=True)
        taxon_list_df["OTUs"] = n_OTUs_list
        taxon_list_df["Reads"] = n_reads_list
        taxon_list_df["Reads (%)"] = n_reads_rel_list
        taxon_list_df["Occurrence (%)"] = rel_occurrence_list
        df_occurrences = pd.DataFrame(abs_occurrences_list, columns=samples)
        taxon_list_df = taxon_list_df.join(df_occurrences)
        if calc_dist == True:
            taxon_list_df["Max. dist."] = max_dist_list
            taxon_list_df["Min. dist."] = min_dist_list
            taxon_list_df["Avg. dist."] = avg_dist_list
        if create_gbif_link == True:
            taxon_list_df["Vernacular name"] = vernacular_name_list
            taxon_list_df["GBIF link"] = gbif_link_list


        ##############################################################################

        taxon_list_df.to_excel(output_xlsx, sheet_name = 'Taxa', index=False)

        answer = sg.PopupYesNo("Open taxon list?")
        if answer == "Yes":
            open_table(output_xlsx)

        closing_text = "Taxon list is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        from taxontabletools.create_log import ttt_log
        ttt_log("taxon list", "analysis", TaXon_table_xlsx.name, output_xlsx.name, "", path_to_outdirs)
