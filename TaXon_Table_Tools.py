import csv, glob, sys, subprocess, os, random
import PySimpleGUI as sg
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import plot, ion, show
from datetime import date
from io import StringIO
from pathlib import Path
from matplotlib.pyplot import cm
import webbrowser

##########################################################################################################################
# required python version: 3.6.8

# update version here (will be displayed on the main layout)
# Support for: u = ubuntu, w = windows, m = macintosh
taxon_tools_version = "Version 1.0.10 uw (r15.02.2020)"

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

##########################################################################################################################
# testing purpose

# TaXon_table_xlsx = "/home/till/Desktop/Projects_Development/MeTa_tools/TaXon_tools_data/TaXon_tables/JAMP_taxonomy_table_taxon_table_derep_no_NC.xlsx"

##########################################################################################################################
# define custom exceptions
class FilesDoNoMatch(Exception):
    pass

class TaxonomyTableError(Exception):
    pass

class ReadTableError(Exception):
    pass

class KronaInstallError(Exception):
    pass

class PCoAImportError(Exception):
    pass

##########################################################################################################################
# define functions for events

# Combine read table and taxonomyigger table
def TaXon_converter(Read_table_xlsx, taxonomy_results_xlsx, TaXon_table_name, path_to_outdirs):

    taxonomy_results_xlsx =  Path(taxonomy_results_xlsx)
    Read_table_xlsx = Path(Read_table_xlsx)
    # create filename and path for output file
    Output_name = TaXon_table_name + ".xlsx"
    Output_file = path_to_outdirs / "TaXon_tables" / Output_name

    # print information to console
    print("\n" + "Taxonomy table:", taxonomy_results_xlsx.name)
    print("Read table:", Read_table_xlsx.name)

    # store the file name for later use
    file_name = taxonomy_results_xlsx.name

    # load the taxonomy file as dataframe
    taxonomy_results_xlsx = pd.ExcelFile(taxonomy_results_xlsx)
    # load the Read table as dataframe
    Read_table_xlsx = pd.ExcelFile(Read_table_xlsx)
    # create datafrmes for both files
    taxonomy_df = pd.read_excel(taxonomy_results_xlsx, 'JAMP results', header=0).set_index('OTUs')
    Read_df = pd.read_excel(Read_table_xlsx, header=0).set_index('OTUs')
    # extract OTU names
    OTU_names_taxonomy_df = taxonomy_df.index.values
    OTU_names_Read_df = Read_df.index.values
    # create a sorting column -> important to match the input files
    taxonomy_df_sorting_col = [OTU.replace("OTU_", "") for OTU in OTU_names_taxonomy_df]
    Read_df_sorting_col = [OTU.replace("OTU_", "") for OTU in OTU_names_Read_df]
    taxonomy_df["sort"] = taxonomy_df_sorting_col
    Read_df["sort"] = Read_df_sorting_col
    # sort both dataframes by the sorting column
    taxonomy_df = taxonomy_df.sort_values(by=['sort'])
    Read_df = Read_df.sort_values(by=['sort'])

    # check if the sorting columns match each other -> if not break the and give an error message
    # raise a custom exception

    if taxonomy_df["sort"].values.tolist() != Read_df["sort"].values.tolist():
        raise FilesDoNoMatch("Files do not match")
        # print("Error!", "\n" + "    The OTU IDs do not match!", "\n" + "    >> Please check your files!")
        # exit()

    # create an empty list that will be used to create the output dataFrame
    dataframe_list = []
    # create a list of available OTUs (does not matter from which dataframe, because they have to be the same at this point)
    OTU_list = taxonomy_df.index.values

    # count the number of OTUs and start a counter for the progress bar
    number_of_OTUs = len(OTU_list) + 1
    progress_count = 0
    progress_counter = 1000 / number_of_OTUs

    # iterate through the OTUs (=row of dataframe) and increase the counter with each step for visual updates
    for i in range(len(OTU_list)):
        # get the name of the OTU
        OTU_name = OTU_list[i]
        # get the taxonomy from the taxonomy file
        taxonomy = taxonomy_df.iloc[i].values.tolist()
        # get all information from the Read table
        read_table = Read_df.iloc[i].values.tolist()
        # get the sequence information from the Read table
        sequence = read_table[-2]
        # get the read numbers from the Read table
        read_table = read_table[:-2]
        # get the sorting row
        sort = taxonomy_df["sort"][i]
        # append all collected information to a list from which the final dataframe will be created
        dataframe_list.append([OTU_name.replace(">", "")] + taxonomy[:-1] + sequence.split() + read_table + [int(sort)])
        # increase the progress bar
        progress_count += progress_counter
        window.Element('progbar').UpdateBar(progress_count)

    # create the column headers for the final file
    # get the names of the taxonomyigger file
    taxonomy_columns = list(taxonomy_df.columns)[:-1]
    # get the names of the Read table
    sample_names_columns = list(Read_df.columns)[:-2]
    # create the column headers
    header_row = ["OTUs"] + taxonomy_columns + ["seq"] + sample_names_columns + ["sort"]

    # create a new dataframe with the just created headers
    df = pd.DataFrame(dataframe_list)
    df.columns = header_row
    # set the index to be the OTUs
    df.set_index('OTUs')
    # sort the dataframe by the sorting column
    df = df.sort_values(by=["sort"])
    # remove the sorting column as it is no longer needed
    df = df.drop(["sort"], axis=1)
    # save the newly created Taxon table in TaXon format as excel file
    df.to_excel(Output_file, sheet_name='TaXon table', index=False)

# Create taxon tables per sample from a main taxon list
def create_TaXon_table_per_sample(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    TaXon_datasheet = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    TaXon_table = TaXon_datasheet.values.tolist()
    samples_to_process = TaXon_datasheet.columns[10:]
    first_ten_columns_header = TaXon_datasheet.columns[:10].values.tolist()
    first_ten_columns = TaXon_datasheet.iloc[:,[0,1,2,3,4,5,6,7,8,9]].values.tolist()
    OTU_list = TaXon_datasheet['OTUs'].values.tolist()

    sample_count = 0

    number_of_samples = len(samples_to_process) + 1
    progress_counter = 1000 / number_of_samples
    progress_count = progress_counter

    for sample in samples_to_process:

        print("Creating TaXon list for sample: ", sample)

        Output_name = Path(sample + ".xlsx")
        Output_file = path_to_outdirs / "TaXon_tables_per_sample" / Output_name

        read_numbers = TaXon_datasheet[sample].values.tolist()
        sample_rows_list = []

        for i, read_number in enumerate(read_numbers):
            if read_number > 0:
                sample_rows_list.append(first_ten_columns[i] + [read_number])

        headers_df = pd.DataFrame([first_ten_columns_header + [sample]])
        sample_df = pd.DataFrame(sample_rows_list)
        sample_df = headers_df.append(sample_df)

        sample_df.to_excel(Output_file, engine='xlsxwriter', sheet_name='TaXon table', index=False, header=False)

        sample_count += 1

        progress_count += progress_counter
        window.Element('progbar').UpdateBar(progress_count)

# Pie charts
def pie_chart(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    count = 0

    allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus", "F_Species"]

    for taxon_query in allowed_taxa:
        taxon_level = taxon_query
        taxon_query = taxon_query[2:]

        print("Calclulating pie chart for:", taxon_query, "level")

        count += 166
        window.Element('progbar').UpdateBar(count)

        dirName = Path(str(path_to_outdirs) + "/" + "Pie_charts" + "/" + str(TaXon_table_file.stem))

        if not os.path.exists(dirName):
            os.mkdir(dirName)

        Output_name = Path(str(dirName) + "/" + taxon_level + ".pdf")


        with PdfPages(Output_name) as pdf:

            taxonomy_results_xlsx = pd.ExcelFile(TaXon_table_file)
            data = pd.read_excel(taxonomy_results_xlsx, 'TaXon table', header=0)

            taxa = data[taxon_query].values.tolist()
            taxa_unique = list(dict.fromkeys(taxa))

            taxa_labels = []
            taxa_sizes = []

            for taxon_name in taxa_unique:
                taxa_labels.append(str(taxon_name))
                taxa_sizes.append(taxa.count(taxon_name))


            patches, texts, autotexts = plt.pie(taxa_sizes, labels=taxa_labels, autopct='%1.1f%%', startangle=140)

            for i in range(len(autotexts)):
                autotexts[i].set_fontsize(4)

            for i in range(len(texts)):
                texts[i].set_fontsize(4)

            plt.axis('equal')
            plt.title("OTU abundance on " + taxon_query + " level", fontsize=10, pad=20)
            fileName= TaXon_table_file.name
            plt.suptitle(fileName, x=0.2, y=0.05, fontsize=5)

            pdf.savefig()
            plt.close()

            #########################

            taxa_labels = []
            taxa_sizes = []

            for taxon_name in taxa_unique:
                if "nan" != str(taxon_name):
                    taxa_labels.append(str(taxon_name))
                    taxa_sizes.append(taxa.count(taxon_name))

            patches, texts, autotexts = plt.pie(taxa_sizes, labels=taxa_labels, autopct='%1.1f%%', startangle=140)

            for i in range(len(autotexts)):
                autotexts[i].set_fontsize(4)

            for i in range(len(texts)):
                texts[i].set_fontsize(4)

            plt.axis('equal')
            plt.title("OTU abundance (excluding nan) on " + taxon_query + " level", fontsize=10, pad=20)

            pdf.savefig()

            if taxon_query == "Species":
                plt.show(block=False)
                print("\n" + "Pie-charts are found in", path_to_outdirs, "Venn_diagrams/")
                sg.Popup("Pie-charts are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

            plt.close()

# Venn diagrams
def venn_diagram(file_a, file_b, path_to_outdirs):

    file_a = Path(file_a)
    file_b = Path(file_b)

    print("\n" + "Input file a:", file_a.stem)
    print("Input file b:", file_b.stem, "\n")
    count = 0

    allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species"]

    venn_dict = {}

    for taxon in allowed_taxa:

        output_name = taxon
        taxon = taxon[2:]

        data_file_a = pd.read_excel(file_a, 'TaXon table', header=0)
        data_file_b = pd.read_excel(file_b, 'TaXon table', header=0)

        file_name_a = file_a.stem
        file_name_b = file_b.stem

        taxa_file_a = data_file_a[taxon].values.tolist()
        taxa_file_b = data_file_b[taxon].values.tolist()

        taxa_unique_a = list(dict.fromkeys(taxa_file_a))
        taxa_unique_b = list(dict.fromkeys(taxa_file_b))

        taxa_labels_a = []
        taxa_labels_b = []
        taxa_sizes_a = []
        taxa_sizes_b = []

        for taxon_name in taxa_unique_a:
            if "nan" != str(taxon_name):
                taxa_labels_a.append(str(taxon_name))
                taxa_sizes_a.append(taxa_file_a.count(taxon_name))

        for taxon_name in taxa_unique_b:
            if "nan" != str(taxon_name):
                taxa_labels_b.append(str(taxon_name))
                taxa_sizes_b.append(taxa_file_b.count(taxon_name))

        taxa_labels_a = sorted(taxa_labels_a)
        taxa_labels_b = sorted(taxa_labels_b)

        a_only = set(taxa_labels_a) - set(taxa_labels_b)
        len_a_only = len(a_only)
        b_only = set(taxa_labels_b) - set(taxa_labels_a)
        len_b_only = len(b_only)
        shared = set(taxa_labels_a) & set(taxa_labels_b)
        len_shared = len(shared)

        venn_dict[taxon + "_a_only"] = a_only
        venn_dict[taxon + "_shared"] = shared
        venn_dict[taxon + "_b_only"] = b_only

        print("Comparing on:", taxon, "level")

        plt.figure(figsize=(20, 10))
        venn2(subsets = (len_a_only, len_b_only, len_shared), set_labels = (file_name_a, file_name_b))

        dirName = Path(str(path_to_outdirs) + "/" + "Venn_diagrams" + "/" + str(file_name_a) + "_VS_" + str(file_name_b) + "/")
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        output_pdf = Path(str(path_to_outdirs) + "/" + "Venn_diagrams" + "/" + file_name_a + "_VS_" + file_name_b + "/" + output_name + ".pdf")
        plt.title(taxon)
        plt.savefig(output_pdf, bbox_inches='tight')

        if taxon == "Species":
            plt.show(block=False)
            print("\n" + "Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/")
            sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

        plt.close()

        count += 166
        window.Element('progbar').UpdateBar(count)

    output_xlsx = Path(str(path_to_outdirs) + "/" + "Venn_diagrams" + "/" + file_name_a + "_VS_" + file_name_b + "/" + "Venn_comparison_results.xlsx")
    df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
    df.to_excel(output_xlsx, index=False)

## TaXon table filtering
#1.a create mask for user input
def create_taxa_mask(TaXon_table_xlsx, mask):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    data = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_taxa = set(data[mask].values.tolist())
    available_taxa = [x for x in available_taxa if str(x) != 'nan']
    available_taxa = sorted(list(available_taxa))

    return available_taxa

# 1.b filter taxa from taxon table
def taxon_table_filter(TaXon_table_xlsx, exclude_taxa, mask, appendix_name, threshold, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    # convert taxa to exclude to a list if only one taxon is given (which is then string)
    if type(exclude_taxa) == str:
        exclude_taxa = [exclude_taxa]

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    # check for taxa to filter
    mask_position = list(df.columns).index(mask)
    df_columns = df.columns

    rows_to_keep = []

    df_rows = df.values.tolist()
    for row in df_rows:
        taxon_to_evaluate = row[mask_position]
        if taxon_to_evaluate not in exclude_taxa:
            if str(taxon_to_evaluate) != 'nan':
                rows_to_keep.append(row)

    df_out = pd.DataFrame(rows_to_keep)

    similarity_position = list(df_columns).index("Similarity")
    threshold = int(threshold)

    filtered_rows = []

    for index, row in df_out.iterrows():
        similarity = list(row)[similarity_position]
        if similarity != 'No Match':
            if int(similarity) >= threshold:
                filtered_rows.append(list(row))

    df_out = pd.DataFrame(filtered_rows)

    if df_out.empty:
        sg.PopupError('Filter theshold were to harsh: Nothing to print', title="Error", keep_on_top=True)
        print('\n' + 'Filter theshold were to harsh: Nothing to print', '\n')

    else:
        df_out.columns = df_columns

        # write output file
        file_name = TaXon_table_file.stem
        output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")
        threshold_output = "Similarity threshold = " + str(threshold)
        exclude_taxa.append(threshold_output)
        df_exclude_taxa = pd.DataFrame(exclude_taxa)
        df_exclude_taxa.columns = ['Filter criteria']
        writer = pd.ExcelWriter(output_name, engine = 'xlsxwriter')
        df_out.to_excel(writer, sheet_name = 'TaXon table', index=False)
        df_exclude_taxa.to_excel(writer, sheet_name = 'Filter criteria', index=False)
        writer.save()
        writer.close()

        # increase the progress bar
        progress_count = 1000
        window.Element('progbar').UpdateBar(progress_count)

        print("\n" + "Filtered lists are found in \'TaXon tables\'/*_filtered")
        sg.Popup("Filtered lists are found in \'TaXon tables\'/*_filtered", title="Finished", keep_on_top=True)

# 2.a create mask for user input
def create_sample_mask(TaXon_table_xlsx):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_samples = df.columns.tolist()[10:]

    return available_samples

# 2.b filter samples from taxon list
def exclude_samples(TaXon_table_xlsx, exluded_samples, appendix_name, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    TaXon_table_xlsx_path = TaXon_table_xlsx
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    if type(exluded_samples) == str:
        exluded_samples = [exluded_samples]

    for sample in exluded_samples:
        df = df.drop(sample, axis=1)

    header = df.columns.values.tolist()

    row_filter_list = []

    for row in df.values.tolist():
        reads = set(row[10:])
        if reads != {0}:
            row_filter_list.append(row)

    df = pd.DataFrame(row_filter_list)
    df.columns = header

    file_name = TaXon_table_file.stem
    output_name = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + file_name + "_" + appendix_name + ".xlsx")
    df.to_excel(output_name, sheet_name = 'TaXon table', index=False)

    # increase the progress bar
    progress_count = 1000
    window.Element('progbar').UpdateBar(progress_count)

    print("\n" + "Filtered lists are found in \'TaXon tables\'/*_your_appendix")
    sg.Popup("Filtered lists are found in \'TaXon tables\'/*_your_appendix", title="Finished", keep_on_top=True)

# create taxon table
def create_taxon_list(TaXon_table_xlsx, taxon_list_output_file_name, language, values, taxon_tools_version, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    if taxon_list_output_file_name == '':
        taxon_list_output_file_name = "my_taxon_list"
    output_txt = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".txt")
    output_xlsx = Path(str(path_to_outdirs) + "/" + "Taxon_lists" + "/" + taxon_list_output_file_name + ".xlsx")

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    df_list = df.values.tolist()

    sample_names = df.columns[10:]
    project_list, waterID_list, site_list = [], [], []

    for sample in sample_names:
        print(sample)
        sample = sample.split("_")
        if len(sample) == 5:
            project_list.append(sample[0])
            waterID_list.append(sample[1])
            site_list.append(sample[2])

    # open document
    document = open(output_txt, "w")

    # add explanation to beginning of the file
    if language == "German":

        project = "Projekt(e): " + ', '.join(sorted(set(project_list)))
        waterID = "Fließgewässer ID(s): " + ', '.join(sorted(set(waterID_list)))
        sites = "Probestelle(n): " + ', '.join(sorted(set(site_list)))

        # if german: add explanation to taxon list
        document.write('Taxonliste ' + "(" + str(date.today()) + ")" + "\n")
        document.write('\n')
        document.write('Datenblatt:' + "\n")
        document.write(project + "\n")
        document.write(waterID + "\n")
        document.write(sites + "\n")
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
        document.write('Folgende Taxonliste enthält die durch DNA-Metabarcoding ermittelten Taxa. Die Anzahl der Operationellen Taxonomischen Einheiten (OTUs) für das jeweilige Taxon ist in Klammern hinter der Taxonomie gegeben.' + "\n")
        document.write('\n')
        document.write('Excel-Tabelle: ' + str(output_xlsx) + "\n")

    if language == "English":

        project = "Project(s): " + ', '.join(sorted(set(project_list)))
        waterID = "Water ID(s): " + ', '.join(sorted(set(waterID_list)))
        sites = "Sample site(s): " + ', '.join(sorted(set(site_list)))

        document.write('Taxon list ' + "(" + str(date.today()) + ")" + "\n")
        document.write('\n')
        document.write('Data sheet:' + "\n")
        document.write(project + "\n")
        document.write(waterID + "\n")
        document.write(sites + "\n")
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
        document.write('The following taxon list contains the via DNA metabarcoding identified taxa. The number of OTUs (operational taxonomic units) for each respective taxon are shown in parantheses behind the taxonomy.' + "\n")
        document.write('\n')
        document.write('Excel-Sheet: ' + str(output_xlsx))

    df_out = df[df.columns[1:7]]
    df_duplicates_removed = DataFrame.drop_duplicates(df_out)
    df_duplicates_removed.to_excel(output_xlsx, sheet_name = 'Taxa', index=False)

    document.close()

    # increase the progress bar
    progress_count = 1000
    window.Element('progbar').UpdateBar(progress_count)

# remove PCR or Extraction Replicates
def combine_replicates(TaXon_table_xlsx, replicate_version, path_to_outdirs):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    # create output file
    output_file = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + str(TaXon_table_file.stem) + "_derep.xlsx")
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    sample_names = df.columns[10:]
    unique_sample_names_list, samples_to_process_list = [], []

    for sample in sample_names:
        sample_name = sample.split("_")[0:4]
        unique_sample_names_list.append("_".join(sample_name))

    unique_sample_names_set = sorted(set(unique_sample_names_list))

    for sample in unique_sample_names_set:
        if len(sample.split("_")) == 4:

            if replicate_version == "extraction":
                rep_a = sample + "_a"
                rep_b = sample + "_b"
            if replicate_version == "PCR":
                rep_a = sample + "_A"
                rep_b = sample + "_B"

            combined = sample + "_comb"

            try:
                df[combined] = df[[rep_a, rep_b]].mean(axis=1)
                df = df.drop([rep_a, rep_b], axis=1)
                abort = False
            except KeyError:
                abort = True
                print("\n" + "Error: Did you choose the correct replicate version?")
                break

    if abort != True:
        df.to_excel(output_file, index=False, sheet_name = 'TaXon table')

        # increase the progress bar
        progress_count = 1000
        window.Element('progbar').UpdateBar(progress_count)

        text = "Combined " + replicate_version + " replicates."
        sg.Popup(text , title="Finished", keep_on_top=True)

    else:
        sg.PopupError("Error: Did you choose the correct replicate version?", keep_on_top=True)

# calculate rarefaction curve
def rarefaction(TaXon_table_xlsx, replicates, error_style, path_to_outdirs):

    TaXon_table_file = Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    output_file = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction.pdf")

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_samples = df.columns.tolist()[10:]
    sample_dict_clean = {}

    # iterate through all available samples
    for sample in available_samples:
        # create a dict for the read numbers of the respective sample
        sample_dict = df[sample].to_dict()
        # create an empty list for keys to store
        key_list = []
        # iterate through the read table of the sample
        for key, value in sample_dict.items():
            # if the read number of the OTU is NOT ZERO
            if value != 0:
                # append the OTU number to the key list
                key_list.append(key)
        # add the key list of each sample to a clean dictionary
        sample_dict_clean[sample] = key_list

    # draw once for each sample
    number_of_draws = len(sample_dict_clean.keys())

    # dictionary to store the drawing results
    draw_dictionary = {}

    # start the progress bar counter
    progress_count = 0
    # calculate the iterations
    progress_counter = 1000 / replicates

    for n_reps in range(0, replicates):
        # store the original dictionary to start over again
        # a copy of the original dictionary is required, because the samples will be removed with each draw
        # thus for each replicate a new dictionary to draw from has to be created
        sample_dict_to_draw = dict(sample_dict_clean)

        OTU_list = []
        OTU_set = []

        for i in range(0, number_of_draws):
            # choose a random sample from the dictionary
            random_choice = random.choice(list(sample_dict_to_draw.keys()))
            # extract the OTU IDs from the chosen sample and add them to the already existing OTU IDs
            OTU_list = OTU_list + sample_dict_clean[random_choice]
            # create a unique set
            OTU_set = set(OTU_list)
            # number of OTUs
            n_OTUs = len(OTU_set)
            # now add the unique OTU list to the output dictionary
            # if the key is not in the dict, create a new entry (= OTU ID plus number of OTUs)
            if i not in draw_dictionary.keys():
                draw_dictionary[i] = [n_OTUs]
            # if the key already exists, calculate the sum of the already existing number of OTUs and the new number of OTUs
            else:
                # create a new list to store the current number of OTUs
                add_OTUs_list = draw_dictionary[i]
                add_OTUs_list.append(n_OTUs)
                draw_dictionary[i] = add_OTUs_list

            # remove the sample to draw only once
            sample_dict_to_draw.pop(random_choice)

        # increase the progress bar
        progress_count += progress_counter
        window.Element('progbar').UpdateBar(progress_count)

    # create a dict to store the average number of OTUs per draw
    rarefaction_dict_average, rarefaction_dict_stdef = {}, {}

    def average(lst):
        return sum(lst) / len(lst)

    # iterate through the draw_dictionary and calculte the average number of OTUs
    for key, value in draw_dictionary.items():
        average_OTUs = average(draw_dictionary[key])
        stdef_OTUs = np.std(draw_dictionary[key], dtype=np.float64)
        rarefaction_dict_average[key] = average_OTUs
        rarefaction_dict_stdef[key] = stdef_OTUs

    if error_style == "a":
        draws = list(rarefaction_dict_average.keys())
        n_OTUs = list(rarefaction_dict_average.values())
        error_bar = list(rarefaction_dict_stdef.values())
        plt.figure(figsize=(20, 10))
        plt.errorbar(draws, n_OTUs, error_bar, linewidth=0.8, color='blue', capsize=1.3, capthick=1, ecolor='lightgrey')
        plt.xlabel('# samples')
        plt.ylabel('# OTUs')
        plt.title('repetitions = ' + str(replicates))

    elif error_style == "b":
        draws = list(rarefaction_dict_average.keys())
        n_OTUs = list(rarefaction_dict_average.values())
        y = np.asarray(n_OTUs)
        error_bar = np.asarray(list(rarefaction_dict_stdef.values()))
        plt.figure(figsize=(20, 10))
        plt.plot(draws, n_OTUs)
        plt.fill_between(draws, y-error_bar, y+error_bar, alpha=0.1)
        plt.xlabel('# samples')
        plt.ylabel('# OTUs')
        plt.title('repetitions = ' + str(replicates))


    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        print("\n" + "Rarefaction curves are found in", path_to_outdirs, "Rarefaction curves/")
        sg.Popup("Rarefaction curves are found in", path_to_outdirs, "Rarefaction curves/", title="Finished", keep_on_top=True)
    else:
        plt.close()

def calculate_taxonomic_resolution(TaXon_table_xlsx, path_to_outdirs, x_tax_res, y_tax_res, font_tax_res):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    output_file = Path(str(path_to_outdirs) + "/" + "Taxonomic_resolution_plots" + "/" + TaXon_table_file.stem + "_taxonomic_resolution.pdf")

    taxonomic_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    statistics_list, statistics_set, statistics_dict, highest_level_dict = [], [], {}, {}
    writer = pd.ExcelWriter(TaXon_table_xlsx, engine = 'xlsxwriter')

    for taxon_to_evaluate in taxonomic_levels:
        taxa_list = [x for x in df[taxon_to_evaluate].values.tolist() if str(x) != 'nan']
        statistics = taxon_to_evaluate, len(taxa_list)
        statistics_set.append(len(set(taxa_list)))
        statistics_list.append(list(statistics))
        statistics_dict[taxon_to_evaluate] = len(taxa_list)

    highest_level_dict["Phylum"] = statistics_dict["Phylum"] - statistics_dict["Class"]
    highest_level_dict["Class"] = statistics_dict["Class"] - statistics_dict["Order"]
    highest_level_dict["Order"] = statistics_dict["Order"] - statistics_dict["Family"]
    highest_level_dict["Family"] = statistics_dict["Family"] - statistics_dict["Genus"]
    highest_level_dict["Genus"] = statistics_dict["Genus"] - statistics_dict["Species"]
    highest_level_dict["Species"] = statistics_dict["Species"]

    taxon_levels = list(highest_level_dict.keys())
    highest_level_OTUs = list(highest_level_dict.values())
    total_OTUs = list(statistics_dict.values())

    plt.figure(figsize=(int(x_tax_res), int(y_tax_res)))

    plt.subplot(2, 1, 1)
    plt.plot(taxon_levels, highest_level_OTUs, '.-')
    #plt.bar(taxon_levels, highest_level_OTUs, align='center', alpha=0.5)
    plt.ylabel('# OTUs')
    plt.title('Taxonomic resolution (highest taxonomic level)')
    for i, v in enumerate(highest_level_OTUs):
        plt.text(i, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_res)

    plt.subplot(2, 1, 2)
    plt.plot(taxon_levels, total_OTUs, '.-')
    plt.ylabel('# OTUs')
    plt.title('Number of OTUs per taxonomic level')
    for i, v in enumerate(total_OTUs):
        plt.text(i, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_res)

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, dpi=500)
        plt.close()
        print("\n" + "Taxonomic resolution plots are found in", path_to_outdirs, "Taxonomic_resolution_plots/")
        sg.Popup("Taxonomic resolution plots are found in", path_to_outdirs, "Taxonomic_resolution_plots/", title="Finished", keep_on_top=True)
    else:
        plt.close()

    # increase the progress bar
    progress_count = 1000
    window.Element('progbar').UpdateBar(progress_count)

def calculate_taxonomic_richness(TaXon_table_xlsx, path_to_outdirs, x_tax_rich, y_tax_rich, font_tax_rich):

    TaXon_table_file =  Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    output_file = Path(str(path_to_outdirs) + "/" + "Taxonomic_richness_plots" + "/" + TaXon_table_file.stem + "_taxonomic_richness.pdf")

    taxonomic_levels = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    statistics_list, statistics_set, statistics_dict, highest_level_dict = [], [], {}, {}
    writer = pd.ExcelWriter(TaXon_table_xlsx, engine = 'xlsxwriter')

    for taxon_to_evaluate in taxonomic_levels:
        taxa_list = [x for x in df[taxon_to_evaluate].values.tolist() if str(x) != 'nan']
        statistics = taxon_to_evaluate, len(taxa_list)
        statistics_set.append(len(set(taxa_list)))
        statistics_list.append(list(statistics))
        statistics_dict[taxon_to_evaluate] = len(taxa_list)

    highest_level_dict["Phylum"] = statistics_dict["Phylum"] - statistics_dict["Class"]
    highest_level_dict["Class"] = statistics_dict["Class"] - statistics_dict["Order"]
    highest_level_dict["Order"] = statistics_dict["Order"] - statistics_dict["Family"]
    highest_level_dict["Family"] = statistics_dict["Family"] - statistics_dict["Genus"]
    highest_level_dict["Genus"] = statistics_dict["Genus"] - statistics_dict["Species"]
    highest_level_dict["Species"] = statistics_dict["Species"]

    taxon_levels = list(highest_level_dict.keys())
    number_of_taxa_per_level = statistics_set

    plt.figure(figsize=(int(x_tax_rich), int(y_tax_rich)))

    plt.plot(taxon_levels, number_of_taxa_per_level, '.-')
    plt.ylabel('# taxa')
    plt.title('Number of taxa per taxonomic level')
    for i, v in enumerate(number_of_taxa_per_level):
        plt.text(i - 0.01, v, str(v), horizontalalignment='center', verticalalignment='center', fontsize=font_tax_rich)

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()

        print("\n" + "Taxonomic richness plots are found in", path_to_outdirs, "taxonomic_richness_plots/")
        sg.Popup("Taxonomic richness plots are found in", path_to_outdirs, "taxonomic_richness_plots/", title="Finished", keep_on_top=True)
    else:
        plt.close()


    # increase the progress bar
    progress_count = 1000
    window.Element('progbar').UpdateBar(progress_count)

def create_metadata_table(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    TaXon_table_xslx_df = pd.read_excel(TaXon_table_xlsx)
    samples_list = TaXon_table_xslx_df.columns.tolist()[10:]
    samples_metadata_list = []

    for sample in samples_list:
        if len(sample.split("_")) == 5:
            project = sample.split("_")[0]
            site = sample.split("_")[1]
            sample_id = sample.split("_")[2]
            date = str(sample.split("_")[3])
            sample_metadata = [sample, project, site, sample_id, date]
            samples_metadata_list.append(sample_metadata)

        else:
            print("Skipping", sample, "- false format")
            sample_metadata = [sample, "nan", "nan", "nan", "nan"]
            samples_metadata_list.append(sample_metadata)

    metadata_df = pd.DataFrame(samples_metadata_list)
    metadata_df = metadata_df.rename({0: 'Samples', 1: 'Project', 2: 'Site', 3: 'Sample', 4: 'Date'}, axis='columns')
    metadata_df.to_excel(Meta_data_table_xlsx, index=False)

    # increase the progress bar
    progress_count = 1000
    window.Element('progbar').UpdateBar(progress_count)

# site occupancy
# 1
def get_available_meta_data(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    meta_data_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    if os.path.exists(meta_data_xlsx):

        meta_data_xlsx = pd.ExcelFile(meta_data_xlsx)
        data = pd.read_excel(meta_data_xlsx, header=0)

        available_meta_data = data.columns.tolist()[1:]

        return available_meta_data

    else:
        return False
# 2
def site_occupancy(TaXon_table_xlsx, meta_data_to_test, taxonomic_level, path_to_outdirs, x_site_occ, y_site_occ, font_site_occ):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header = 0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]

    TaXon_table_n_samples = len(TaXon_table_samples)

    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    answer = "Ask"
    output_message = "No"

    try:
        if os.path.exists(Meta_data_table_xlsx):

            Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header = 0)
            Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

            n_sites = len(set(Meta_data_table_df[meta_data_to_test].tolist()))

            if (TaXon_table_samples == Meta_data_table_samples and TaXon_table_n_samples != n_sites):

                site_occupancy_dict = {}

                sites = set(Meta_data_table_df[meta_data_to_test].tolist())

                for site in sites:
                    # this can either be a species name or the above specified taxonomic level
                    present_OTU_list = []

                    # extract samples that belong to the site from the Meta data file
                    included_samples_list = Meta_data_table_df[Meta_data_table_df.values  == site]['Samples'].values.tolist()

                    # count the number of samples per site to calculate the site occupancy
                    n_samples = len(included_samples_list)

                    # create a list of all species (or the specified taxonomic level)
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

                    # store occupancy of each species in a list, will be accessed by position in list
                    species_list = []
                    occupancy_list = []

                    # count the number of occurences for each species and calculate the occpancy based on the number of samples
                    for species in overall_included_species_set:
                        count = present_OTU_list_flattened.count(species)
                        occupancy = count / n_samples
                        species_list.append(species)
                        occupancy_list.append(occupancy)

                    # create a simple bar plot for each site
                    plt.figure(figsize=(int(x_site_occ), int(y_site_occ)))
                    plt.bar(species_list, occupancy_list)
                    plt.title(site + " (" + taxonomic_level + ")")
                    plt.xticks(fontsize=font_site_occ, rotation=90)
                    plt.yticks(fontsize=font_site_occ)

                    occupancy_plot_directory = Path(str(path_to_outdirs) + "/" + "Site_occupancy_plots" + "/" + TaXon_table_xlsx.stem)
                    if not os.path.exists(occupancy_plot_directory):
                        os.mkdir(occupancy_plot_directory)

                    occupancy_plot_pdf = Path(str(occupancy_plot_directory) + "/" + site + "_" + taxonomic_level + ".pdf")

                    if answer == "Ask":
                        # show the first plot
                        plt.show(block=False)
                        # ask for answer
                        answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
                        # if answer is yes, save all upcoming plots
                        if answer == "Yes":
                            output_message = "Yes"
                            plt.savefig(occupancy_plot_pdf)
                            plt.close()
                        else:
                            # answer is "No"
                            plt.close()
                    elif answer == "Yes":
                            plt.savefig(occupancy_plot_pdf)
                            plt.close()
                    else:
                            plt.close()

                    # increase the progress bar
                    progress_count = 1000
                    window.Element('progbar').UpdateBar(progress_count)

                if output_message == "Yes":
                    print("\n" + "Site occupancy plots are found in", path_to_outdirs, "Site_occupancy_plots/")
                    sg.Popup("Site occupancy plots are found in", path_to_outdirs, "Site_occupancy_plots/", title="Finished", keep_on_top=True)

            else:
                sg.PopupError("Please check your Meta data file and Taxon table file: The samples do not match or the meta data is unique for all samples!", keep_on_top=True)
                print("Please check your Meta data file and Taxon table file: The samples do not match or the meta data is unique for all samples!")

        else:
            sg.PopupError("Missing meta data file!", keep_on_top=True)
            print("Missing meta data file!")

    except:
        sg.PopupError("Something went wrong! Check your meta data file and read the manual", title="Error", keep_on_top=True)
        print("Something went wrong! Check your meta data file and read the manual")

# read read_proportions
def read_proportions(TaXon_table_xlsx, taxonomic_level, path_to_outdirs, x_read_props_scatter, y_read_props_scatter, font_read_props_scatter):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)
    samples_list = TaXon_table_df.columns.tolist()[10:]
    Species_read_proportion_dict = {}

    ##############################################################################
    # scatter plot

    read_proportion_plot_directory = Path(str(path_to_outdirs) + "/" + "Read_proportions_plots" + "/" + TaXon_table_xlsx.stem)
    scatter_pdf = Path(str(read_proportion_plot_directory) + "_" + taxonomic_level + "_scatter.pdf")

    # start the progress bar counter
    progress_count = 0
    # calculate the iterations
    progress_counter = 1000 / len(samples_list)

    for sample in samples_list:

        df = TaXon_table_df[['OTUs', taxonomic_level, sample]]
        df_2 = df[[sample]]/df[[sample]].sum()
        df = df.assign(perc=df_2.values)

        TaXon_table_list = df.values.tolist()

        for line in TaXon_table_list:

            species = line[1]

            reads = line[3] * 100
            dict_key = sample, species

            if dict_key not in Species_read_proportion_dict.keys():
                Species_read_proportion_dict[dict_key] = reads
            else:
                Species_read_proportion_dict[dict_key] = Species_read_proportion_dict[dict_key] + reads

        # increase the progress bar
        progress_count += progress_counter
        window.Element('progbar').UpdateBar(progress_count)

    sample_list = []
    species_list = []
    read_proportion_list = []

    for dict_key, value in list(Species_read_proportion_dict.items()):
        sample_list.append(dict_key[0])
        species_list.append(dict_key[1])
        read_proportion_list.append(value)

    # create data
    x = sample_list #['a', 'b', 'c']
    y = species_list #['d', 'f', 'g']-
    z = read_proportion_list #[199, 290, 290]

    # calculate the number of taxa
    n_taxa = len(set(species_list))
    # calculate the number of scatter entries for the plot
    n_species_list = len(species_list)

    # create an even color list if the number of taxa is even
    if (n_taxa % 2) == 0:
       print("even")
       clist = (["#52af2fff", "#4472c4ff"] * int(n_species_list / 2))
    # create an odd color list if the number of taxa is odd
    # this every sample will start with the same color
    else:
       print("odd")
       clist = (["#52af2fff", "#4472c4ff"] * (int((n_taxa - 1) / 2)) + ["#52af2fff"]) * n_taxa
       clist = clist[:n_species_list]

    title = "Read proportions (" + taxonomic_level + ")"
    fig, ax = plt.subplots(figsize=(int(x_read_props_scatter), int(y_read_props_scatter)))
    ax.scatter(x, y, s=z, color=clist, alpha=0.9)
    plt.xticks(fontsize=font_read_props_scatter, rotation=90)
    plt.yticks(fontsize=font_read_props_scatter)
    plt.title(title, fontsize=8)
    ax.set_axisbelow(True)
    ax.grid(color='gray', alpha=0.1)

    ########################

    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        print("Saving...")
        plt.savefig(scatter_pdf)
        plt.close()
        print("\n" + "Read proportion plots are found in", path_to_outdirs, "Read_proportion_plots/")
        sg.Popup("Read proportion plots are found in", path_to_outdirs, "Read_proportion_plots/", title="Finished", keep_on_top=True)
    else:
        plt.close()

    ########################

    # ##############################################################################
    # # stacked barchart
    #
    # barchart_pdf = Path(str(read_proportion_plot_directory) + "_" + taxonomic_level + "_barchart.pdf")
    #
    # relative_read_abundances_dict = {}
    # samples_set = set(sample_list)
    # included_samples_list = []
    #
    # for sample in samples_set:
    #         entry_list = TaXon_table_df[[sample, taxonomic_level]].values.tolist()
    #         entry_dict = {}
    #         for entry in entry_list:
    #             try:
    #                 entry_dict[entry[1]] = entry_dict[entry[1]] + entry[0]
    #             except KeyError:
    #                 entry_dict[entry[1]] = entry[0]
    #
    #         absolute_read_abundances_list = list(entry_dict.values())
    #         relative_read_abundances_list = []
    #         absolute_read_abundances_sum = sum(absolute_read_abundances_list)
    #
    #         if absolute_read_abundances_sum != 0:
    #             for entry in absolute_read_abundances_list:
    #                 relative_read_abundance = entry / absolute_read_abundances_sum * 100
    #                 relative_read_abundances_list.append(relative_read_abundance)
    #
    #             relative_read_abundances_dict[sample] = relative_read_abundances_list
    #
    #             included_samples_list.append(sample)
    #
    #         else:
    #             print("Skipping: ", sample)
    #
    #         # increase the progress bar
    #         progress_count += progress_counter
    #         window.Element('progbar').UpdateBar(progress_count)
    #
    #
    # taxa_list = list(entry_dict.keys())
    #
    # relative_read_abundances_df = pd.DataFrame.from_dict(relative_read_abundances_dict)
    # relative_read_abundances_df.index = taxa_list
    #
    # data = np.array(relative_read_abundances_df.values.tolist())
    # n = len(relative_read_abundances_df.values.tolist())
    # color=iter(cm.tab20_r(np.linspace(0,1,n)))
    # X = np.arange(data.shape[1])
    #
    # plt.figure(figsize=(int(x_read_props_scatter)*2, int(y_read_props_scatter)*2))
    #
    # for i in range(data.shape[0]):
    #   c=next(color)
    #   plt.bar(X, data[i], bottom = np.sum(data[:i], axis = 0), color = c, edgecolor='white')
    #
    # plt.legend(taxa_list, loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.xticks(X, included_samples_list, rotation=90)
    # plt.savefig(barchart_pdf)
    # plt.close()

# basic statistics
def basic_stats(TaXon_table_xlsx, path_to_outdirs):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, 'nan', regex=True)

    # number of samples
    n_samples = len(TaXon_table_df.columns[10:].tolist())

    # number of OTUs
    n_OTUs = len(TaXon_table_df['OTUs'].tolist())

    # number of taxa per taxon level
    n_Phyla = len(set(TaXon_table_df['Phylum'].tolist()))
    n_Classes = len(set(TaXon_table_df['Class'].tolist()))
    n_Orders = len(set(TaXon_table_df['Order'].tolist()))
    n_Families = len(set(TaXon_table_df['Family'].tolist()))
    n_Genera = len(set(TaXon_table_df['Genus'].tolist()))
    n_Species = len(set(TaXon_table_df['Species'].tolist()))

    # number of respective status
    status_dict = {}
    status_entries_set = set(TaXon_table_df['Status'].tolist())
    for status in status_entries_set:
        count = TaXon_table_df['Status'].tolist().count(status)
        status_dict[status] = count

    # sequence lengths
    sequence_list = TaXon_table_df['seq'].tolist()
    sequence_len_list = []
    for sequence in sequence_list:
        sequence_len_list.append(len(sequence))
    sequence_len_set = set(sequence_len_list)
    min_len_seq = min(sequence_len_set)
    max_len_seq = max(sequence_len_set)
    avg_len_seq = round(sum(sequence_len_set) / len(sequence_len_set))

    # read stats per sample
    samples = TaXon_table_df.columns[10:].tolist()
    reads_dict = {}
    for sample in samples:
        reads_list = TaXon_table_df[sample].tolist()
        reads_sum = sum(reads_list)
        reads_avg = round(sum(reads_list) / len(reads_list))
        reads_dict[sample] = [reads_sum, reads_avg]

    # read stats total
    read_sum_total = 0
    for read_sum in reads_dict.values():
        read_sum_total += read_sum[0]

    #####################################################################################
    output_list = []
    #df = pd.DataFrame(simple_list,columns=['col1','col2'])

    output_list.append(['Samples',str(n_samples),''])
    output_list.append(['OTUs',str(n_OTUs),''])
    output_list.append(['Number of taxa per taxon level','',''])
    output_list.append(['Phyla',str(n_Phyla),''])
    output_list.append(['Classes',str(n_Classes),''])
    output_list.append(['Orders',str(n_Orders),''])
    output_list.append(['Families',str(n_Families),''])
    output_list.append(['Genera',str(n_Genera),''])
    output_list.append(['Species',str(n_Species),''])
    output_list.append(['Database status','',''])
    for status, count in status_dict.items():
        output_list.append([status,count,''])
    output_list.append(['Sequence length (bp)','',''])
    output_list.append(['Min', min_len_seq,''])
    output_list.append(['Avg', avg_len_seq,''])
    output_list.append(['Max', max_len_seq,''])
    output_list.append(['Read statistics (sample, avg, sum)','',''])
    for sample, reads_stats in reads_dict.items():
        output_list.append([sample, str(reads_stats[1]), str(reads_stats[0])])
    output_list.append(['Total reads', '', str(read_sum_total)])

    df = pd.DataFrame(output_list,columns=['A','B','C'])

    basic_stats_directory = Path(str(path_to_outdirs) + "/" + "Basic_stats" + "/" + TaXon_table_xlsx.stem)
    basic_stats_xlsx = Path(str(basic_stats_directory) + "_basic_stats.xlsx")
    df.to_excel(basic_stats_xlsx, index=False)

    sg.Popup(df.to_string(index=False, header=False), title="Basic stats", keep_on_top=True)

# check the input format
def check_taxon_table_format(taxonomy_results_xlsx):

    try:
        taxonomy_table_df = pd.read_excel(Path(taxonomy_results_xlsx), 'JAMP results')
        taxonomy_table_df = taxonomy_table_df.replace(np.nan, 'nan', regex=True)
    except:
        sg.PopupError("Taxonomy table must be in excel sheet format (.xlsx)", keep_on_top=True)
        raise TaxonomyTableError("Taxonomy table must be in excel sheet format (.xlsx)")

    ###################################
    # A) header prompt
    taxonomy_table_df_header = taxonomy_table_df.columns.tolist()
    header_prompt = ["OTUs", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Similarity", "Status"]

    if taxonomy_table_df_header != header_prompt:
        sg.PopupError("Oops! Something is wrong with the header!", title="Error", keep_on_top=True)
        raise TaxonomyTableError("Oops! Something is wrong with the header!")

    ###################################
    # B) OTUs prompt

    OTU_list = taxonomy_table_df['OTUs'].values.tolist()

    # loop through all available OTUs
    for OTU in OTU_list:

        # define a variable Error message
        ErrorMessage = "Oops! Something is wrong with the OTU ID: " + OTU

        try:
            # split the OTU, which should be containing of two elements
            p1 = OTU.split("_")[0]
            p2 = OTU.split("_")[1]

        except:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise TaxonomyTableError(ErrorMessage)

        # check if the first part in named "OTU"
        if p1 != "OTU":
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise TaxonomyTableError(ErrorMessage)

        # check if the second part is an integer
        # convert p2 to an integer
        p2 = int(p2)
        # now check if it is an integer
        if type(int(p2)) is not int:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise TaxonomyTableError(ErrorMessage)

    ###################################
    # C) Species prompt
    # does not necessarily need to include both the Genus and epithet, but it is recommended

    species_list = taxonomy_table_df['Species'].values.tolist()
    epithet_only_list = []

    for species in species_list:
        if (species != 'nan' and len(species.split()) < 2):
            epithet_only_list.append(species)

    if epithet_only_list != []:
        WarningMessage = "Warning: There are species that do not fit the binomial nomenclature:" + "\n" + "\n" + ', '.join(epithet_only_list) + "\n" + "\n" + "This error message will be ignored, but it's recommended to use binomial nomenclature!"
        sg.Popup(WarningMessage, title="Warning", keep_on_top=True)

    ###################################
    # Wrap up

    sg.Popup("Your file looks great and is ready to use!", title="Taxonomy table check", keep_on_top=True)

def check_read_table_format(Read_table_xlsx):

    try:
        read_table_df = pd.read_excel(Path(Read_table_xlsx))
        read_table_df = read_table_df.replace(np.nan, 'nan', regex=True)
    except:
        sg.PopupError("Taxonomy table must be in excel sheet format (.xlsx)", keep_on_top=True)
        raise ReadTableError("Taxonomy table must be in excel sheet format (.xlsx)")

    ###################################
    # A) header prompt
    read_table_df_header = read_table_df.columns.tolist()
    header_prompt_OTUs = read_table_df_header[0]
    header_prompt_sequences = read_table_df_header[-1]

    if header_prompt_OTUs != "OTUs":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_OTUs + "\n" + "\n" + "Prompt: OTUs"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise ReadTableError(ErrorMessage)

    if header_prompt_sequences != "Sequences":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_sequences + "\n" + "\n" + "Prompt: Sequences"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise ReadTableError(ErrorMessage)

    ###################################
    # B) OTUs prompt
    OTU_list = read_table_df['OTUs'].values.tolist()

    # loop through all available OTUs
    for OTU in OTU_list:

        # define a variable Error message
        ErrorMessage = "Oops! Something is wrong with the OTU ID: " + OTU

        try:
            # split the OTU, which should be containing of two elements
            p1 = OTU.split("_")[0]
            p2 = OTU.split("_")[1]

        except:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise ReadTableError(ErrorMessage)

        # check if the first part in named "OTU"
        if p1 != "OTU":
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise ReadTableError(ErrorMessage)

        # check if the second part is an integer
        # convert p2 to an integer
        p2 = int(p2)
        # now check if it is an integer
        if type(int(p2)) is not int:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise ReadTableError(ErrorMessage)

    ###################################
    # C) Samples prompt

    samples_list = read_table_df.columns[1:-1].tolist()
    too_short_list = []
    wrong_appendix_list = []

    for sample in samples_list:
        try:
            sample_split = sample.split("_")
        except:
            print(sample)

        if len(sample_split) != 5:
            too_short_list.append(sample)

        elif sample_split[-1] not in ["a", "b", "A", "B"]:
            wrong_appendix_list.append(sample)

    if too_short_list != [] or wrong_appendix_list != []:
        WarningMessage = "Samples that do not fit the prompt: " + "\n" + ', '.join(too_short_list) + "\n" + "\n" + "Samples that do not fit replicates prompt: " + "\n" + ', '.join(wrong_appendix_list) + "\n" + "\n" + "This error message will be ignored, but it's recommended to use the samples prompt!"
        sg.Popup(WarningMessage, title="Warning", keep_on_top=True)

    ###################################
    # Wrap up
    sg.Popup("Your file looks great and is ready to use!", title="Read table check", keep_on_top=True)

# krona chart
def create_krona_chart(TaXon_table_xlsx, path_to_outdirs):

    try:
        subprocess.call(["ktImportText"], stdout=open(os.devnull, 'wb'))
    except:
        sg.PopupError("You need to install Krona tools first!" + "\n" * 2 + "Krona tools is currently not supported on Windows!" + "\n", title="Error")
        raise KronaInstallError("Krona tools needs to be installed")

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.replace(np.nan, '__', regex=True)

    row1 = ["sample-ID", "", "", "", "", "", ""]
    row2 = ["count", "phylum", "class", "order", "family", "genus", "species"]
    krona_taxonomy_list = []
    krona_taxonomy_list.append(row1)
    krona_taxonomy_list.append(row2)

    for OTU in TaXon_table_df.values.tolist():
        taxonomy = OTU[1:7]
        reads = sum(OTU[10:])
        krona_taxonomy_list.append([reads] + taxonomy)

    krona_taxonomy_df = pd.DataFrame(krona_taxonomy_list)

    krona_chart_directory = Path(str(path_to_outdirs) + "/" + "Krona_charts" + "/" + TaXon_table_xlsx.stem)
    krona_table_tsv = Path(str(krona_chart_directory) + "_krona_table.tsv")
    krona_chart_html = Path(str(krona_chart_directory) + "_krona_table.html")

    # write krona table to tsv
    krona_taxonomy_df.to_csv(krona_table_tsv, sep="\t", header=False, index=False)

    os.system("ktImportText " + str(krona_table_tsv) + " -o " + str(krona_chart_html))

    print("\n" + "Krona charts are found in", path_to_outdirs, "Krona_charts/")
    sg.Popup("Krona charts are found in", path_to_outdirs, "Krona_charts/", title="Finished", keep_on_top=True)

# PCoA analysis
def PCoA_analysis(TaXon_table_xlsx, meta_data_to_test, path_to_outdirs):

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    Output_PCoA_plot_pdf = Path(str(path_to_outdirs) + "/" + "PCoA_plots" + "/" + TaXon_table_xlsx.stem + ".pdf")

    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header = 0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]


    if TaXon_table_samples == Meta_data_table_samples:

        samples = TaXon_table_samples
        OTU_abundances = []
        samples_metadata_list = []
        remove_list= []

        # remove samples that do not fit the format
        for i, sample in enumerate(samples):
            meta_data = str(Meta_data_table_df.loc[i][meta_data_to_test])
            samples_metadata_list.append([meta_data])

        if len(set(tuple(i) for i in samples_metadata_list)) > 1:

            try:

                # calculate the OTU abundances
                for sample in samples:
                    OTU_abundances.append(TaXon_table_df[sample].values.tolist())

                #################################
                # First, we’ll compute observed OTUs, an alpha diversity metric, for each sample using the alpha_diversity driver function:

                adiv_obs_otus = alpha_diversity('observed_otus', OTU_abundances, samples)

                #################################
                # calculate jaccard distances

                j_dm = beta_diversity("jaccard", OTU_abundances, samples)

                #################################
                # calculate braycurtis distances

                bc_dm = beta_diversity("braycurtis", OTU_abundances, samples)

                #################################
                # Next we’ll do some work with these beta diversity distance matrices. First, we’ll determine if the Jaccard and
                # Bray-Curtis distance matrices are significantly correlated by computing the Mantel correlation between them.
                # Then we’ll determine if the p-value is significant based on an alpha of 0.05.

                r, p_value, n = mantel(j_dm, bc_dm)
                print("Mantel test: Correlation between Jaccard and Bray-Curtis distance matrices")
                print(" r =", round(r, 4))
                print(" p =", p_value)

                print()

                #################################
                # Next, we’ll perform principal coordinates analysis (PCoA) on our Jaccard distance matrix.
                j_pc = pcoa(j_dm)

                # create a new dataframe that includes the sample name (as index), the site and the sample number samples = ["A", "B", "C", "D", "E"]
                sample_md = pd.DataFrame(samples_metadata_list, index=samples, columns=[meta_data_to_test])

                figure = j_pc.plot(sample_md, meta_data_to_test, axis_labels=('PC 1', 'PC 2', 'PC 3'), title=meta_data_to_test, cmap='rainbow', s=50)

                print(samples_metadata_list)

                #################################
                # Let’s use ANOSIM to quantify the strength of the clustering we see in the ordination plots above, using our weighted Jaccard distance matrix and sample metadata.
                results = anosim(j_dm, sample_md, column=meta_data_to_test, permutations=999)
                print("PCoA clustering strength test")
                ANOSIM = round(results['test statistic'], 4)
                p_value = results['p-value']
                print(" ANOSIM =", ANOSIM)
                print(" p =", p_value)

                # add the test to the plot
                text = "ANOSIM = " + str(ANOSIM) + "\n" + "p = " + str(p_value)
                figure.text(0.8, 0, text)
                # figure.savefig(Output_PCoA_plot_pdf, bbox_inches='tight')
                # plt.close(figure)
                plt.draw()
                plt.pause(0.001)
                answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
                if answer == "Yes":
                    figure.savefig(Output_PCoA_plot_pdf, bbox_inches='tight')
                    plt.close(figure)
                    print("\n" + "Taxonomic richness plots are found in", path_to_outdirs, "taxonomic_richness_plots/")
                    sg.Popup("Taxonomic richness plots are found in", path_to_outdirs, "taxonomic_richness_plots/", title="Finished", keep_on_top=True)
                else:
                    plt.close(figure)

            except ValueError:
                raise
                sg.PopupError("Error: All values in the grouping vector are the same or there are floats!")
                print("Error: All values in the grouping vector are the same or there are floats")

        else:
            sg.PopupError("Error: All values in the grouping vector are the same or there are floats!")
            print("Error: All values in the grouping vector are the same or there are floats")
    else:
        sg.PopupError("Error: Please check your Meta data file and Taxon table file: The samples do not match!")
        print("Error: Please check your Meta data file and Taxon table file: The samples do not match!")

# taxon table consistency filtered
def replicate_consistency_filter(TaXon_table_xlsx, path_to_outdirs):

    try:

        TaXon_table_xlsx = Path(TaXon_table_xlsx)
        TaXon_table_df = pd.read_excel(TaXon_table_xlsx)

        samples = TaXon_table_df.columns[10:].tolist()
        OTUs = TaXon_table_df["OTUs"].values.tolist()

        derep_sample_names_dict =  {}

        # get all replicates from the file and write them to a dict
        # skip all samples that don't fit the format

        for sample in samples:
            sample = sample.split("_")
            sample_stem = "_".join(sample[:-1])
            sample_name = "_".join(sample)

            if len(sample) == 5:
                if "_".join(sample[:-1]) not in derep_sample_names_dict:
                    derep_sample_names_dict[sample_stem] = sample_name
                else:
                    values = [derep_sample_names_dict[sample_stem]] + [sample_name]
                    derep_sample_names_dict[sample_stem] = values
            else:
                print("Skipping:", "_".join(sample))

        for key, values in derep_sample_names_dict.items():
            new_df = TaXon_table_df[values]
            header = new_df.columns.tolist()
            processed_reads = []

            for n_reads in new_df.values.tolist():
                if 0 in n_reads:
                    if len(set(n_reads)) > 1:
                        n_reads = len(n_reads) * [0]
                processed_reads.append(n_reads)

            df_out = pd.DataFrame(processed_reads)
            df_out.columns = header
            TaXon_table_df = TaXon_table_df.drop(values, axis=1)
            TaXon_table_df[values] = df_out

        taxon_tables_directory = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem)
        taxon_tables_cons_xlsx = Path(str(taxon_tables_directory) + "_cons.xlsx")

        TaXon_table_df.to_excel(taxon_tables_cons_xlsx, sheet_name='TaXon table', index=False)

        print("\n" + "Taxon tables are found in", path_to_outdirs, "TaXon_tables/")
        sg.Popup("Taxon tables are found in", path_to_outdirs, "TaXon_tables/", title="Finished", keep_on_top=True)

    except:
        sg.PopupError("Something went wrong! Did you use the correct sample promt?", keep_on_top=True)
        print("Something went wrong! Did you use the correct sample promt?")

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

sg.ChangeLookAndFeel('Reddit')

##########################################################################################################################
# start Popup window

# check for already existing projects
projects = glob.glob('Projects/*')
projects_list = []

for project in projects:
    projects_list.append(Path(project).stem)

projects_radio = list(slices([sg.Radio(name, "projects", default=True) for name in sorted(projects_list)], 2))

start_window_layout = [
            [sg.Text('',size=(1,1))],
            [sg.Text('Create new project folder:')],
            [sg.Input('', key='new_project_folder', size=(40,1))],
            [sg.Text('',size=(1,1))],
            [sg.Button('Create new')],
            [sg.Text('',size=(1,4))],
            [sg.Text('Load existing project folder')],
            [sg.Frame(layout = projects_radio, title = '')],
            [sg.Text('',size=(1,1))],
            [sg.Button('Load')],
            [sg.Text('',size=(1,4))],
            [sg.Button('Exit', button_color=('black', 'red'))],
            ]

start_window = sg.Window('Projects', start_window_layout, keep_on_top=True)
event, values = start_window.read()

while True:

    if event == 'Create new':
        if values["new_project_folder"] != '':
            project_folder = values["new_project_folder"].replace(" ", "_")
            break
        else:
            project_folder = "Default_project"
            break

    if event == 'Load':
        project_folder = ''
        for key, value in values.items():
            if value == True:
                project_folder = sorted(projects_list)[key]
        if project_folder == '':
            project_folder = "Default_project"
        break

    if event == 'Exit':
        sys.exit()

start_window.close()

##########################################################################################################################
# check folders and create folders if neccessary

if not os.path.exists("Projects"):
    os.mkdir(Path("Projects"))

path_to_outdirs = Path("Projects" + "/" + project_folder)
if not os.path.exists(path_to_outdirs):
    os.mkdir(path_to_outdirs)

directories_to_create = ["Pie_charts", "Pie_charts", "Venn_diagrams","TaXon_tables", "TaXon_tables_per_sample", "Taxon_lists", "Rarefaction_curves",
"Taxonomic_richness_plots", "PCoA_plots", "Taxonomic_resolution_plots", "Meta_data_table", "Site_occupancy_plots", "Read_proportions_plots", "Basic_stats", "Krona_charts"]

for directory in directories_to_create:
    dirName = Path("Projects" + "/" + project_folder + "/" + directory + "/")
    if not os.path.exists(dirName):
        os.mkdir(Path(dirName))

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

data_conversion = [
					[sg.Text('',size=(1,1))],
					[sg.Text('TaXon converter', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Enter path to taxonomy table and read table to create a taxon table in TaXon format', size=(85,2))],
					[sg.Text('Taxonomy table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxonomy_table_path'), sg.Button("Check table", key = 'Taxon_table_check', button_color=('black', 'white'))],
					[sg.Text('Read table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'Read_table_path'), sg.Button("Check table", key = 'Read_table_check', button_color=('black', 'white'))],
                    [sg.Text('Name of TaXon table:', size=(20, 1)), sg.Input(project_folder + "_taxon_table", key = 'TaXon_table_name'), sg.Text('.xlsx')],
					[sg.Button("Convert", key = 'taXon_converter')],
					[sg.Text('_'*98)],
					[sg.Image(r'_source/ttt_title_image.png')],
					]

TaXon_table_processing_layout = [
					[sg.Text('',size=(1,1))],
					[sg.Text('Replicates processing', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'TaXon_table_processing_table')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('Combine replicates',size=(40,2), font=('Arial', 11, "bold"))],
                    [sg.Text("Replicate:"), sg.Radio('Extraction replicates (_a/_b)', 'Replicates', default=True, key='extraction_replicates'), sg.Radio('PCR replicates (_A/_B)', 'Replicates', key='PCR_replicates')],
					[sg.Text('',size=(1,1))],
					[sg.Text("Combine replicates", size=(30, 1)), sg.Button("Run", key='replicate_processing')],
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon table per sample', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text("Create taxon tables per-sample", size=(30,1)), sg.Button("Run", key = 'taxon_table_per_sample')],
					[sg.Text('',size=(1,1))],
					[sg.Text('Meta data table', size=(50,2), font=('Arial', 11, "bold")), sg.Text("", size=(1,1))],
					[sg.Text("Create meta data table", size=(30,1)), sg.Button("Run", key = 'run_meta_data'), sg.Text("", size=(5,1)),
                    sg.Frame(layout=[[sg.Text("Please adjust the meta data table to your needs.")]], title="Options")],
					[sg.Text('',size=(1,1))],
					]

TaXon_table_filtering_layout = [
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon filter', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_list_to_filter_path')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon-based filter',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Taxonomic level: "), sg.Radio('Phylum', "tt_filter", default=True, key='filter_phylum'), sg.Radio('Class', "tt_filter", key='filter_class'),
					sg.Radio('Order', "tt_filter", key='filter_order'), sg.Radio('Family', "tt_filter", key='filter_family'),
					sg.Radio('Genus', "tt_filter", key='filter_genus'), sg.Radio('Species', "tt_filter", key='filter_species')],
					[sg.Text("", size=(20, 1))],
					[sg.Text("Taxon based filter:", size=(25, 1)), sg.Button("Filter", key = 'taxon_filtering')],
					[sg.Text('',size=(1,1))],
					[sg.Text('Sample-based filter', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text("Sample based filter:", size=(25, 1)), sg.Button("Filter", key = 'sample_filtering')],
					[sg.Text('',size=(1,1))],
					[sg.Text('Replicate consistency filter',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Filter replicates", size=(25, 1)), sg.Button("Filter", key='replicate_consistency_filter'), sg.Text("", size=(3,1)),
                    sg.Frame(layout=[[sg.Text("Sample prompt for replicates is required.")]], title="Information")],
					[sg.Text('',size=(1,1))],
					]

statistics_1_layout = [
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon table analyses #1', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_1_taxon_table_path')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('Basic statistics',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Calculate basic statistics", size=(30,1)), sg.Button("Run", key = 'run_basic_stats')],
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxonomic resolution',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Create taxonomic resolution plot", size=(30,1)), sg.Button("Run", key = 'run_taxonomic_resolution'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("20", size=(3,1), key="x_tax_res"),
                    sg.Input("10", size=(3,1), key="y_tax_res"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("10", size=(3,1), key="font_tax_res")]], title="Options")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxonomic richness',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Create taxonomic richness plot", size=(30,1)), sg.Button("Run", key = 'run_taxonomic_richness'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("20", size=(3,1), key="x_tax_rich"),
                    sg.Input("10", size=(3,1), key="y_tax_rich"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("10", size=(3,1), key="font_tax_rich")]], title="Options")],
					[sg.Text('',size=(1,1))],
					]

statistics_2_layout =  [
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon table analyses #2', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_2_taxon_table_path')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('OTU abundance pie-charts',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Create pie-charts", size=(25,1)), sg.Button("Run", key = 'pie_chart'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Text("20", size=(3,1)),
                    sg.Text("10", size=(3,1)), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Text("10", size=(3,1))]], title="Options")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Venn diagrams',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text('Compare to table:', size=(16, 1)), sg.Input(), sg.FileBrowse(key = 'venn_taxon_table_2_path')],
					[sg.Text("", size=(1,1))],
					[sg.Text("Compare taxon tables", size=(25,1)), sg.Button("Run", key = 'venn_diagram'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Text("20", size=(3,1)),
                    sg.Text("10", size=(3,1)), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Text("10", size=(3,1))]], title="Options")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Rarefaction curve',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Calculate rarefaction curve with"), sg.Input('15000', size=(6,1), key='repetitions'), sg.Text('repetitions'), sg.Button("Run", key = 'rarefaction_curve'), sg.Text("", size=(1,1)),
                    sg.Frame(layout=[[sg.Text("Error bar style:"), sg.Radio('a', "error_bar_style", key='error_bar_style_a', default=True), sg.Radio('b', "error_bar_style", key='error_bar_style_b')]], title="Options")],
					[sg.Text('',size=(1,1))],
					]

statistics_3_layout =  [
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon table analyses #3', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_3_taxon_table_path')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('Read proportions',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Taxonomic level: "), sg.Radio('Phylum', "ra_filter", key='ra_phylum'), sg.Radio('Class', "ra_filter", key='ra_class'),
					sg.Radio('Order', "ra_filter", key='ra_order'), sg.Radio('Family', "ra_filter", key='ra_family'),
					sg.Radio('Genus', "ra_filter", key='ra_genus'), sg.Radio('Species', "ra_filter", default=True, key='ra_species')],
					[sg.Text("", size=(20, 1))],
					[sg.Text("Calculate read proportions", size=(25,1)), sg.Button("Run", key = 'run_read_proportions'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("10", size=(3,1), key="x_read_props_scatter"),
                    sg.Input("5", size=(3,1), key="y_read_props_scatter"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("4", size=(3,1), key="font_read_props_scatter")]], title="Options")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Site occupancy', size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Calculate site occupancy", size=(25,1)), sg.Button("Run", key = 'run_site_occupancy'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("10", size=(3,1), key="x_site_occ"),
                    sg.Input("5", size=(3,1), key="y_site_occ"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("8", size=(3,1), key="font_site_occ")]], title="Options")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Krona charts',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Create Krona chart", size=(25,1)), sg.Button("Run", key = 'create_krona_chart'), sg.Text("", size=(1,1)),
                    sg.Frame(layout=[[sg.Text("Charts are based on read proportions. Krona tools is required.")]], title="Information")],
					[sg.Text('',size=(1,1))],
					]

statistics_4_layout =  [
					[sg.Text('',size=(1,1))],
					[sg.Text('Taxon table analyses #4', size=(50,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_4_taxon_table_path')],
					[sg.Text('_'*98)],
					[sg.Text('',size=(1,1))],
					[sg.Text('PCoA analysis',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text("Calculate PCoA: "), sg.Button("Run", key = 'run_PCoA_analysis'), sg.Text("", size=(2,1)),
                    sg.Frame(layout=[[sg.Text("Scikit-bio dependecies are required.")]], title="Information")],
					[sg.Text('',size=(1,1))],
					]

create_taxon_list_layout = [
					[sg.Text('',size=(1,1))],
					[sg.Text('Create taxon list',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_list_taxon_table_path')],
					[sg.Text('_'*98)],
                    [sg.Text("Enter name of output file:"), sg.Input(project_folder + "_taxon_list", size=(48,1), key='taxon_list_output_file_name'), sg.Text(".txt + .xlsx")],
                    [sg.Text('Choose output language:'), sg.Radio('English', 'Language', key='English'), sg.Radio('German', 'Language', key='German', default=True)],
                    [sg.Text('', size=(1,1))],
                    [sg.Text('Additional information', font=('Arial', 10, "bold"))],
                    [sg.Text('Description:', size=(20,1)), sg.Input('', size=(48,1), key='TL_description')], # TL = taxon list -> makes it easier to find
                    [sg.Text('Author(s):', size=(20,1)), sg.Input('', size=(48,1), key='TL_authors')],
                    [sg.Text('Lab protocol:', size=(20,1)), sg.Input('', size=(48,2), key='TL_lab_protocol')],
                    [sg.Text('No. of replicates:', size=(20,1)), sg.Input('X (extraction); X (1-st-step); X (2nd-step)', size=(48,1), key='TL_replicates')],
                    [sg.Text('No. of neg. contr.:', size=(20,1)), sg.Input('', size=(48,1), key='TL_negative_controls')],
                    [sg.Text('Primers:', size=(20,1)), sg.Input('', size=(48,1), key='TL_primers')],
                    [sg.Text('Sequencing run:', size=(20,1)), sg.Input('', size=(48,1), key='TL_sequencing_run')],
                    [sg.Text('Bionf. pipeline:', size=(20,1)), sg.Input('', size=(48,1), key='TL_pipeline')],
                    [sg.Text('', size=(1,1))],
                    [sg.Text("Create taxon list"), sg.Button("Run", key = 'taxon_list')],
					]

create_asterics_input_layout = [
					[sg.Text('',size=(1,1))],
					[sg.Text('Perlodes',size=(40,2), font=('Arial', 11, "bold"))],
					[sg.Text('Taxon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'asterics_taxon_table_path')],
					[sg.Text('_'*98)],
					[sg.Text('Create Perlodes input file',size=(40,2), font=('Arial', 11, "bold"))],
#                    [sg.Frame(layout=[[sg.Text("Scikit-bio dependecies" + "\n" + "are required.", size=(92,5))]], title="About Perlodes")],
					[sg.Text('',size=(1,1))],
					[sg.Text('Operational taxon list:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'operational_taxon_list_path')],
					[sg.Text('',size=(1,1))],
					[sg.Text("Create Asterics input file: "), sg.Button("Run", key = 'create_asterics_input_file')],
					[sg.Button("More information", key = 'open_perlodes')],
                    ]

layout = [  [sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
			[sg.Text('',size=(1,1))],
			[sg.TabGroup([[
            sg.Tab('Data conversion', data_conversion),
            sg.Tab('Table processing', TaXon_table_processing_layout),
            sg.Tab('Filtering', TaXon_table_filtering_layout),
            sg.Tab('Analysis #1', statistics_1_layout),
            sg.Tab('Analysis #2', statistics_2_layout),
            sg.Tab('Analysis #3', statistics_3_layout),
            sg.Tab('Analysis #4', statistics_4_layout),
            sg.Tab('Taxon list', create_taxon_list_layout),
            sg.Tab('Asterics', create_asterics_input_layout)]])],
			[sg.Text('',size=(1,1))],
			[sg.Output(size=(100, 6))],
			[sg.Text('Progress:'), sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progbar')],
			[sg.Exit(button_color=('black', 'red')), sg.Text('', size=(65,1)), sg.Text(taxon_tools_version, font=('Arial', 8))]]


# Create the Window
window = sg.Window('TaXonTableTools', layout)
win2_active=False

##########################################################################################################################

while True:
    try:

        event, values = window.Read()

        i = 0
        # define variables
        TaXon_table_name = values['TaXon_table_name']
        TaXon_table_processing_table = values['TaXon_table_processing_table']
        extraction_replicates = values['extraction_replicates']
        PCR_replicates = values['PCR_replicates']
        taxonomy_table_path = values['taxonomy_table_path']
        Read_table_path = values['Read_table_path']
        venn_taxon_table_2_path = values['venn_taxon_table_2_path']
        taxon_list_taxon_table_path = values['taxon_list_taxon_table_path']
        taxon_list_to_filter_path = values['taxon_list_to_filter_path']
        filter_phylum = values['filter_phylum']
        filter_order = values['filter_order']
        filter_genus = values['filter_genus']
        filter_class = values['filter_class']
        filter_family = values['filter_family']
        filter_species = values['filter_species']
        ra_phylum = values['ra_phylum']
        ra_order = values['ra_order']
        ra_genus = values['ra_genus']
        ra_class = values['ra_class']
        ra_family = values['ra_family']
        ra_species = values['ra_species']
        statistics_1_taxon_table_path = values['statistics_1_taxon_table_path']
        statistics_2_taxon_table_path = values['statistics_2_taxon_table_path']
        statistics_3_taxon_table_path = values['statistics_3_taxon_table_path']
        statistics_4_taxon_table_path = values['statistics_4_taxon_table_path']
        taxon_list_output_file_name = values['taxon_list_output_file_name']
        English = values['English']
        German = values['German']
        repetitions = values['repetitions']
        x_tax_res = values["x_tax_res"]
        y_tax_res = values["y_tax_res"]
        font_tax_res = values["font_tax_res"]
        x_tax_rich = values["x_tax_rich"]
        y_tax_rich = values["y_tax_rich"]
        font_tax_rich = values["font_tax_rich"]
        x_read_props_scatter = values["x_read_props_scatter"]
        y_read_props_scatter = values["y_read_props_scatter"]
        font_read_props_scatter = values["font_read_props_scatter"]
        x_site_occ = values["x_site_occ"]
        y_site_occ = values["y_site_occ"]
        font_site_occ = values["font_site_occ"]
        error_bar_style_a = values["error_bar_style_a"]
        error_bar_style_b = values["error_bar_style_b"]

        if event is None or event == 'Exit':
            break

        print("########", "\n")

        count=0
        window.Element('progbar').UpdateBar(count)

        if event == 'taXon_converter':
            print("Converting to TaXon table")

            if (taxonomy_table_path == '' or Read_table_path == ''):
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                TaXon_converter(Read_table_path, taxonomy_table_path, TaXon_table_name, path_to_outdirs)
                print("Taxon tables are found in" , path_to_outdirs, "TaXon_tables/")
                sg.Popup("Taxon tables are found in", path_to_outdirs, "TaXon_tables/", title="Finished", keep_on_top=True)

        if event == 'Taxon_table_check':
            print("Checking taxonomy table format")

            if (taxonomy_table_path == ''):
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                print(taxonomy_table_path)
                check_taxon_table_format(taxonomy_table_path)

        if event == 'Read_table_check':
            print("Checking read table format")

            if (Read_table_path == ''):
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                print(Read_table_path)
                check_read_table_format(Read_table_path)

        if event == "taxon_table_per_sample":
            print("Creating OTU table for each sample")

            if TaXon_table_processing_table == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                create_TaXon_table_per_sample(TaXon_table_processing_table, path_to_outdirs)
                print("\n" + "Taxon tables are found in", path_to_outdirs, "TaXon_tables_per_sample/")
                sg.Popup("Taxon tables are found in", path_to_outdirs, "TaXon_tables_per_sample/", title="Finished", keep_on_top=True)

        if event == 'pie_chart':
            print("Pie chart")

            if statistics_2_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                pie_chart(statistics_2_taxon_table_path, path_to_outdirs)

        if event == 'venn_diagram':
            print("Venn diagrams")

            if (statistics_2_taxon_table_path == '' or venn_taxon_table_2_path == ''):
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                venn_diagram(statistics_2_taxon_table_path, venn_taxon_table_2_path, path_to_outdirs)

        if event == 'rarefaction_curve':
            print("Rarefaction curve")

            if statistics_2_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                if error_bar_style_a == True:
                    error_bar_style = "a"
                if error_bar_style_b == True:
                    error_bar_style = "b"

                rarefaction(statistics_2_taxon_table_path, int(repetitions), error_bar_style, path_to_outdirs)

        if event == 'taxon_list':
            print("Creating taxon list")

            if English == True:
                language = "English"

            if German == True:
                language = "German"

            if taxon_list_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                print("Taxon list")
                create_taxon_list(taxon_list_taxon_table_path, taxon_list_output_file_name, language, values, taxon_tools_version, path_to_outdirs)
                print("\n" + "Taxon lists are found in", path_to_outdirs, "Taxon_lists")
                sg.Popup("Taxon lists are found in", path_to_outdirs, "Taxon_lists", title="Finished", keep_on_top=True)

        if event == 'taxon_filtering' and not win2_active:
            print("Taxon filtering")
            if filter_phylum == True:
                mask = "Phylum"
            elif filter_class == True:
                mask = "Class"
            elif filter_order == True:
                mask = "Order"
            elif filter_family == True:
                mask = "Family"
            elif filter_genus == True:
                mask = "Genus"
            elif filter_species == True:
                mask = "Species"

            if taxon_list_to_filter_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")

            else:
                available_taxa = create_taxa_mask(taxon_list_to_filter_path, mask)
                win2_active = True
                window.Hide()
                available_taxa_list = list(slices([sg.CB(name, default=True) for name in sorted(available_taxa)], 5))
                layout2 = [[sg.Text("Taxon filter", size=(20,1))],
                [sg.Frame(layout = available_taxa_list, title = 'Check taxa to exlude')],
                [sg.Text('Enter appendix for output file: '), sg.InputText('file_appendix', key = 'appendix_name', size=(35,1))],
                [sg.Text("Filter taxa below a similarity of: "), sg.Input('0', size=(4,1), key='threshold'), sg.Text("%")],
                [sg.Button('Filter')],
                [sg.Button('Back')]]

                win2 = sg.Window('Taxon table filtering', layout2, keep_on_top=True)
                while True:
                    event2, values2 = win2.Read()

                    if event2 == 'Filter':

                        i = 0
                        check_dict = {}
                        check_list = []
                        appendix_name = values2['appendix_name']
                        threshold = values2['threshold']
                        del values2['appendix_name']
                        del values2['threshold']
                        for key, value in values2.items():
                            check_dict[available_taxa[i]] = value
                            i += 1
                        for key, value in check_dict.items():
                            if value == True:
                                check_list.append(key)

                        if check_list == []:
                            check_list = ['no taxa to filter']
                            taxon_table_filter(taxon_list_to_filter_path, check_list, mask, appendix_name, threshold, path_to_outdirs)

                        else:
                            taxon_table_filter(taxon_list_to_filter_path, check_list, mask, appendix_name, threshold, path_to_outdirs)

                    if event2 is None or event2 == 'Back':
                        win2.Close()
                        win2_active = False
                        window.UnHide()
                        break

        if event == 'sample_filtering' and not win2_active:
            print("Sample filtering")

            if taxon_list_to_filter_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                available_samples = sorted(create_sample_mask(taxon_list_to_filter_path))
                win2_active = True
                window.Hide()
                available_samples_list = list(slices([sg.CB(name, default=False) for name in sorted(available_samples)], 3))
                layout2 = [[sg.Text("Sample filter", size=(20,1))],
                [sg.Frame(layout = available_samples_list, title = 'Check samples to exlude')],
                [sg.Text('Enter appendix for output file: '), sg.InputText('file_appendix', key = 'appendix_name', size=(35,1))],
                [sg.Button('Filter')],
                [sg.Button('Back')]]

                win2 = sg.Window('Taxon table filtering', layout2, keep_on_top=True)

                while True:
                    event2, values2 = win2.Read()

                    if event2 == 'Filter':
                        i = 0
                        check_dict = {}
                        check_list = []
                        appendix_name = values2['appendix_name']
                        del values2['appendix_name']
                        for key, value in values2.items():
                            check_dict[available_samples[i]] = value
                            i += 1
                        for key, value in check_dict.items():
                            if value == True:
                                check_list.append(key)

                        print("Excluded samples:", ', '.join(check_list))
                        exclude_samples(taxon_list_to_filter_path, check_list, appendix_name, path_to_outdirs)

                    if event2 is None or event2 == 'Back':
                        win2.Close()
                        win2_active = False
                        window.UnHide()
                        break

        if event == 'replicate_processing':
            print("Removing duplicates")

            if PCR_replicates == True:
                replicate_version = "PCR"
            elif extraction_replicates == True:
                replicate_version = "extraction"

            if TaXon_table_processing_table == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                combine_replicates(TaXon_table_processing_table, replicate_version, path_to_outdirs)

        if event == 'run_taxonomic_richness':
            print("Taxon statistics")

            if statistics_1_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                calculate_taxonomic_richness(statistics_1_taxon_table_path, path_to_outdirs, x_tax_rich, y_tax_rich, font_tax_rich)

        if event == 'run_taxonomic_resolution':
            print("OTU statistics")

            if statistics_1_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                calculate_taxonomic_resolution(statistics_1_taxon_table_path, path_to_outdirs, x_tax_res, y_tax_res, font_tax_res)

        if event == 'run_meta_data':
            print("Creating meta data table")

            if TaXon_table_processing_table == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                create_metadata_table(TaXon_table_processing_table, path_to_outdirs)
                print("\n" + "Meta data tables are found in", path_to_outdirs, "Meta_data_table/")
                sg.Popup("Meta data tables are found in", path_to_outdirs, "Meta_data_table/", title="Finished", keep_on_top=True)

        if event == 'run_site_occupancy':
            print("Calculate site occupancy")

            if statistics_3_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")

            else:
                meta_data_to_test = get_available_meta_data(statistics_3_taxon_table_path, path_to_outdirs)

                if meta_data_to_test == False:
                             sg.PopupError("Meta data file not found!", keep_on_top=True)
                             print("Meta data file not found!")
                             print("\n" + "########")
                else:
                             win2_active = True
                             window.Hide()
                             meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))
                             taxonmic_levels_to_choose = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                             taxonmic_levels_list = list(slices([sg.Radio(name, "taxa", key=name, default=True) for name in taxonmic_levels_to_choose], 6))
                             layout2 = [[sg.Text("Site occupancy", size=(20,1))],
                             [sg.Frame(layout = meta_data_to_test_list, title = 'Check taxa to exlude')],
                             [sg.Text('',size=(1,1))],
                             [sg.Frame(layout = taxonmic_levels_list, title = 'Check taxonomic level to use')],
                             [sg.Text('',size=(1,1))],
                             [sg.Button('Calculate')],
                             [sg.Button('Back')]] # test = True

                             win2 = sg.Window('Site occpancy', layout2, keep_on_top=True)

                             while True:
                                          event2, values2 = win2.Read()

                                          if event2 == 'Calculate':
                                                       for input_value, test in values2.items():
                                                                    if (test == True and input_value not in taxonmic_levels_to_choose):
                                                                                 meta_data_to_test = input_value
                                                                    if (test == True and input_value not in meta_data_to_test):
                                                                                 taxon_level = input_value

                                                       site_occupancy(statistics_3_taxon_table_path, meta_data_to_test, taxon_level, path_to_outdirs, x_site_occ, y_site_occ, font_site_occ)
                                                       win2.Close()
                                                       win2_active = False
                                                       window.UnHide()

                                          if event2 is None or event2 == 'Back':
                                                       win2.Close()
                                                       win2_active = False
                                                       window.UnHide()
                                                       break

                             print("\n" + "########")

        if event == 'run_read_proportions':

            print("Creating read abundance plots")

            if statistics_3_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")

            else:
                if ra_phylum == True:
                    taxonomic_level = "Phylum"
                elif ra_class == True:
                    taxonomic_level = "Class"
                elif ra_order == True:
                    taxonomic_level = "Order"
                elif ra_family == True:
                    taxonomic_level = "Family"
                elif ra_genus == True:
                    taxonomic_level = "Genus"
                elif ra_species == True:
                    taxonomic_level = "Species"

                read_proportions(statistics_3_taxon_table_path, taxonomic_level, path_to_outdirs, x_read_props_scatter, y_read_props_scatter, font_read_props_scatter)

        if event == 'run_basic_stats':
            print("Basic stats")

            if statistics_1_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                basic_stats(statistics_1_taxon_table_path, path_to_outdirs)
                print("\n" + "Basic stats are found in", path_to_outdirs, "Basic_stats/")
                sg.Popup("Basic stats are found in", path_to_outdirs, "Basic_stats/", title="Finished", keep_on_top=True)

        if event == 'create_krona_chart':
            print("Krona chart")

            if statistics_3_taxon_table_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                create_krona_chart(statistics_3_taxon_table_path, path_to_outdirs)

        if event == 'run_PCoA_analysis':
                print("Creating meta data table")

                # try to import the required dependecies
                try:
                    import_abort = "False"
                    from scipy.spatial import distance
                    from skbio.diversity import alpha_diversity
                    from skbio.diversity import beta_diversity
                    from skbio.stats.distance import mantel
                    from skbio.stats.distance import anosim
                    from skbio.stats.ordination import pcoa
                    from skbio import TreeNode
                except ModuleNotFoundError:
                    pass
                    import_abort = "True"
                except ImportError:
                    pass
                    import_abort = "True"

                if statistics_4_taxon_table_path == '':
                    sg.PopupError("Please provide a file")
                    print("Error: Please provide a file")

                elif import_abort == "True":
                    sg.PopupError("You need to install skbio dependecies")
                    print("You need to install skbio dependecies")

                else:
                    meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Meta data file not found!")
                            print("Meta data file not found!")
                            print("\n" + "########")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                            layout2 = [[sg.Text("PCoA analysis", size=(20,1))],
                            [sg.Frame(layout = meta_data_to_test_list, title = 'Check taxa to exlude')],
                            [sg.Button('Calculate')],
                            [sg.Button('Back')]] # test = True

                            win2 = sg.Window('PCoA analysis', layout2, keep_on_top=True)

                            while True:

                                event2, values2 = win2.Read()

                                if event2 == 'Calculate':
                                        for input_value, test in values2.items():
                                                if test == True:

                                                        meta_data_to_test = input_value
                                                        PCoA_analysis(statistics_4_taxon_table_path, meta_data_to_test, path_to_outdirs)
                                                        win2.Close()
                                                        win2_active = False
                                                        window.UnHide()
                                                        break

                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

                    print("\n" + "########")

        if event == 'replicate_consistency_filter':
            print("Creating meta data table")

            if taxon_list_to_filter_path == '':
                sg.PopupError("Please provide a file", keep_on_top=True)
                print("Error: Please provide a file")
            else:
                replicate_consistency_filter(taxon_list_to_filter_path, path_to_outdirs)

        if event == 'open_perlodes':
            webbrowser.open('https://www.gewaesser-bewertung-berechnung.de')

        # close all windows to prevent crazy overlaps
        if plt.get_fignums():
            plt.close()

        print("\n" + "########")

    # end the script when KeyboardInterrupt is raised
    except KeyboardInterrupt:
        exit()

    except FilesDoNoMatch:
        print("Error!", "\n" + "    The OTU IDs do not match!", "\n" + "    >> Please check your files!")
        sg.PopupError("Error!", "\n" + "    The OTU IDs do not match!", "\n" + "    >> Please check your files!", keep_on_top=True)
        print("\n" + "########")

    except TaxonomyTableError:
        print("\n" + "########")

    except ReadTableError:
        pass
        print("\n" + "########")

    except KronaInstallError:
        pass
        print("\n" + "########")

    except PCoAImportError:
        pass
        print("\n" + "########")

    # if there are unexpected errors print a message and continue the script!
    except:
        #raise
        layout = [
                    [sg.Image(r'_source/crash.png'), sg.Text(" You've been awarded with the gold medal in program crashing!")],
                    [sg.Text("", size=(1,1))],
                    [sg.Text("Unexpected error: " + str(sys.exc_info()[0]))],
                    [sg.Text("", size=(1,1))],
                    [sg.Text("An unexpected error occured:")],
                    [sg.Text("Please refer to the manual.")],
                    [sg.Text("", size=(1,1))],
                    [sg.Button('Return', button_color=('black', 'red'))]
                 ]

        window_error = sg.Window('Error', layout, keep_on_top=True)

        while (True):

            event, values = window_error.Read()

            if event in ('Return', None):
                break

        window_error.Close()

        print("\n" + "########")

window.Close()
