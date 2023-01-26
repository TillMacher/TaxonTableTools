import os
import time
import PySimpleGUI as sg
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3
from matplotlib.pyplot import plot, ion, show
from pathlib import Path
from taxontabletools.taxontable_manipulation import strip_metadata
from taxontabletools.create_log import ttt_log

file_a = '/Users/tillmacher/Desktop/TTT_projects/Projects/Lippe_eRNA_fwh_invertebrates/TaXon_tables/Lippe_eRNA_fwh_taxon_table_cons_NCsub_invertebrates_normalized_eDNA.xlsx'
file_b = '/Users/tillmacher/Desktop/TTT_projects/Projects/Lippe_eRNA_fwh_invertebrates/TaXon_tables/Lippe_eRNA_fwh_taxon_table_cons_NCsub_invertebrates_normalized_eRNA.xlsx'
file_c = ''
venn_diagram_name = 'test'
path_to_outdirs = 'test'
clustering_unit = 'OTUs'

def venn_diagram(file_a, file_b, file_c, venn_diagram_name, path_to_outdirs, clustering_unit):

    file_a = Path(file_a)
    file_b = Path(file_b)
    venn_font = 40

    if file_c == False:
    ############################################################################
    # use venn2

        count = 0

        G = "G_" + clustering_unit
        allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species", G]
        venn_dict = {}
        data_file_a = pd.read_excel(file_a, header=0).fillna('')
        data_file_b = pd.read_excel(file_b, header=0).fillna('')

        for taxon in allowed_taxa:

            time.sleep(1)

            output_name = taxon
            taxon = taxon[2:]
            col_name = taxon

            if taxon in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
                col_name = taxon
                taxon="ID"

            file_name_a = file_a.stem
            file_name_b = file_b.stem

            taxa_file_a = set([i for i in data_file_a[taxon].values.tolist() if i != ''])
            taxa_file_b = set([i for i in data_file_b[taxon].values.tolist() if i != ''])

            a_only = taxa_file_a - taxa_file_b
            len_a_only = len(a_only)
            b_only = taxa_file_b - taxa_file_a
            len_b_only = len(b_only)
            shared = taxa_file_a & taxa_file_b
            len_shared = len(shared)

            venn_dict[col_name + "_a_only"] = a_only
            venn_dict[col_name + "_shared"] = shared
            venn_dict[col_name + "_b_only"] = b_only

            plt.figure(figsize=(20, 10))
            out = venn2(subsets = (len_a_only, len_b_only, len_shared), set_labels = (file_name_a, file_name_b))
            for text in out.set_labels:
                text.set_fontsize(venn_font)
            for x in range(len(out.subset_labels)):
                if out.subset_labels[x] is not None:
                    out.subset_labels[x].set_fontsize(venn_font)

            dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
            if not os.path.exists(dirName):
                os.mkdir(dirName)

            output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
            plt.title(output_name[2:])
            plt.savefig(output_pdf, bbox_inches='tight')
            plt.close()

        output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
        df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
        df.to_excel(output_xlsx, index=False)

        sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

        ttt_log("venn diagram", "analysis", file_a.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_b.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)

    else:
        ############################################################################
        # use venn3

        if file_c == '':
            sg.PopupError("Please provide a file", keep_on_top=True)
            raise RuntimeError()

        file_c = Path(file_c)

        count = 0

        G = "G_" + clustering_unit
        allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species", G]

        venn_dict = {}

        data_file_a = pd.read_excel(file_a, header=0).fillna('')
        data_file_b = pd.read_excel(file_b, header=0).fillna('')
        data_file_c = pd.read_excel(file_c, header=0).fillna('')

        for taxon in allowed_taxa:

            time.sleep(1)

            output_name = taxon
            taxon = taxon[2:]
            col_name = taxon

            if taxon in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
                col_name = taxon
                taxon="ID"

            file_name_a = file_a.stem
            file_name_b = file_b.stem
            file_name_c = file_c.stem

            taxa_file_a = data_file_a[taxon].values.tolist()
            taxa_file_b = data_file_b[taxon].values.tolist()
            taxa_file_c = data_file_c[taxon].values.tolist()

            taxa_unique_a = list(dict.fromkeys(taxa_file_a))
            taxa_unique_b = list(dict.fromkeys(taxa_file_b))
            taxa_unique_c = list(dict.fromkeys(taxa_file_c))

            taxa_labels_a = []
            taxa_labels_b = []
            taxa_labels_c = []
            taxa_sizes_a = []
            taxa_sizes_b = []
            taxa_sizes_c = []

            for taxon_name in taxa_unique_a:
                if "nan" != str(taxon_name):
                    taxa_labels_a.append(str(taxon_name))
                    taxa_sizes_a.append(taxa_file_a.count(taxon_name))

            for taxon_name in taxa_unique_b:
                if "nan" != str(taxon_name):
                    taxa_labels_b.append(str(taxon_name))
                    taxa_sizes_b.append(taxa_file_b.count(taxon_name))

            for taxon_name in taxa_unique_c:
                if "nan" != str(taxon_name):
                    taxa_labels_c.append(str(taxon_name))
                    taxa_sizes_c.append(taxa_file_c.count(taxon_name))

            taxa_labels_a = sorted(taxa_labels_a)
            taxa_labels_b = sorted(taxa_labels_b)
            taxa_labels_c = sorted(taxa_labels_c)

            a_only = set(taxa_labels_a) - set(taxa_labels_b) - set(taxa_labels_c)
            len_a_only = len(a_only)
            b_only = set(taxa_labels_b) - set(taxa_labels_a) - set(taxa_labels_c)
            len_b_only = len(b_only)
            c_only = set(taxa_labels_c) - set(taxa_labels_a) - set(taxa_labels_b)
            len_c_only = len(c_only)

            shared_all = set(taxa_labels_a) & set(taxa_labels_b) & set(taxa_labels_c)
            len_shared_all = len(shared_all)
            shared_a_b = set(taxa_labels_a) & set(taxa_labels_b) - set(taxa_labels_c)
            len_shared_a_b = len(shared_a_b)
            shared_a_c = set(taxa_labels_a) & set(taxa_labels_c) - set(taxa_labels_b)
            len_shared_a_c = len(shared_a_c)
            shared_b_c = set(taxa_labels_b) & set(taxa_labels_c) - set(taxa_labels_a)
            len_shared_b_c = len(shared_b_c)

            venn_dict[col_name + "_a_only"] = a_only
            venn_dict[col_name + "_b_only"] = b_only
            venn_dict[col_name + "_c_only"] = c_only
            venn_dict[col_name + "_shared_all"] = shared_all
            venn_dict[col_name + "_shared_a_b"] = shared_a_b
            venn_dict[col_name + "_shared_a_c"] = shared_a_c
            venn_dict[col_name + "_shared_b_c"] = shared_b_c

            plt.figure(figsize=(20, 10))
            out = venn3(subsets = (len_a_only, len_b_only, len_shared_a_b, len_c_only, len_shared_a_c, len_shared_b_c, len_shared_all), set_labels = (file_name_a, file_name_b, file_name_c))
            for text in out.set_labels:
                text.set_fontsize(venn_font)
            for x in range(len(out.subset_labels)):
                if out.subset_labels[x] is not None:
                    out.subset_labels[x].set_fontsize(venn_font)

            dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
            if not os.path.exists(dirName):
                os.mkdir(dirName)

            output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
            plt.title(output_name[2:])
            plt.savefig(output_pdf, bbox_inches='tight')
            plt.close()

        output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
        df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
        df.to_excel(output_xlsx, index=False)

        sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

        ttt_log("venn diagram", "analysis", file_a.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_b.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_c.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)

def venn_diagram_metadata(TaXon_table_xlsx, venn_diagram_name, meta_data_to_test, path_to_outdirs, clustering_unit):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('unidentified')
    TaXon_table_df = strip_metadata(TaXon_table_df)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]

    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
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
    u_categories = list(set(categories))
    len_categories = len(u_categories)

    categories_dict = {}
    for sample, metadata in Meta_data_table_df[["Samples", meta_data_to_test]].values.tolist():
        if metadata in categories_dict.keys():
            categories_dict[metadata] = categories_dict[metadata] + [sample]
        else:
            categories_dict[metadata] = [sample]

    ## create a y axis title text
    taxon_title = clustering_unit.lower()

    ## adjust taxonomic level if neccessary
    if clustering_unit in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = clustering_unit
        clustering_unit = "ID"

    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.PopupError("The meta data has to differ between samples!", title="Error")
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):
        venn_font = 40

        if len_categories == 2:
        ############################################################################
        # use venn2

            count = 0

            G = "G_" + clustering_unit
            allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species", G]

            venn_dict = {}

            for taxon in allowed_taxa:

                time.sleep(1)

                output_name = taxon
                taxon = taxon[2:]
                col_name = taxon

                if taxon in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
                    col_name = taxon
                    taxon="ID"

                ## search present taxa for each metadata
                samples_a = categories_dict[u_categories[0]]
                category_a = u_categories[0]
                taxa_labels_a = sorted(list(set([i[0] for i in TaXon_table_df[[taxon] + samples_a].values.tolist() if i[0] != 'unidentified' and sum(i[1::]) != 0])))

                samples_b = categories_dict[u_categories[1]]
                taxa_labels_b = sorted(list(set([i[0] for i in TaXon_table_df[[taxon] + samples_b].values.tolist() if i[0] != 'unidentified' and sum(i[1::]) != 0])))
                category_b = u_categories[1]

                a_only = set(taxa_labels_a) - set(taxa_labels_b)
                len_a_only = len(a_only)
                b_only = set(taxa_labels_b) - set(taxa_labels_a)
                len_b_only = len(b_only)
                shared = set(taxa_labels_a) & set(taxa_labels_b)
                len_shared = len(shared)

                venn_dict[col_name + "_a_only"] = a_only
                venn_dict[col_name + "_shared"] = shared
                venn_dict[col_name + "_b_only"] = b_only

                plt.figure(figsize=(20, 10))
                out = venn2(subsets = (len_a_only, len_b_only, len_shared), set_labels = (category_a, category_b))
                for text in out.set_labels:
                    text.set_fontsize(venn_font)
                for x in range(len(out.subset_labels)):
                    if out.subset_labels[x] is not None:
                        out.subset_labels[x].set_fontsize(venn_font)

                dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
                if not os.path.exists(dirName):
                    os.mkdir(dirName)

                output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
                plt.title(output_name[2:])
                plt.savefig(output_pdf, bbox_inches='tight')

                plt.close()

            output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
            df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
            df.to_excel(output_xlsx, index=False)

            sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

            ttt_log("venn diagram", "analysis", category_a, output_xlsx.name, venn_diagram_name, path_to_outdirs)
            ttt_log("venn diagram", "analysis", category_b, output_xlsx.name, venn_diagram_name, path_to_outdirs)

        elif len_categories == 3:
        ############################################################################
        # use venn3

            count = 0

            G = "G_" + clustering_unit
            allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species", G]

            venn_dict = {}

            for taxon in allowed_taxa:

                time.sleep(1)

                output_name = taxon
                taxon = taxon[2:]
                col_name = taxon

                if taxon in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
                    col_name = taxon
                    taxon="ID"

                ## search present taxa for each metadata
                samples_a = categories_dict[u_categories[0]]
                category_a = u_categories[0]
                taxa_labels_a = sorted(list(set([i[0] for i in TaXon_table_df[[taxon] + samples_a].values.tolist() if i[0] != 'unidentified' and sum(i[1::]) != 0])))

                samples_b = categories_dict[u_categories[1]]
                category_b = u_categories[1]
                taxa_labels_b = sorted(list(set([i[0] for i in TaXon_table_df[[taxon] + samples_b].values.tolist() if i[0] != 'unidentified' and sum(i[1::]) != 0])))

                samples_c = categories_dict[u_categories[2]]
                category_c = u_categories[2]
                taxa_labels_c = sorted(list(set([i[0] for i in TaXon_table_df[[taxon] + samples_c].values.tolist() if i[0] != 'unidentified' and sum(i[1::]) != 0])))

                a_only = set(taxa_labels_a) - set(taxa_labels_b) - set(taxa_labels_c)
                len_a_only = len(a_only)
                b_only = set(taxa_labels_b) - set(taxa_labels_a) - set(taxa_labels_c)
                len_b_only = len(b_only)
                c_only = set(taxa_labels_c) - set(taxa_labels_a) - set(taxa_labels_b)
                len_c_only = len(c_only)

                shared_all = set(taxa_labels_a) & set(taxa_labels_b) & set(taxa_labels_c)
                len_shared_all = len(shared_all)
                shared_a_b = set(taxa_labels_a) & set(taxa_labels_b) - set(taxa_labels_c)
                len_shared_a_b = len(shared_a_b)
                shared_a_c = set(taxa_labels_a) & set(taxa_labels_c) - set(taxa_labels_b)
                len_shared_a_c = len(shared_a_c)
                shared_b_c = set(taxa_labels_b) & set(taxa_labels_c) - set(taxa_labels_a)
                len_shared_b_c = len(shared_b_c)

                venn_dict[col_name + "_a_only"] = a_only
                venn_dict[col_name + "_b_only"] = b_only
                venn_dict[col_name + "_c_only"] = c_only
                venn_dict[col_name + "_shared_all"] = shared_all
                venn_dict[col_name + "_shared_a_b"] = shared_a_b
                venn_dict[col_name + "_shared_a_c"] = shared_a_c
                venn_dict[col_name + "_shared_b_c"] = shared_b_c

                plt.figure(figsize=(20, 10))
                out = venn3(subsets = (len_a_only, len_b_only, len_shared_a_b, len_c_only, len_shared_a_c, len_shared_b_c, len_shared_all), set_labels = (category_a, category_b, category_c))
                for text in out.set_labels:
                    text.set_fontsize(venn_font)
                for x in range(len(out.subset_labels)):
                    if out.subset_labels[x] is not None:
                        out.subset_labels[x].set_fontsize(venn_font)

                dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
                if not os.path.exists(dirName):
                    os.mkdir(dirName)

                output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
                plt.title(output_name[2:])
                plt.savefig(output_pdf, bbox_inches='tight')

                plt.close()

            output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
            df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
            df.to_excel(output_xlsx, index=False)

            sg.Popup("Venn diagrams are found in", str(path_to_outdirs), "Venn_diagrams/", title="Finished", keep_on_top=True)

            ttt_log("venn diagram", "analysis", category_a, output_xlsx.name, venn_diagram_name, path_to_outdirs)
            ttt_log("venn diagram", "analysis", category_b, output_xlsx.name, venn_diagram_name, path_to_outdirs)
            ttt_log("venn diagram", "analysis", category_c, output_xlsx.name, venn_diagram_name, path_to_outdirs)

        else:
            sg.PopupError("You cannot use more than 3 metadata categories for Venn diagrams!")

    else:
        sg.PopupError("Error: The samples do not match", title="Error")
