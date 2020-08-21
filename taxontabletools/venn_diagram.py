def venn_diagram(file_a, file_b, file_c, venn_diagram_name, path_to_outdirs):

    import os
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2
    from matplotlib_venn import venn3
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path

    file_a = Path(file_a)
    file_b = Path(file_b)

    if file_c == False:
    ############################################################################
    # use venn2

        print("\n" + "Input file a:", file_a.stem)
        print("Input file b:", file_b.stem, "\n")
        count = 0

        allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species"]

        venn_dict = {}

        ############################################################################
        ## create the progress bar window
        layout = [[sg.Text('Progress bar')],
                  [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
                  [sg.Cancel()]]
        window_progress_bar = sg.Window('Progress bar', layout)
        progress_bar = window_progress_bar['progressbar']
        progress_update = 167*2
        ############################################################################

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

            dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
            if not os.path.exists(dirName):
                os.mkdir(dirName)

            output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
            plt.title(taxon)
            plt.savefig(output_pdf, bbox_inches='tight')

            if taxon == "Species":
                plt.show(block=False)
                print("\n" + "Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/")
                sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

            plt.close()

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                print('Cancel')
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += 167
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
        df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
        df.to_excel(output_xlsx, index=False)

        from taxontabletools.create_log import ttt_log
        ttt_log("venn diagram", "analysis", file_a.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_b.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)

    else:
    ############################################################################
    # use venn3

        if file_c == '':
            sg.PopupError("Please provide a file", keep_on_top=True)
            print("Error: Please provide a file")
            raise RuntimeError()

        file_c = Path(file_c)

        print("\n" + "Input file a:", file_a.stem)
        print("Input file b:", file_b.stem, "\n")
        print("Input file c:", file_c.stem, "\n")

        count = 0

        allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus","F_Species"]

        venn_dict = {}

        ############################################################################
        ## create the progress bar window
        layout = [[sg.Text('Progress bar')],
                  [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
                  [sg.Cancel()]]
        window_progress_bar = sg.Window('Progress bar', layout)
        progress_bar = window_progress_bar['progressbar']
        progress_update = 167*2
        ############################################################################

        for taxon in allowed_taxa:

            output_name = taxon
            taxon = taxon[2:]

            data_file_a = pd.read_excel(file_a, 'TaXon table', header=0)
            data_file_b = pd.read_excel(file_b, 'TaXon table', header=0)
            data_file_c = pd.read_excel(file_c, 'TaXon table', header=0)

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

            venn_dict[taxon + "_a_only"] = a_only
            venn_dict[taxon + "_b_only"] = b_only
            venn_dict[taxon + "_c_only"] = c_only
            venn_dict[taxon + "_shared_all"] = shared_all
            venn_dict[taxon + "_shared_a_b"] = shared_a_b
            venn_dict[taxon + "_shared_a_c"] = shared_a_c
            venn_dict[taxon + "_shared_b_c"] = shared_b_c

            print("Comparing on:", taxon, "level")

            plt.figure(figsize=(20, 10))
            venn3(subsets = (len_a_only, len_b_only, len_shared_a_b, len_c_only, len_shared_a_c, len_shared_b_c, len_shared_all), set_labels = (file_name_a, file_name_b, file_name_c))

            dirName = Path(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name)
            if not os.path.exists(dirName):
                os.mkdir(dirName)

            output_pdf = Path(str(dirName) + "/" + output_name + ".pdf")
            plt.title(taxon)
            plt.savefig(output_pdf, bbox_inches='tight')

            if taxon == "Species":
                plt.show(block=False)
                print("\n" + "Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/")
                sg.Popup("Venn diagrams are found in", path_to_outdirs, "Venn_diagrams/", title="Finished", keep_on_top=True)

            plt.close()

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                print('Cancel')
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += 167
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        output_xlsx = Path(str(dirName) + "/" + "Venn_comparison_results.xlsx")
        df = pd.DataFrame.from_dict(venn_dict, orient='index').transpose()
        df.to_excel(output_xlsx, index=False)

        from taxontabletools.create_log import ttt_log
        ttt_log("venn diagram", "analysis", file_a.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_b.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)
        ttt_log("venn diagram", "analysis", file_c.name, output_xlsx.name, venn_diagram_name, path_to_outdirs)

    ############################################################################
    # finish script
    closing_text = "Venn diagram is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
    print(closing_text)
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
