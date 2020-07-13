def pie_charts(TaXon_table_xlsx, path_to_outdirs, pc_label_font_size):

    ## import functions
    import os
    import PySimpleGUI as sg
    from pathlib import Path
    import pandas as pd
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 167*2
    ############################################################################

    TaXon_table_file =  Path(TaXon_table_xlsx)
    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_file)
    allowed_taxa = ["A_Phylum","B_Class","C_Order","D_Family","E_Genus", "F_Species"]

    print("\n" + "Input file:", TaXon_table_file.name)

    for taxon_query in allowed_taxa:
        taxon_level = taxon_query
        taxon_query = taxon_query[2:]

        print("Calclulating pie chart for:", taxon_query, "level")

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

            patches, texts, autotexts = plt.pie(taxa_sizes, startangle=140, autopct='%1.1f%%', wedgeprops=dict(width=0.5))

            bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.12)
            kw = dict(arrowprops=dict(arrowstyle='-', linewidth=0.1), bbox=bbox_props, zorder=0, va="center")

            for i, p in enumerate(patches):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = np.sin(np.deg2rad(ang))
                x = np.cos(np.deg2rad(ang))
                horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
                connectionstyle = "angle,angleA=0,angleB={}".format(ang)
                kw["arrowprops"].update({"connectionstyle": connectionstyle})
                plt.annotate(taxa_labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y), horizontalalignment=horizontalalignment, **kw, fontsize=pc_label_font_size)

            for i in range(len(texts)):
                texts[i].set_fontsize(4)

            for i in range(len(autotexts)):
                autotexts[i].set_fontsize(4)

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

            patches, texts, autotexts = plt.pie(taxa_sizes, startangle=140, autopct='%1.1f%%', wedgeprops=dict(width=0.5))

            bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.12)
            kw = dict(arrowprops=dict(arrowstyle='-', linewidth=0.1), bbox=bbox_props, zorder=0, va="center")

            for i, p in enumerate(patches):
                ang = (p.theta2 - p.theta1)/2. + p.theta1
                y = np.sin(np.deg2rad(ang))
                x = np.cos(np.deg2rad(ang))
                horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
                connectionstyle = "angle,angleA=0,angleB={}".format(ang)
                kw["arrowprops"].update({"connectionstyle": connectionstyle})
                plt.annotate(taxa_labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y), horizontalalignment=horizontalalignment, **kw, fontsize=pc_label_font_size)

            for i in range(len(texts)):
                texts[i].set_fontsize(4)

            for i in range(len(autotexts)):
                autotexts[i].set_fontsize(4)

            plt.axis('equal')
            plt.title("OTU abundance (excluding nan) on " + taxon_query + " level", fontsize=10, pad=20)

            pdf.savefig()

            if taxon_query == "Species":
                plt.show(block=False)
                closing_text = "\n" + "Pie charts are found in: " + str(path_to_outdirs) + "/Pie_charts/"
                print(closing_text)
                sg.Popup(closing_text, title="Finished", keep_on_top=True)

            plt.close()

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            print('Cancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_bar.UpdateBar(progress_update)
        progress_update += 167
        ############################################################################
    window_progress_bar.close()
