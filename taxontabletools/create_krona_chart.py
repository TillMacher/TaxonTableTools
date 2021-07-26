def create_krona_chart_single(TaXon_table_xlsx, path_to_outdirs):

    import subprocess, os, webbrowser
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    try:
        subprocess.call(["ktImportText"], stdout=open(os.devnull, 'wb'))
    except:
        sg.PopupError("Krona tools must be manually installed first!" + "\n" * 2 + "Note: Krona tools is currently not supported on Windows!" + "\n", title="Error")
        raise RuntimeError("Krona tools needs to be installed")

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    TaXon_table_df = TaXon_table_df.replace(np.nan, '__', regex=True)

    ## create an output folder
    krona_chart_name = Path(TaXon_table_xlsx).name.replace(".xlsx", "")
    dirName = Path(str(path_to_outdirs) + "/Krona_charts/" + krona_chart_name)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        pa_data = True
    else:
        pa_data = False

    row1 = ["sample-ID", "", "", "", "", "", ""]
    row2 = ["count", "phylum", "class", "order", "family", "genus", "species"]
    krona_taxonomy_list = []
    krona_taxonomy_list.append(row1)
    krona_taxonomy_list.append(row2)

    for OTU in TaXon_table_df.values.tolist():
        taxonomy = OTU[1:7]
        reads = sum(OTU[10:])
        if pa_data == True:
            krona_taxonomy_list.append([1] + taxonomy)
        else:
            krona_taxonomy_list.append([reads] + taxonomy)

    krona_taxonomy_df = pd.DataFrame(krona_taxonomy_list)

    krona_chart_directory = Path(str(path_to_outdirs) + "/" + "Krona_charts" + "/" + TaXon_table_xlsx.stem)
    krona_table_tsv = Path(str(dirName) + "/single_krona_table.tsv")
    krona_chart_html = Path(str(dirName) + "_krona_single.html")

    # write krona table to tsv
    krona_taxonomy_df.to_csv(krona_table_tsv, sep="\t", header=False, index=False)

    os.system("ktImportText " + str(krona_table_tsv) + " -o " + str(krona_chart_html))

    # finish script
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(krona_chart_html))

    closing_text = "Krona chart is found under:\n" + '/'.join(str(krona_chart_html).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("krona chart", "analysis", TaXon_table_xlsx.name, krona_chart_html.name, "nan", path_to_outdirs)

def create_krona_chart_multi(TaXon_table_xlsx, path_to_outdirs):

    import subprocess, os, webbrowser
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    try:
        subprocess.call(["ktImportText"], stdout=open(os.devnull, 'wb'))
    except:
        sg.PopupError("Krona tools must be manually installed first!" + "\n" * 2 + "Note: Krona tools is currently not supported on Windows!" + "\n", title="Error")
        raise RuntimeError("Krona tools needs to be installed")

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    TaXon_table_df = TaXon_table_df.replace(np.nan, '__', regex=True)

    samples = TaXon_table_df.columns.tolist()[10:]
    columns = TaXon_table_df.columns.tolist()[:10]

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test == {1,0}:
        pa_data = True
    else:
        pa_data = False

    ## create an output folder
    krona_chart_name = Path(TaXon_table_xlsx).name.replace(".xlsx", "")
    dirName = Path(str(path_to_outdirs) + "/Krona_charts/" + krona_chart_name)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    ## store the names of the sample tsv files
    sample_tsv_path = []

    ## write a seperate tsv file for each sample in the TaXon table
    for sample in samples:
        row1 = ["sample-ID", "", "", "", "", "", ""]
        row2 = ["count", "phylum", "class", "order", "family", "genus", "species"]
        krona_taxonomy_list = []
        krona_taxonomy_list.append(row1)
        krona_taxonomy_list.append(row2)

        for OTU in TaXon_table_df[columns+[sample]].values.tolist():
            taxonomy = OTU[1:7]
            reads = sum(OTU[10:])
            if reads != 0:
                if pa_data == True:
                    krona_taxonomy_list.append([1] + taxonomy)
                else:
                    krona_taxonomy_list.append([reads] + taxonomy)

        ## store the data in df
        krona_taxonomy_df = pd.DataFrame(krona_taxonomy_list)
        krona_table_tsv = Path(str(dirName) + "/" + sample.replace(" ", "_") + "_krona_table.tsv")
        sample_tsv_path.append(str(krona_table_tsv))
        # write krona table to tsv
        krona_taxonomy_df.to_csv(krona_table_tsv, sep="\t", header=False, index=False)

    krona_chart_html = Path(str(dirName) + "_krona_multi.html")
    os.system("ktImportText " + ' '.join(sample_tsv_path) + " -o " + str(krona_chart_html))

    # finish script
    answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
    if answer == "Yes":
        webbrowser.open('file://' + str(krona_chart_html))

    closing_text = "Krona chart is found under:\n" + '/'.join(str(krona_chart_html).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("krona chart", "analysis", TaXon_table_xlsx.name, krona_chart_html.name, "nan", path_to_outdirs)
