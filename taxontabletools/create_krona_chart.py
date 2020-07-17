# krona chart
def create_krona_chart(TaXon_table_xlsx, path_to_outdirs):

    import subprocess, os
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
    krona_table_tsv = Path(str(krona_chart_directory) + "_krona_table.tsv")
    krona_chart_html = Path(str(krona_chart_directory) + "_krona_table.html")

    # write krona table to tsv
    krona_taxonomy_df.to_csv(krona_table_tsv, sep="\t", header=False, index=False)

    os.system("ktImportText " + str(krona_table_tsv) + " -o " + str(krona_chart_html))

    closing_text = "Krona chart is found under:\n" + '/'.join(str(krona_chart_html).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)
