def beta_diversity(TaXon_table_xlsx, beta_w, beta_h, beta_cmap, meta_data_to_test, path_to_outdirs, template):

    import pandas as pd
    import numpy as np
    from skbio.diversity import beta_diversity
    from skbio.stats.distance import anosim
    import plotly.express as px
    from pathlib import Path
    import PySimpleGUI as sg
    import webbrowser

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    Meta_data_table_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    TaXon_table_samples = TaXon_table_df.columns.tolist()[10:]
    Meta_data_table_df = pd.read_excel(Meta_data_table_xlsx, header=0)
    Meta_data_table_samples = Meta_data_table_df['Samples'].tolist()
    metadata_list = Meta_data_table_df[meta_data_to_test].values.tolist()

    # check for presence absence data
    # otherwise abort and print error message
    pa_test = set([val for sublist in TaXon_table_df[TaXon_table_samples].values.tolist() for val in sublist])
    if pa_test != {1,0}:
        sg.Popup("Please use presence absence data!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == len(Meta_data_table_df['Samples'].tolist()):
        sg.Popup("The meta data is unique for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    # check if the meta data differs
    if len(set(Meta_data_table_df[meta_data_to_test])) == 1:
        sg.Popup("The meta data is similar for all samples. Please adjust the meta data table!", title=("Error"))
        raise RuntimeError

    if sorted(TaXon_table_samples) == sorted(Meta_data_table_samples):

        samples = Meta_data_table_samples
        data = TaXon_table_df[TaXon_table_samples].transpose().values.tolist()
        jaccard_dm = beta_diversity("jaccard", data, samples)

        anosim_results = anosim(jaccard_dm, metadata_list, permutations=999)
        anosim_r = round(anosim_results['test statistic'], 5)
        anosim_p = anosim_results['p-value']
        textbox = "Anosim (" + meta_data_to_test + ")\n" + "R = " + str(anosim_r) + "\n" + "p = " + str(anosim_p)

        matrix = jaccard_dm.data
        matrix_df = pd.DataFrame(matrix)
        matrix_df.columns = samples
        matrix_df.index = samples

        # create plot
        fig = px.imshow(matrix, x=samples,y=samples, labels=dict(color="Jaccard distance"))
        fig.update_layout(height=int(beta_h), width=int(beta_w), template=template, showlegend=True, title=textbox)

        # finish script
        output_pdf = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_jc.pdf")
        output_html = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_jc.html")
        output_xlsx = Path(str(path_to_outdirs) + "/" + "Beta_diversity" + "/" + TaXon_table_xlsx.stem + "_" + meta_data_to_test + "_jc.xlsx")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))
        matrix_df.to_excel(output_xlsx)

        ## ask to show plot
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## write to log file
        sg.Popup("Beta diversity estimate are found in", path_to_outdirs, "/Beta_diversity/", title="Finished", keep_on_top=True)
        from taxontabletools.create_log import ttt_log
        ttt_log("beta diversity", "analysis", TaXon_table_xlsx.name, output_pdf.name, meta_data_to_test, path_to_outdirs)

    else:
        sg.PopupError("Error: The samples between the taxon table and meta table do not match!", keep_on_top=True)
