def replicate_correlation_analysis(TaXon_table_xlsx, suffix_list, path_to_outdirs, width, height, repcorr_options, template, theme, font_size, clustering_unit):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import numpy as np
    import scipy.stats, webbrowser, os
    import matplotlib.pyplot as plt
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import statsmodels.api as sm

    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]

    width = int(width)
    height = int(height)
    x_zero = repcorr_options[0]
    y_zero = repcorr_options[1]

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)

    ## create an output folder
    replicate_analysis_name = Path(TaXon_table_xlsx).name.replace(".xlsx", "")
    dirName = Path(str(path_to_outdirs) + "/Replicate_analysis/" + replicate_analysis_name)
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    sample_names = TaXon_table_df.columns[10:].tolist()
    OTUs = TaXon_table_df["ID"].values.tolist()

    derep_sample_names_dict =  {}
    unique_sample_names_list = []
    replicates_dict = {}

    for sample in sample_names:
        sample_name = sample.split("_")[0:-1]
        unique_sample_names_list.append("_".join(sample_name))

    unique_sample_names_set = sorted(set(unique_sample_names_list))

    if len(suffix_list) == 2:

        ##############################
        # number of reads

        x1 = []
        y1 = []

        for sample in unique_sample_names_set:
            rep1, rep2 = [], []
            try:
                rep1 = (sample + "_" + suffix_list[0])
                rep2 = (sample + "_" + suffix_list[1])
                if (rep1 != [] and rep2 != []):
                    reads_rep1 = sum(TaXon_table_df[str(rep1)].values.tolist())
                    reads_rep2 = sum(TaXon_table_df[str(rep2)].values.tolist())
                    y1.append(reads_rep1)
                    x1.append(reads_rep2)
            except:
                print("No replicates found for:", sample)

        ##############################
        # number of OTUs

        x2 = []
        y2 = []

        for sample in unique_sample_names_set:
            rep1, rep2 = [], []
            try:
                rep1 = (sample + "_" + suffix_list[0])
                rep2 = (sample + "_" + suffix_list[1])
                if (rep1 != [] and rep2 != []):
                    OTUs_rep1 = len([OTU for OTU in TaXon_table_df[str(rep1)].values.tolist() if OTU != 0])
                    OTUs_rep2 = len([OTU for OTU in TaXon_table_df[str(rep2)].values.tolist() if OTU != 0])
                    y2.append(OTUs_rep1)
                    x2.append(OTUs_rep2)
            except:
                print("No replicates found for:", sample)

        # spearman's rho
        # reads
        spearman = scipy.stats.spearmanr(x1, y1)
        spearman_p_reads = spearman[1]
        if spearman_p_reads <= 0.05:
            spearman_rho_reads = str(round(spearman[0], 3)) + "*"
        else:
            spearman_rho_reads = str(round(spearman[0], 3))
        # OTUs
        spearman = scipy.stats.spearmanr(x2, y2)
        spearman_p_OTUs = spearman[1]
        if spearman_p_OTUs <= 0.05:
            spearman_rho_OTUs = str(round(spearman[0], 3)) + "*"
        else:
            spearman_rho_OTUs = str(round(spearman[0], 3))

        # create subplots
        fig = make_subplots(rows=1, cols=2, subplot_titles=("A)", "B)"))
        # calculate line if best fit for reads and add the scatter plot
        df = pd.DataFrame({'X':x1, 'Y':y1})
        df['bestfit'] = sm.OLS(df['Y'],sm.add_constant(df['X'])).fit().fittedvalues
        fig.add_trace(go.Scatter(name="Reads",x=x1, y=y1, mode='markers', marker=dict(color=color1), showlegend=False),row=1, col=1)
        fig.add_trace(go.Scatter(name="rho(reads)=" + spearman_rho_reads, x=x1, y=df['bestfit'], mode='lines', marker=dict(color=color2)),row=1, col=1)
        fig.update_xaxes(title_text = "# reads (rep2)", row=1, col=1)
        fig.update_yaxes(title_text = "# reads (rep1)", row=1, col=1)
        # calculate line if best fit for OTUs and add the scatter plot
        df = pd.DataFrame({'X':x2, 'Y':y2})
        df['bestfit'] = sm.OLS(df['Y'],sm.add_constant(df['X'])).fit().fittedvalues
        fig.add_trace(go.Scatter(name=clustering_unit, x=x2, y=y2, mode='markers', marker=dict(color=color1), showlegend=False),row=1, col=2)
        title1 = "rho(" + clustering_unit + ")="+  spearman_rho_OTUs
        fig.add_trace(go.Scatter(name=title1, x=x2, y=df['bestfit'], mode='lines', marker=dict(color=color2)),row=1, col=2)
        title2 = "# " + clustering_unit + " (rep2)"
        fig.update_xaxes(title_text = title2, row=1, col=2)
        title3 = "# " + clustering_unit + " (rep3)"
        fig.update_yaxes(title_text = title3, row=1, col=2)
        # update the layouts
        fig.update_layout(height=int(height), width=int(width), template=template, font_size=font_size, title_font_size=font_size, showlegend=True)
        if y_zero == True:
            fig.update_yaxes(rangemode="tozero")
        if x_zero == True:
            fig.update_xaxes(rangemode="tozero")

        ## write files
        output_pdf = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_repcorr_OTUs.pdf")
        output_html = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_repcorr_OTUs.html")
        output_text = Path(str(dirName) + "/" + TaXon_table_xlsx.stem + "_repcorr_OTUs.txt")

        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))
        f = open(output_text, "w")
        f.write("Spearman correlation results\n")
        f.write("Reads\n" + "rho = " + str(spearman_rho_reads) + "\np = " + str(spearman_p_reads) + "\n")
        f.write(clustering_unit + "\n" + "rho = " + str(spearman_rho_OTUs) + "\np = " + str(spearman_p_OTUs) + "\n")
        f.close()

        ## ask to show file
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## print closing text
        closing_text = "Plots are found under:\n" + "Projects/Replicate_analysis/"
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        ttt_log("replicate correlation analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)

    else:
        sg.Popup("Correlation analyses can only be performed with two replicates.")
