def per_taxon_analysis(TaXon_table_xlsx, height, width, taxonomic_level, path_to_outdirs, template, theme, font_size):

    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import itertools, webbrowser

    ## collect plot variables
    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]
    height = int(height)
    width = int(width)

    ## load taxon table
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df = TaXon_table_df.fillna("nan")

    ## collect the taxa to test on
    taxa = sorted(list(set([taxon for taxon in TaXon_table_df[taxonomic_level].values.tolist() if taxon != "nan"])))

    ## check if there are more than 8 taxa
    answer = "Yes"
    if len(taxa) > 8:
        answer = sg.PopupYesNo("There are more than 8 taxa detected. This can render the plot difficult to read. Continue anyway?")

    if answer == "Yes":
        ## collect the OTUs
        OTUs = TaXon_table_df["ID"].values.tolist()

        ## count the number of OTUs per taxon
        n_OTUs = [TaXon_table_df[taxonomic_level].values.tolist().count(taxon) for taxon in taxa]

        ## collect all OTUs on species level
        OTU_species = [OTU for OTU in TaXon_table_df[[taxonomic_level, "Species"]].values.tolist() if OTU[1] != "nan"]
        OTU_species.sort()
        OTU_species = list(k for k,_ in itertools.groupby(OTU_species))
        OTU_species = [OTU[0] for OTU in OTU_species]
        n_species = [OTU_species.count(taxon) for taxon in taxa]

        ## count reads for each taxon
        n_reads = []
        for taxon in taxa:
            n_reads.append(sum([sum(OTU[10:]) for OTU in TaXon_table_df[TaXon_table_df[taxonomic_level]==taxon].values.tolist()]))

        ## calculate the read proportions
        reads_sum = sum(n_reads)
        n_reads = [round(reads / reads_sum * 100, 2) for reads in n_reads]

        ## create subplots
        fig = make_subplots(rows=1, cols=2, subplot_titles=("A)", "B)"))

        ## percentage of reads per taxonomic level
        hovertext = 'Taxon: %{x}, Reads: %{y}'
        fig.add_trace(go.Bar(hovertemplate=hovertext, name="",x=taxa, y=n_reads),row=1, col=1)
        fig.update_yaxes(title_text = "reads (%)", title_standoff=5, row=1, col=1)
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value, showlegend=False, row=1, col=1)

        ## Number of OTUs
        hovertext = 'Taxon: %{x}, OTUs: %{y}'
        fig.add_trace(go.Bar(hovertemplate=hovertext, name="",x=taxa, y=n_OTUs, text=n_OTUs, showlegend=False),row=1, col=2)
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value, row=1, col=2)
        fig.update_yaxes(title_text = "# OTUs", title_standoff=5, row=1, col=2, rangemode="tozero")

        ## Number of OTUs on species level
        hovertext = 'Taxon: %{x}, Species: %{text}'
        fig.add_trace(go.Scatter(textposition = "top center", hovertemplate=hovertext, text=n_species, name="Species",x=taxa, y=n_OTUs, showlegend=True, mode='text'),row=1, col=2)
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value, row=1, col=2)

        ## fig.add_annotation( text='─ Species', align='left', showarrow=False, xref='paper', yref='paper', x=1.05, y=0.5, bordercolor='black', borderwidth=1)

        ## update the layout
        fig.update_layout(barmode='stack', height=int(height), width=int(width), template=template, showlegend=True, font_size=font_size, title_font_size=font_size)

        ## write ouput files
        output_pdf = Path(str(path_to_outdirs) + "/Per_taxon_statistics/" + TaXon_table_xlsx.stem + "_" + taxonomic_level + ".pdf")
        output_html = Path(str(path_to_outdirs) + "/Per_taxon_statistics/" + TaXon_table_xlsx.stem + "_" + taxonomic_level + ".html")
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

        ## ask to show file
        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html))

        ## print closing text
        closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Per_taxon_statistics/"
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        ttt_log("per taxon statistics", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)
