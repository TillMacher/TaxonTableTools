import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import itertools, webbrowser
from taxontabletools.taxontable_manipulation import strip_metadata

def per_taxon_analysis(TaXon_table_xlsx, height, width, taxonomic_level, path_to_outdirs, template, theme, font_size, clustering_unit):

    ## save the taxon title
    taxon_title = clustering_unit

    ## collect plot variables
    color1 = theme[0]
    color2 = theme[1]
    opacity_value = theme[2]
    height = int(height)
    width = int(width)

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('nan')
    TaXon_table_df = strip_metadata(TaXon_table_df)

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
        n_OTUs = {taxon:TaXon_table_df[taxonomic_level].values.tolist().count(taxon) for taxon in taxa}
        n_OTUs = dict(sorted(n_OTUs.items(), key=lambda item: item[1], reverse=True))
        taxa = list(n_OTUs.keys())

        ## collect all OTUs on species level
        OTU_species = [OTU for OTU in TaXon_table_df[[taxonomic_level, "Species"]].values.tolist() if OTU[1] != "nan"]
        OTU_species.sort()
        OTU_species = list(k for k,_ in itertools.groupby(OTU_species))
        OTU_species = [OTU[0] for OTU in OTU_species]
        n_species = {taxon:OTU_species.count(taxon) for taxon in taxa}

        ## count reads for each taxon
        n_reads = []
        for taxon in taxa:
            n_reads.append(sum([sum(OTU[10:]) for OTU in TaXon_table_df[TaXon_table_df[taxonomic_level]==taxon].values.tolist()]))
        if (taxonomic_level == "Species" or taxonomic_level == "Genus"):
            x_values = ["<i>" + taxon + "</i>" for taxon in taxa]
        else:
            x_values = taxa

        ## calculate the read proportions
        reads_sum = sum(n_reads)
        n_reads = [round(reads / reads_sum * 100, 2) for reads in n_reads]

        ## create subplots
        fig = go.Figure()
        ## percentage of reads per taxonomic level
        hovertext = 'Taxon: %{x}, Reads: %{y}'
        fig.add_trace(go.Bar(hovertemplate=hovertext, name="",x=x_values, y=n_reads))
        fig.update_yaxes(title_text = "reads (%)", title_standoff=5)
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value, showlegend=False)

        ## update the layout
        fig.update_layout(barmode='stack', uniformtext_mode='hide', height=int(height), width=int(width), template=template, showlegend=False, font_size=font_size, title_font_size=font_size)
        fig.update_xaxes(tickmode='linear')
        fig.update_xaxes(tickangle=-90)

        ## write ouput files
        output_pdf = Path(path_to_outdirs).joinpath('Per_taxon_statistics', '{}_{}_reads.pdf'.format(TaXon_table_xlsx.stem, taxonomic_level))
        output_html1 = Path(path_to_outdirs).joinpath('Per_taxon_statistics', '{}_{}_reads.html'.format(TaXon_table_xlsx.stem, taxonomic_level))
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html1))


        ## create subplots
        fig = go.Figure()
        ## Number of OTUs
        hovertext = 'Taxon: %{x}, OTUs: %{y}'
        title_text = "# " + taxon_title
        fig.add_trace(go.Bar(hovertemplate=hovertext, name="",x=x_values, y=list(n_OTUs.values()), showlegend=False))
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value)
        fig.update_yaxes(title_text=title_text, title_standoff=5, rangemode="tozero")

        ## Number of OTUs on species level
        hovertext = 'Taxon: %{x}, Species: %{text}'
        fig.add_trace(go.Scatter(textposition = "top center", hovertemplate=hovertext, text=list(n_species.values()), name="Species",x=x_values, y=list(n_OTUs.values()), showlegend=False, mode='text'))
        fig.update_traces(marker_color=color1, marker_line_color=color2, marker_line_width=1, opacity=opacity_value)

        ## update the layout
        fig.update_layout(barmode='stack', uniformtext_mode='hide', height=int(height), width=int(width), template=template, showlegend=False, font_size=font_size, title_font_size=font_size)
        fig.update_xaxes(tickmode='linear')
        fig.update_xaxes(tickangle=-90)

        ## write ouput files
        output_pdf = Path(path_to_outdirs).joinpath('Per_taxon_statistics', '{}_{}_OTUs.pdf'.format(TaXon_table_xlsx.stem, taxonomic_level))
        output_html2 = Path(path_to_outdirs).joinpath('Per_taxon_statistics', '{}_{}_OTUs.html'.format(TaXon_table_xlsx.stem, taxonomic_level))
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html2))

        ## ask to show file
        answer = sg.PopupYesNo('Show plots?', keep_on_top=True)
        if answer == "Yes":
            webbrowser.open('file://' + str(output_html1))
            webbrowser.open('file://' + str(output_html2))

        ## print closing text
        closing_text = "Plots are found under: " + str(path_to_outdirs) + "/Per_taxon_statistics/"
        sg.Popup(closing_text, title="Finished", keep_on_top=True)

        ## write to log
        from taxontabletools.create_log import ttt_log
        ttt_log("per taxon statistics", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)
