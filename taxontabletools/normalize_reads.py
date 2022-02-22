import random
import PySimpleGUI as sg
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import math, webbrowser
import plotly.express as px
from statistics import mean
from collections import Counter
from taxontabletools.taxontable_manipulation import strip_metadata
from taxontabletools.taxontable_manipulation import collect_metadata
from taxontabletools.taxontable_manipulation import add_metadata

def collect_reads(TaXon_table_xlsx):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df = strip_metadata(TaXon_table_df)
    samples = TaXon_table_df.columns.tolist()[10:]

    reads_per_sample = [sum(TaXon_table_df[sample].values.tolist()) for sample in samples]
    reads_sum = sum(reads_per_sample)
    reads_max = max(reads_per_sample)
    reads_mean = int(mean(reads_per_sample))
    reads_min = min(reads_per_sample)

    return [reads_sum, reads_max, reads_mean, reads_min]

def normalize_reads(TaXon_table_xlsx, path_to_outdirs, sub_sample_size):

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)
    samples = TaXon_table_df.columns.tolist()[10:]

    ## count the lowest number of reads to use as threshold
    reads_list = []
    for sample in samples:
        reads_list.append(sum(TaXon_table_df[sample].values.tolist()))

    ## convert sub_sample_size to int
    sub_sample_size = int(sub_sample_size)

    ## collect all OTUs
    OTU_list = TaXon_table_df["ID"].values.tolist()

    df_out = pd.DataFrame()
    required_columns = TaXon_table_df.columns.tolist()[0:10]
    df_out[required_columns] = TaXon_table_df[required_columns]

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(samples)
    ############################################################################

    warnings = []

    for sample in samples:
        ## filter sample from data
        read_df = TaXon_table_df[[sample, "ID"]]
        ## drop empty OTUs
        read_df = read_df[read_df[sample] != 0]
        ## check if sample can be normalized, otherwise just keep all reads and OTUs
        if sub_sample_size <= sum(read_df[sample].values.tolist()):
            ## create read list to draw the subsamples from
            read_list = pd.Series(np.repeat(read_df['ID'].to_list(), read_df[sample].to_list()))
            ## create a random subsample
            sub_sample = read_list.sample(n = sub_sample_size)
            ## count the number of reads per OTU
            sub_sample_reads = dict(pd.Series(sub_sample).value_counts())
            ## create a sorted OTU list
            OTU_sample_list = []
            for OTU in OTU_list:
                if OTU in sub_sample_reads.keys():
                    OTU_sample_list.append(sub_sample_reads[OTU])
                else:
                    OTU_sample_list.append(0)
        else:
            OTU_sample_list = TaXon_table_df[[sample]]
            warnings.append(sample)

        ## add OTUs to dataframe
        df_out[sample] = OTU_sample_list

        ############################################################################
        event, values = window_progress_bar.read(timeout=1)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    ## add already existing metadata back to the df
    if len(TaXon_table_df_metadata.columns) != 1:
        df_out = add_metadata(df_out, TaXon_table_df_metadata)

    ## write output file
    output_pdf = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem + "_normalized.xlsx")
    df_out.to_excel(output_pdf, index=False)

    ## also count stats if selected
    OTU_loss = []
    for sample in samples:
        original_OTUs = len([i[1] for i in TaXon_table_df[[sample, "ID"]].values.tolist() if i[0] != 0])
        subsample_OTUs = len([i[1] for i in df_out[[sample, "ID"]].values.tolist() if i[0] != 0])

        if subsample_OTUs == 0:
            d_OTUs = 100
        else:
            d_OTUs = round(100 - subsample_OTUs / original_OTUs * 100, 4)
        OTU_loss.append(d_OTUs)

    closing_text = "All samples were normalized to maximum {} reads!\n\nMean discarded OTUs: {} %.\n\nSamples below the normalization theshold: {}".format(str(sub_sample_size), str(round(mean(OTU_loss), 4)), len(warnings))
    sg.Popup(closing_text, title='Finished')

    ## write to log
    from taxontabletools.create_log import ttt_log
    ttt_log("normalization", "processing", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)









            #
