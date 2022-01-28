import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path

def replicate_consistency_filter(TaXon_table_xlsx, suffix_list, path_to_outdirs, consistency):

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    TaXon_table_df_orr = TaXon_table_df

    sample_names = TaXon_table_df.columns[10:].tolist()
    OTUs = TaXon_table_df["ID"].values.tolist()

    derep_sample_names_dict =  {}
    unique_sample_names_list = []
    replicates_dict = {}

    for sample in sample_names:
        sample_name = sample.split("_")[0:-1]
        unique_sample_names_list.append("_".join(sample_name))

    unique_sample_names_set = sorted(set(unique_sample_names_list))

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(unique_sample_names_set) + 1
    ############################################################################

    ## merge and replicate consistency version
    if consistency == True:
        no_replicates_list = []
        for sample in unique_sample_names_set:

            for i, suffix in enumerate(suffix_list):
                replicates_dict["rep_" + str(i)] = sample + "_" + str(suffix_list[i])

            replicate_names_list = list(replicates_dict.values())

            try:
                new_df = TaXon_table_df[replicate_names_list]
                header = new_df.columns.tolist()
                processed_reads = []

                for n_reads in new_df.values.tolist():
                    if 0 in n_reads:
                        if len(set(n_reads)) > 1:
                            n_reads = len(n_reads) * [0]
                    processed_reads.append(n_reads)

                df_out = pd.DataFrame(processed_reads)
                df_out.columns = header
                TaXon_table_df = TaXon_table_df.drop(replicate_names_list, axis=1)
                TaXon_table_df[sample] = df_out.sum(axis=1)

            except:
                no_replicates_list.append(sample)

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                print('Cancel')
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        if len(no_replicates_list) == len(unique_sample_names_set):
            sg.PopupError("No replicates found. Please check your replicate suffixes.")
        else:
            dropped_OTUs_list = []
            # filter for 0 hit OTUs (can happen after consistency filtering)
            columns = TaXon_table_df.columns.tolist()
            TaXon_table_list = TaXon_table_df.values.tolist()
            TaXon_table_list_final = []
            for entry in TaXon_table_list:
                if sum(entry[10:]) != 0:
                    TaXon_table_list_final.append(entry)
                else:
                    dropped_OTUs_list.append(entry[0])

            ## Create a Pandas Excel writer using XlsxWriter as the engine.
            taxon_tables_directory = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem)
            output_xlsx = Path(str(taxon_tables_directory) + "_cons.xlsx")
            writer = pd.ExcelWriter(output_xlsx, engine='xlsxwriter')

            ## create a dataframe and write the main table
            TaXon_table_df = pd.DataFrame(TaXon_table_list_final, columns=columns)
            TaXon_table_df.to_excel(writer, sheet_name='TaXon table', index=False)

            answer = sg.PopupYesNo("Write dropped OTUs to a seperate sheet?")
            if answer == "Yes":
                dropped_OTUs_df = pd.DataFrame([i for i in TaXon_table_df_orr.values.tolist() if i[0] in dropped_OTUs_list], columns=TaXon_table_df_orr.columns.tolist())
                dropped_OTUs_df.to_excel(writer, sheet_name='Dropped OTUs', index=False)

            original_reads = sum([sum(i[10:]) for i in TaXon_table_df_orr.values.tolist()])
            remaining_reads = original_reads - sum([sum(i[10:]) for i in TaXon_table_df.values.tolist()])
            rel_remaining_reads = 100 - round(remaining_reads / original_reads * 100, 2)


            closing_text = "Taxon table is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
            reads_info = str(rel_remaining_reads) + "% of reads were kept."
            OTU_info = str(len(dropped_OTUs_list)) + " OTUs were removed from the dataset.\n"

            sg.Popup(closing_text + "\n\n" + reads_info + "\n" + OTU_info, title="Finished", keep_on_top=True)

            writer.save()

            from taxontabletools.create_log import ttt_log
            ttt_log("replicate consistency", "processing", TaXon_table_xlsx.name, output_xlsx.name, "consistency merged", path_to_outdirs)

    ## merge only version
    else:
        no_replicates_list = []
        for sample in unique_sample_names_set:

            for i, suffix in enumerate(suffix_list):
                replicates_dict["rep_" + str(i)] = sample + "_" + str(suffix_list[i])

            replicate_names_list = list(replicates_dict.values())

            try:
                new_df = TaXon_table_df[replicate_names_list]
                TaXon_table_df = TaXon_table_df.drop(replicate_names_list, axis=1)
                TaXon_table_df[sample] = new_df.sum(axis=1)
            except:
                no_replicates_list.append(sample)

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                print('Cancel')
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        if len(no_replicates_list) == len(unique_sample_names_set):
            sg.PopupError("No replicates found. Please check your replicate suffixes.")

        else:
            taxon_tables_directory = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem)
            output_xlsx = Path(str(taxon_tables_directory) + "_merged.xlsx")

            TaXon_table_df.to_excel(output_xlsx, sheet_name='TaXon table', index=False)

            closing_text = "Taxon table is found under:\n" + '/'.join(str(output_xlsx).split("/")[-4:])
            sg.Popup(closing_text, title="Finished", keep_on_top=True)

            from taxontabletools.create_log import ttt_log
            ttt_log("replicate merging", "processing", TaXon_table_xlsx.name, output_xlsx.name, "merged", path_to_outdirs)
