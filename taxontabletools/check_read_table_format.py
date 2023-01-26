import PySimpleGUI as sg
import pandas as pd
import numpy as np
from pathlib import Path

def check_read_table_format_TTT(read_table_xlsx):

    try:
        read_table_df = pd.read_excel(Path(read_table_xlsx))
        read_table_df = read_table_df.replace(np.nan, 'nan', regex=True)
    except:
        sg.PopupError("Read table must be in excel sheet format (.xlsx)", keep_on_top=True)
        raise RuntimeError("Read table must be in excel sheet format (.xlsx)")

    ###################################
    # A) header prompt
    read_table_df_header = read_table_df.columns.tolist()
    header_prompt_OTUs = read_table_df_header[0]
    header_prompt_sequences = read_table_df_header[-1]

    if header_prompt_OTUs != "ID":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_OTUs + "\n" + "\n" + "Prompt: ID"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise RuntimeError(ErrorMessage)

    if header_prompt_sequences != "Seq":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_sequences + "\n" + "\n" + "Prompt: Sequence or Seq"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise RuntimeError(ErrorMessage)

    ###################################
    # B) OTUs prompt

    OTU_list = read_table_df['ID'].values.tolist()

    # loop through all available OTUs
    for OTU in OTU_list:

        # define a variable Error message
        ErrorMessage = "Oops! Something is wrong with the OTU ID: " + OTU

        try:
            # split the OTU, which should be containing of two elements
            p1 = OTU.split("_")[0]
            p2 = OTU.split("_")[1]

        except:
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

        # check if the second part is an integer
        # convert p2 to an integer
        try:
            p2 = int(p2)
        except :
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # D) Samples prompt

    samples_list = read_table_df.columns.tolist()[10:]

    for sample in samples_list:
        if " " in sample:
            ErrorMessage = "Please do not use spaces in the sample names:\n" + sample
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # E) Reads prompt

    OTU_reads_list = read_table_df.values.tolist()

    for OTU_reads in OTU_reads_list:
        read_numbers = OTU_reads[1:-1]
        for read_number in read_numbers:
            try:
                read_number = int(read_number)
            except:
                OTU = OTU_reads[0]
                ErrorMessage = "Please check your read numbers in " + OTU + " -> " + str(read_number)
                sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
                raise RuntimeError(ErrorMessage)

    ###################################
    # D) Samples prompt

    samples_list = read_table_df.columns.tolist()[1:-1]

    for sample in samples_list:
        if " " in sample:
            ErrorMessage = "Please do not use spaces in the sample names:\n" + sample
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # Wrap up
    sg.Popup("Your file looks great and is ready to use!", title="Read table check", keep_on_top=True)

def check_read_table_format_qiime2(read_table_xlsx):

    try:
        read_table_df = pd.read_csv(Path(read_table_xlsx), sep="\t")
    except:
        sg.PopupError("Read table must be in tsv format (.tsv)", keep_on_top=True)
        raise RuntimeError("Read table must be in tsv format (.tsv)")

    ###################################
    # A) header prompt
    read_table_df_header = read_table_df.columns.tolist()
    header_prompt_OTUs = read_table_df_header[0]
    header_prompt_sequences = read_table_df_header[-1]

    if header_prompt_OTUs != "id":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_OTUs + "\n" + "\n" + "Prompt: id"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise RuntimeError(ErrorMessage)

    if header_prompt_sequences != "Seq":
        ErrorMessage = "Oops! Something is wrong with the header: " + header_prompt_sequences + "\n" + "\n" + "Prompt: Sequence or Seq"
        sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
        raise RuntimeError(ErrorMessage)

    ###################################
    # B) OTUs prompt

    ## removed in v 1.2.0

    ###################################
    # D) Samples prompt

    samples_list = read_table_df.columns.tolist()[1:-1]

    for sample in samples_list:
        if " " in sample:
            ErrorMessage = "Please do not use spaces in the sample names:\n" + sample
            sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
            raise RuntimeError(ErrorMessage)

    ###################################
    # E) Reads prompt

    OTU_reads_list = read_table_df.values.tolist()

    for OTU_reads in OTU_reads_list:
        if OTU_reads[0] != "#q2:types":
            read_numbers = OTU_reads[1:-1]
            for read_number in read_numbers:
                try:
                    read_number = int(read_number)
                except:
                    OTU = OTU_reads[0]
                    ErrorMessage = "Please check your read numbers in " + OTU + " -> " + str(read_number)
                    sg.PopupError(ErrorMessage, title="Error", keep_on_top=True)
                    raise RuntimeError(ErrorMessage)

    ###################################
    # Wrap up
    sg.Popup("Your file looks great and is ready to use!", title="Read table check", keep_on_top=True)
